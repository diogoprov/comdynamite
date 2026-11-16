# Integrating over posterior configurations in the DAMOCLES paper (dev/plano-artigo-damocles.md,
# S3): how to summarise fits over samples of regimes or branch rates from a method such as BAMM
# or ClaDS. Synthetic posteriors with controlled uncertainty (validation with real ClaDS output
# to come).
# Trees: sim_shift_tree(112, 1, 2) redrawn until the shifted regime holds 0.25 of the tips (factor
# 0.8-1.2), height one, mu = 1, n = 14 or 40 communities, K configurations per replicate.
# - "regimes": gamma = 5 in the background and 5 x ratio (1: none; 1.5) below the true shift.
#   Each posterior sample puts the one shift on the true edge with probability p_right (0.7 or
#   0.4), otherwise on a neighbouring edge (parent, child or sister) drawn at random; the regime
#   is "1" on the shift edge and below it, as in sim_shift_tree(). The rate of colonisation varies
#   (vary = "gamma"); the effect is the ratio gamma_1 / gamma_0.
# - "covariate": a branch rate x with log x Brownian along the tree (unit increments per edge,
#   rescaled to standard deviation 0.5 across edges), gamma = 5 (x / geometric mean of x)^beta,
#   beta = 0 (none) or 1. Posterior: a centre m = log x + N(0, s) per edge and K samples
#   m + N(0, s), with s = 0.3 or 0.7 (0.7 is the mean posterior SD of the log branch rates of
#   the ClaDS example of data-raw/damocles_clads.R, against an SD of 0.56 across edges).
# Estimators:
# - oracle: the true configuration;
# - (a) per sample: damocles_fit() on each configuration; estimate the median across samples,
#   rejection when the median p-value is below 0.05 (half the samples reject);
# - (b) one summary configuration: the regime of each edge by majority over the samples, or the
#   posterior mean of log x;
# - (c) integrated likelihood: L(theta) = mean over k of L(theta | configuration k), maximised
#   once; LRT against constant rates (which do not depend on the configuration), 1 df;
# - (d) SIMEX on (b), covariate only: noise of variance lambda v added to the posterior mean
#   (v: posterior variance of each edge), lambda = 1 and 2, three draws each, quadratic
#   extrapolation to lambda = -1 (estimate only).
# REPS replicates per setting (default 100), K (default 20), WORKERS (default 16), JOBS (default
# data-raw/s3_jobs, ignored by git), OUT (default data-raw/s3_integration.rds).
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE))
workers <- as.integer(Sys.getenv("WORKERS", "16"))
reps <- as.integer(Sys.getenv("REPS", "100"))
K <- as.integer(Sys.getenv("K", "20"))
jobs_dir <- Sys.getenv("JOBS", "data-raw/s3_jobs")
out <- Sys.getenv("OUT", "data-raw/s3_integration.rds")
RNGkind("L'Ecuyer-CMRG")
dir.create(jobs_dir, showWarnings = FALSE, recursive = TRUE)

quiet <- function(x) suppressWarnings(x)
focal_tree <- function(n_tips, fraction, max_tries = 1000) {
  for (i in seq_len(max_tries)) {
    sim <- sim_shift_tree(n_tips, 1, 2)
    f <- mean(sim$tip_regime > 0)
    if (f >= 0.8 * fraction && f <= 1.2 * fraction) return(sim)
  }
  stop("no tree with the focal fraction in ", max_tries, " tries", call. = FALSE)
}
# edges below (and including) edge e
clade_edges <- function(tree, e) {
  below <- e
  frontier <- tree$edge[e, 2]
  while (length(frontier)) {
    kids <- which(tree$edge[, 1] %in% frontier)
    below <- c(below, kids)
    frontier <- tree$edge[kids, 2]
  }
  below
}
shift_regimes <- function(tree, e) {
  r <- rep("0", nrow(tree$edge))
  r[clade_edges(tree, e)] <- "1"
  factor(r, levels = c("0", "1"))
}
neighbours <- function(tree, e) {
  up <- tree$edge[e, 1]
  c(which(tree$edge[, 2] == up), which(tree$edge[, 1] == tree$edge[e, 2]),
    setdiff(which(tree$edge[, 1] == up), e))
}
log_rate_edges <- function(tree) {
  ord <- ape::reorder.phylo(tree, "cladewise", index.only = TRUE)
  z <- numeric(nrow(tree$edge))
  for (e in ord) {
    parent_edge <- which(tree$edge[, 2] == tree$edge[e, 1])
    z[e] <- (if (length(parent_edge)) z[parent_edge] else 0) + stats::rnorm(1)
  }
  0.5 * (z - mean(z)) / stats::sd(z)
}

fit_one <- function(comm, tree, start, regimes = NULL, covariate = NULL) {
  quiet(damocles_fit(comm, tree, shared = TRUE, regimes = regimes, edge_covariate = covariate,
                     vary = "gamma", start = start))
}
effect_of <- function(fit, kind) {
  if (kind == "regimes") fit$gamma_1[1] / fit$gamma_0[1] else fit$beta_gamma[1]
}
integrated_fit <- function(comm, tree, configs, kind, start) {
  presence <- damocles_presence(comm, tree)
  post <- damocles_postorder(tree)
  design_of <- function(cfg) if (kind == "regimes") {
    damocles_design(tree, cfg, NULL, "gamma", 0)
  } else {
    damocles_design(tree, NULL, cfg, "gamma", 0)
  }
  objectives <- lapply(configs, function(cfg) {
    damocles_objective(presence, tree, post, "both", design_of(cfg))
  })
  design <- design_of(configs[[1]])
  starts <- damocles_par_starts(damocles_starts(post$tree, start), design)
  f <- function(par) {
    ll <- -vapply(objectives, function(o) o(par), numeric(1))
    ll[ll <= -1e10] <- -Inf
    if (all(!is.finite(ll))) return(1e10)
    top <- max(ll)
    -(top + log(mean(exp(ll - top))))
  }
  best <- damocles_optim(f, starts)
  values <- damocles_par_values(best$par, design)
  list(loglik = -best$value,
       effect = if (kind == "regimes") values[["gamma_1"]] / values[["gamma_0"]] else
         values[["beta_gamma"]], convergence = best$convergence)
}

one <- function(s) {
  sim <- focal_tree(112, 0.25)
  tree <- sim$tree
  n_edge <- nrow(tree$edge)
  if (s$kind == "regimes") {
    truth <- sim$regimes
    gamma <- ifelse(truth != "0", 5 * s$effect, 5)
    true_edge <- which(truth != "0" & !(tree$edge[, 1] %in% tree$edge[truth != "0", 2]))
    near <- neighbours(tree, true_edge)
    configs <- lapply(seq_len(K), function(k) {
      e <- if (stats::runif(1) < s$noise) true_edge else near[sample.int(length(near), 1)]
      shift_regimes(tree, e)
    })
    share <- rowMeans(vapply(configs, function(r) r == "1", logical(n_edge)))
    summary_cfg <- factor(ifelse(share > 0.5, "1", "0"), levels = c("0", "1"))
    if (nlevels(droplevels(summary_cfg)) < 2) summary_cfg <- configs[[1]]
    arg <- function(cfg) list(regimes = cfg)
    oracle_cfg <- factor(ifelse(truth != "0", "1", "0"), levels = c("0", "1"))
  } else {
    log_x <- log_rate_edges(tree)
    gamma <- 5 * exp(s$effect * (log_x - mean(log_x)))
    centre <- log_x + stats::rnorm(n_edge, 0, s$noise)
    draws <- matrix(centre + stats::rnorm(n_edge * K, 0, s$noise), n_edge, K)
    configs <- lapply(seq_len(K), function(k) exp(draws[, k]))
    summary_cfg <- exp(rowMeans(draws))
    arg <- function(cfg) list(covariate = cfg)
    oracle_cfg <- exp(log_x)
  }
  comm <- simulate_damocles(tree, gamma, 1, n = s$n_comm)
  null <- quiet(damocles_fit(comm, tree, shared = TRUE))
  start <- c(null$gamma[1], null$mu[1])
  ll0 <- attr(null, "loglik_total")
  p_of <- function(ll) stats::pchisq(pmax(2 * (ll - ll0), 0), 1, lower.tail = FALSE)
  fit_cfg <- function(cfg) {
    f <- do.call(fit_one, c(list(comm = comm, tree = tree, start = start), arg(cfg)))
    c(effect = effect_of(f, s$kind), p = p_of(attr(f, "loglik_total")), conv = f$convergence[1])
  }
  oracle <- fit_cfg(oracle_cfg)
  per <- vapply(configs, fit_cfg, numeric(3))
  summ <- fit_cfg(summary_cfg)
  int <- integrated_fit(comm, tree, configs, s$kind, start)
  simex <- NA_real_
  if (s$kind == "covariate") {
    v <- apply(draws, 1, stats::var)
    b <- vapply(c(1, 2), function(lambda) stats::median(vapply(1:3, function(i) {
      fit_cfg(exp(rowMeans(draws) + stats::rnorm(n_edge, 0, sqrt(lambda * v))))[["effect"]]
    }, numeric(1))), numeric(1))
    simex <- 3 * summ[["effect"]] - 3 * b[1] + b[2]
  }
  data.frame(oracle_effect = oracle[["effect"]], oracle_p = oracle[["p"]],
             per_effect = stats::median(per["effect", ]), per_p = stats::median(per["p", ]),
             per_share = mean(per["p", ] < 0.05),
             summary_effect = summ[["effect"]], summary_p = summ[["p"]],
             int_effect = int$effect, int_p = p_of(int$loglik), simex_effect = simex,
             convergence = max(oracle[["conv"]], per["conv", ], summ[["conv"]], int$convergence))
}

settings <- rbind(
  expand.grid(kind = "regimes", n_comm = c(14, 40), noise = c(0.7, 0.4), effect = c(1, 1.5),
              stringsAsFactors = FALSE),
  expand.grid(kind = "covariate", n_comm = c(14, 40), noise = c(0.3, 0.7), effect = c(0, 1),
              stringsAsFactors = FALSE))
only <- Sys.getenv("SETTINGS", "")
run <- seq_len(nrow(settings))
if (nzchar(only)) run <- intersect(run, as.integer(strsplit(only, ",")[[1]]))
jobs <- expand.grid(setting = run, rep = seq_len(reps))
jobs$file <- file.path(jobs_dir, sprintf("i%02d_r%04d.rds", jobs$setting, jobs$rep))
todo <- which(!file.exists(jobs$file))
todo <- todo[order(-settings$n_comm[jobs$setting[todo]])]
cat("S3 integration:", length(todo), "of", nrow(jobs), "replicates on", workers, "workers, K =",
    K, "\n")
started <- Sys.time()
invisible(parallel::mclapply(todo, function(j) {
  s <- settings[jobs$setting[j], ]
  set.seed(303 + 10000 * jobs$setting[j] + jobs$rep[j])
  t0 <- Sys.time()
  r <- tryCatch(one(s), error = function(e) data.frame(error = conditionMessage(e)))
  if (is.null(r$error)) r$error <- NA_character_
  r$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(cbind(s, rep = jobs$rep[j], r, row.names = NULL), jobs$file[j])
  NULL
}, mc.cores = workers, mc.preschedule = FALSE))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
cols <- c("oracle_effect", "oracle_p", "per_effect", "per_p", "per_share", "summary_effect",
          "summary_p", "int_effect", "int_p", "simex_effect", "convergence", "error", "seconds")
files <- jobs$file[file.exists(jobs$file)]
results <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f)
  for (k in setdiff(cols, names(x))) x[[k]] <- NA
  x[, c(names(settings), "rep", cols)]
}))

# Summary: median estimate and RMSE of the effect (the ratio on the log scale for regimes, beta
# for the covariate), and rejection at 5%.
keys <- c("kind", "n_comm", "noise", "effect")
rate <- function(p) mean(p < 0.05, na.rm = TRUE)
summary <- do.call(rbind, lapply(split(results, results[, keys], drop = TRUE), function(d) {
  truth <- if (d$kind[1] == "regimes") log(d$effect[1]) else d$effect[1]
  sc <- if (d$kind[1] == "regimes") log else identity
  est <- function(x) c(stats::median(x, na.rm = TRUE),
                       sqrt(mean((sc(x) - truth)^2, na.rm = TRUE)))
  e <- rbind(oracle = est(d$oracle_effect), per = est(d$per_effect),
             summary = est(d$summary_effect), int = est(d$int_effect),
             simex = est(d$simex_effect))
  data.frame(d[1, keys],
             med_oracle = e["oracle", 1], med_per = e["per", 1], med_summary = e["summary", 1],
             med_int = e["int", 1], med_simex = e["simex", 1],
             rmse_oracle = e["oracle", 2], rmse_per = e["per", 2],
             rmse_summary = e["summary", 2], rmse_int = e["int", 2], rmse_simex = e["simex", 2],
             rej_oracle = rate(d$oracle_p), rej_per = rate(d$per_p),
             rej_summary = rate(d$summary_p), rej_int = rate(d$int_p),
             errors = sum(!is.na(d$error)), not_converged = sum(d$convergence != 0, na.rm = TRUE))
}))
summary <- summary[order(summary$kind, summary$effect, summary$noise, summary$n_comm), ]
cat("S3 integration:", nrow(results), "replicates in", round(elapsed, 1), "min\n")
options(width = 250)
print(summary, row.names = FALSE, digits = 3)
saveRDS(list(settings = settings, results = results, summary = summary, K = K,
             minutes = elapsed), out)
