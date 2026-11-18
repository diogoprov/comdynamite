# Validation of the S3 posterior-integration results (dev/plano-artigo-damocles.md, S3) with
# real ClaDS posteriors instead of synthetic ones. RPANDA 2.5 (GPL, only used here; see
# data-raw/damocles_clads.R for its installation).
# Two stages, chosen by MODE:
# - MODE=mcmc: TREES trees (default 16) from RPANDA::sim_ClaDS() (lambda_0 = 0.1, turnover 0.5,
#   sigma = 0.2, alpha = 0.95, 113 extant tips; the SD of the log branch rates across edges is
#   about 0.5, as in the synthetic posteriors and the Caprimulgidae example), each fitted with
#   ClaDS2 (three chains, sample_fraction = 1, thin 20) in blocks of BLOCK iterations (default
#   2500), checkpointed in JOBS (default data-raw/s3c_jobs, ignored by git), until the 95th
#   percentile of the R-hat of the parameters (second half of the chains) is below 1.05 or MAX_IT
#   iterations (default 30000). Rerunning resumes from the checkpoints.
# - MODE=communities: on each tree (rescaled to height one) with the true branch rates x of the
#   simulation, gamma = 5 (x / geometric mean of x)^beta, mu = 1, beta = 0 or 1, n = 14 or 40
#   communities, REPS replicates per tree and setting (default 10). Estimators as in
#   data-raw/s3_integration.R with K = 20 posterior samples of the branch rates (second half of
#   the chains, transformed as in RPANDA::getMAPS_ClaDS()): oracle (true rates), per sample,
#   posterior mean of the log rates, the MAP rates of getMAPS_ClaDS(), integrated likelihood and
#   SIMEX on the posterior mean. Also per tree: R-hat, correlation and slope of the true log
#   rates on their posterior means, and coverage of the 95% intervals.
# WORKERS (default 16), OUT (default data-raw/s3_clads_validation.rds).
mode <- Sys.getenv("MODE", "mcmc")
src <- readLines("data-raw/s3_integration.R")
eval(parse(text = src[seq_len(grep("^one <- function", src) - 1)]))
suppressMessages(library(RPANDA))
workers <- as.integer(Sys.getenv("WORKERS", "16"))
n_trees <- as.integer(Sys.getenv("TREES", "16"))
block <- as.integer(Sys.getenv("BLOCK", "2500"))
max_it <- as.integer(Sys.getenv("MAX_IT", "30000"))
reps <- as.integer(Sys.getenv("REPS", "10"))
K <- as.integer(Sys.getenv("K", "20"))
jobs_dir <- Sys.getenv("JOBS", "data-raw/s3c_jobs")
out <- Sys.getenv("OUT", "data-raw/s3_clads_validation.rds")
dir.create(jobs_dir, showWarnings = FALSE, recursive = TRUE)
tree_file <- function(i) file.path(jobs_dir, sprintf("tree_%02d.rds", i))
mcmc_file <- function(i) file.path(jobs_dir, sprintf("mcmc_%02d.rds", i))

# log branch rates of the kept draws (second half), parameters and R-hat
clads_draws <- function(sampler) {
  n_par <- ncol(sampler$chains[[1]]) - 2
  n <- nrow(sampler$chains[[1]])
  kept <- seq(floor(n / 2), n)
  mc <- coda::mcmc.list(lapply(sampler$chains, function(ch) coda::mcmc(ch[kept, 1:n_par])))
  rhat <- coda::gelman.diag(mc, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
  draws <- do.call(rbind, lapply(sampler$chains, function(ch) ch[kept, , drop = FALSE]))
  log_rates <- t(apply(draws, 1, function(d) d[4] + d[2] * sampler$alpha_effect + d[5:n_par]))
  list(log_rates = log_rates, rhat = rhat)
}

if (mode == "mcmc") {
  for (i in seq_len(n_trees)) if (!file.exists(tree_file(i))) {
    set.seed(3030 + i)
    repeat {
      s <- sim_ClaDS(0.1, 0.5, sigma_lamb = 0.2, alpha_lamb = 0.95, condition = "taxa",
                     taxa_stop = 113, prune_extinct = TRUE)
      if (!is.null(s$tree) && ape::Ntip(s$tree) >= 100) break
    }
    saveRDS(list(tree = s$tree, log_rate = log(s$lamb[s$rates])), tree_file(i))
  }
  cat("S3 ClaDS MCMC:", n_trees, "trees, blocks of", block, "up to", max_it, "iterations\n")
  invisible(parallel::mclapply(seq_len(n_trees), function(i) {
    set.seed(30300 + i)
    tree <- readRDS(tree_file(i))$tree
    state <- if (file.exists(mcmc_file(i))) readRDS(mcmc_file(i)) else
      list(sampler = NULL, iterations = 0, minutes = 0, rhat95 = Inf)
    while (state$rhat95 >= 1.05 && state$iterations < max_it) {
      t0 <- Sys.time()
      state$sampler <- suppressWarnings(fit_ClaDS(tree, sample_fraction = 1, iterations = block,
                                                  thin = 20, it_save = block, nCPU = 1,
                                                  model_id = "ClaDS2", verbose = FALSE,
                                                  mcmcSampler = state$sampler))
      state$iterations <- state$iterations + block
      state$minutes <- state$minutes + as.numeric(difftime(Sys.time(), t0, units = "mins"))
      state$rhat95 <- stats::quantile(clads_draws(state$sampler)$rhat, 0.95)
      saveRDS(state, mcmc_file(i))
      message(sprintf("tree %02d: %d iterations, %.0f min, R-hat 95%% %.3f", i,
                      state$iterations, state$minutes, state$rhat95))
    }
    NULL
  }, mc.cores = min(workers, n_trees), mc.preschedule = FALSE))
  quit(save = "no")
}

# MODE=communities ------------------------------------------------------------------------------
posteriors <- lapply(seq_len(n_trees), function(i) {
  sim <- readRDS(tree_file(i))
  state <- readRDS(mcmc_file(i))
  d <- clads_draws(state$sampler)
  lr <- d$log_rates
  pm <- colMeans(lr)
  lo <- apply(lr, 2, stats::quantile, 0.025)
  hi <- apply(lr, 2, stats::quantile, 0.975)
  tree <- sim$tree
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  list(tree = tree, log_rate = sim$log_rate, draws = lr,
       map = log(getMAPS_ClaDS(state$sampler, thin = 1)[-(1:4)]),
       diagnostics = data.frame(tree = i, iterations = state$iterations,
                                minutes = state$minutes, rhat95 = unname(state$rhat95),
                                rhat_max = max(d$rhat), cor = stats::cor(pm, sim$log_rate),
                                slope = unname(stats::coef(stats::lm(sim$log_rate ~ pm))[2]),
                                sd_true = stats::sd(sim$log_rate), sd_mean = stats::sd(pm),
                                post_sd = mean(apply(lr, 2, stats::sd)),
                                coverage = mean(sim$log_rate >= lo & sim$log_rate <= hi)))
})
diagnostics <- do.call(rbind, lapply(posteriors, `[[`, "diagnostics"))

one_clads <- function(s) {
  p <- posteriors[[s$tree]]
  tree <- p$tree
  n_edge <- nrow(tree$edge)
  gamma <- 5 * exp(s$effect * (p$log_rate - mean(p$log_rate)))
  draws <- t(p$draws[sort(sample.int(nrow(p$draws), K)), , drop = FALSE])
  configs <- lapply(seq_len(K), function(k) exp(draws[, k]))
  comm <- simulate_damocles(tree, gamma, 1, n = s$n_comm)
  null <- quiet(damocles_fit(comm, tree, shared = TRUE))
  start <- c(null$gamma[1], null$mu[1])
  ll0 <- attr(null, "loglik_total")
  p_of <- function(ll) stats::pchisq(pmax(2 * (ll - ll0), 0), 1, lower.tail = FALSE)
  fit_cfg <- function(cfg) {
    f <- fit_one(comm, tree, start, covariate = cfg)
    c(effect = f$beta_gamma[1], p = p_of(attr(f, "loglik_total")), conv = f$convergence[1])
  }
  oracle <- fit_cfg(exp(p$log_rate))
  per <- vapply(configs, fit_cfg, numeric(3))
  mean_cfg <- exp(colMeans(p$draws))
  summ <- fit_cfg(mean_cfg)
  map <- fit_cfg(exp(p$map))
  int <- integrated_fit(comm, tree, configs, "covariate", start)
  v <- apply(p$draws, 2, stats::var)
  b <- vapply(c(1, 2), function(lambda) stats::median(vapply(1:3, function(i) {
    fit_cfg(exp(colMeans(p$draws) + stats::rnorm(n_edge, 0, sqrt(lambda * v))))[["effect"]]
  }, numeric(1))), numeric(1))
  data.frame(oracle_effect = oracle[["effect"]], oracle_p = oracle[["p"]],
             per_effect = stats::median(per["effect", ]), per_p = stats::median(per["p", ]),
             summary_effect = summ[["effect"]], summary_p = summ[["p"]],
             map_effect = map[["effect"]], map_p = map[["p"]],
             int_effect = int$effect, int_p = p_of(int$loglik),
             simex_effect = 3 * summ[["effect"]] - 3 * b[1] + b[2],
             convergence = max(oracle[["conv"]], per["conv", ], summ[["conv"]], map[["conv"]],
                               int$convergence))
}

settings <- expand.grid(tree = seq_len(n_trees), n_comm = c(14, 40), effect = c(0, 1))
jobs <- expand.grid(setting = seq_len(nrow(settings)), rep = seq_len(reps))
jobs$file <- file.path(jobs_dir, sprintf("v%03d_r%03d.rds", jobs$setting, jobs$rep))
todo <- which(!file.exists(jobs$file))
cat("S3 ClaDS communities:", length(todo), "of", nrow(jobs), "replicates on", workers,
    "workers, K =", K, "\n")
started <- Sys.time()
invisible(parallel::mclapply(todo, function(j) {
  s <- settings[jobs$setting[j], ]
  set.seed(3131 + 10000 * jobs$setting[j] + jobs$rep[j])
  t0 <- Sys.time()
  r <- tryCatch(one_clads(s), error = function(e) data.frame(error = conditionMessage(e)))
  if (is.null(r$error)) r$error <- NA_character_
  r$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(cbind(s, rep = jobs$rep[j], r, row.names = NULL), jobs$file[j])
  NULL
}, mc.cores = workers, mc.preschedule = FALSE))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
cols <- c("oracle_effect", "oracle_p", "per_effect", "per_p", "summary_effect", "summary_p",
          "map_effect", "map_p", "int_effect", "int_p", "simex_effect", "convergence", "error",
          "seconds")
results <- do.call(rbind, lapply(jobs$file[file.exists(jobs$file)], function(f) {
  x <- readRDS(f)
  for (k in setdiff(cols, names(x))) x[[k]] <- NA
  x[, c(names(settings), "rep", cols)]
}))
rate <- function(p) mean(p < 0.05, na.rm = TRUE)
rmse <- function(x, truth) sqrt(mean((x - truth)^2, na.rm = TRUE))
summary <- do.call(rbind, lapply(split(results, results[, c("n_comm", "effect")]), function(d) {
  est <- c("oracle", "per", "summary", "map", "int", "simex")
  med <- vapply(est, function(e) stats::median(d[[paste0(e, "_effect")]], na.rm = TRUE), 1)
  err <- vapply(est, function(e) rmse(d[[paste0(e, "_effect")]], d$effect[1]), 1)
  rej <- vapply(setdiff(est, "simex"), function(e) rate(d[[paste0(e, "_p")]]), 1)
  data.frame(n_comm = d$n_comm[1], effect = d$effect[1], t(setNames(med, paste0("med_", est))),
             t(setNames(err, paste0("rmse_", est))), t(setNames(rej, paste0("rej_", names(rej)))),
             errors = sum(!is.na(d$error)), not_converged = sum(d$convergence != 0, na.rm = TRUE))
}))
cat("S3 ClaDS communities:", nrow(results), "replicates in", round(elapsed, 1), "min\n")
options(width = 200)
print(diagnostics, row.names = FALSE, digits = 3)
print(summary, row.names = FALSE, digits = 3)
saveRDS(list(settings = settings, diagnostics = diagnostics, results = results,
             summary = summary, minutes = elapsed), out)
