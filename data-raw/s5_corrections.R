# Corrections for spatial dependence in the DAMOCLES paper (dev/plano-artigo-damocles.md, S5,
# layer B). Communities in blocks that share their colonisation history down to the separation
# age a_s (simulate_blocks() and the trees of data-raw/s5_spatial.R, gamma = 5, mu = 1, effect
# none or gamma x 1.5 in the focal regime), n = 40 in blocks of 2, 4 or 10 (20, 10 or 4 blocks).
# For the 1-df tests of gamma and of mu among regimes:
# - naive: the LRT against chi-squared(1) with all communities, as in S2;
# - subsample: one community drawn at random from each block, naive LRT on those;
# - Chandler & Bate (2007, Biometrika, doi:10.1093/biomet/asm015): the independence
#   log-likelihood with a sandwich adjustment by blocks. At the estimates of the alternative
#   model, H is minus the Hessian of the summed log-likelihood and J = G / (G - 1) times the sum
#   over the G blocks of the outer product of each block's score (numDeriv). For the contrast c
#   of the log rates of the two regimes, the LRT is divided by
#   lambda = c' H^-1 J H^-1 c / c' H^-1 c and compared with chi-squared(1) ("cb") and with
#   F(1, G - 1) ("cb_f", for few blocks). A fit whose H cannot be inverted (a rate at the
#   boundary) gives NA, counted in `cb_na`.
# REPS replicates per setting (default 200), WORKERS (default 16), JOBS (default
# data-raw/s5c_jobs, ignored by git), OUT (default data-raw/s5_corrections.rds).
src <- readLines("data-raw/s5_spatial.R")
eval(parse(text = src[seq_len(grep("^one <- function", src) - 1)]))
workers <- as.integer(Sys.getenv("WORKERS", "16"))
reps <- as.integer(Sys.getenv("REPS", "200"))
jobs_dir <- Sys.getenv("JOBS", "data-raw/s5c_jobs")
out <- Sys.getenv("OUT", "data-raw/s5_corrections.rds")
dir.create(jobs_dir, showWarnings = FALSE, recursive = TRUE)

# log-likelihood of each community at the working parameters `par` of a shared fit
ll_rows <- function(fit, par) {
  a <- attr(fit, "damocles")
  post <- damocles_postorder(a$tree)
  r <- damocles_par_rates(par, a$design)
  n_edge <- nrow(a$tree$edge)
  g <- rep_len(r$gamma, n_edge)
  m <- rep_len(r$mu, n_edge)
  damocles_prune(a$presence, post$tree, g[post$index], m[post$index], a$root, r$sympatry,
                 damocles_root_rates(a$tree, g, m))
}
cb_lambda <- function(fit, block, contrast) {
  par <- attr(fit, "damocles")$par[[1]]
  total <- function(p) sum(ll_rows(fit, p))
  H <- -numDeriv::hessian(total, par)
  S <- numDeriv::jacobian(function(p) as.numeric(rowsum(ll_rows(fit, p), block)), par)
  G <- nrow(S)
  J <- G / (G - 1) * crossprod(S)
  Hi <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(Hi)) return(NA_real_)
  v_naive <- drop(t(contrast) %*% Hi %*% contrast)
  v_robust <- drop(t(contrast) %*% Hi %*% J %*% Hi %*% contrast)
  if (!is.finite(v_naive) || v_naive <= 0 || !is.finite(v_robust) || v_robust <= 0) return(NA_real_)
  v_robust / v_naive
}
tests <- function(comm, tree, regimes, block = NULL) {
  null <- quiet(damocles_fit(comm, tree, shared = TRUE))
  start <- c(null$gamma[1], null$mu[1])
  fit <- function(vary) quiet(damocles_fit(comm, tree, shared = TRUE, regimes = regimes,
                                           vary = vary, start = start))
  g <- fit("gamma")
  m <- fit("mu")
  ll <- function(x) attr(x, "loglik_total")
  lrt <- c(gamma = max(2 * (ll(g) - ll(null)), 0), mu = max(2 * (ll(m) - ll(null)), 0))
  out <- list(lrt = lrt, p = stats::pchisq(lrt, 1, lower.tail = FALSE))
  if (!is.null(block)) {
    lambda <- c(gamma = cb_lambda(g, block, c(-1, 1, 0)), mu = cb_lambda(m, block, c(0, -1, 1)))
    G <- length(unique(block))
    out$lambda <- lambda
    out$p_cb <- stats::pchisq(lrt / lambda, 1, lower.tail = FALSE)
    out$p_cb_f <- stats::pf(lrt / lambda, 1, G - 1, lower.tail = FALSE)
  }
  out
}
one <- function(s) {
  sim <- focal_tree(112, 0.25)
  focal <- sim$regimes != "0"
  gamma <- rep(5, length(focal))
  mu <- rep(1, length(focal))
  if (s$effect == "gamma") gamma[focal] <- 5 * 1.5
  block <- rep(seq_len(s$n_comm / s$block_size), each = s$block_size)
  comm <- simulate_blocks(sim$tree, gamma, mu, block, s$a_s)
  full <- tests(comm, sim$tree, sim$regimes, block)
  pick <- vapply(split(seq_along(block), block), function(i) i[sample.int(length(i), 1)], 1L)
  sub <- tests(comm[pick, , drop = FALSE], sim$tree, sim$regimes)
  data.frame(naive_gamma = full$p[["gamma"]], naive_mu = full$p[["mu"]],
             sub_gamma = sub$p[["gamma"]], sub_mu = sub$p[["mu"]],
             cb_gamma = full$p_cb[["gamma"]], cb_mu = full$p_cb[["mu"]],
             cbf_gamma = full$p_cb_f[["gamma"]], cbf_mu = full$p_cb_f[["mu"]],
             lambda_gamma = full$lambda[["gamma"]], lambda_mu = full$lambda[["mu"]])
}

settings <- expand.grid(n_comm = 40, block_size = c(2, 4, 10), a_s = c(Inf, 0.1, 0.05, 0.02),
                        effect = c("none", "gamma"), stringsAsFactors = FALSE)
only <- Sys.getenv("SETTINGS", "")
run <- seq_len(nrow(settings))
if (nzchar(only)) run <- intersect(run, as.integer(strsplit(only, ",")[[1]]))
jobs <- expand.grid(setting = run, rep = seq_len(reps))
jobs$file <- file.path(jobs_dir, sprintf("c%02d_r%04d.rds", jobs$setting, jobs$rep))
todo <- which(!file.exists(jobs$file))
cat("S5 corrections:", length(todo), "of", nrow(jobs), "replicates on", workers, "workers\n")
started <- Sys.time()
invisible(parallel::mclapply(todo, function(j) {
  s <- settings[jobs$setting[j], ]
  set.seed(5050 + 10000 * jobs$setting[j] + jobs$rep[j])
  t0 <- Sys.time()
  r <- tryCatch(one(s), error = function(e) data.frame(error = conditionMessage(e)))
  if (is.null(r$error)) r$error <- NA_character_
  r$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(cbind(s, rep = jobs$rep[j], r, row.names = NULL), jobs$file[j])
  NULL
}, mc.cores = workers, mc.preschedule = FALSE))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
files <- jobs$file[file.exists(jobs$file)]
cols <- c("naive_gamma", "naive_mu", "sub_gamma", "sub_mu", "cb_gamma", "cb_mu", "cbf_gamma",
          "cbf_mu", "lambda_gamma", "lambda_mu", "error", "seconds")
results <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f)
  for (k in setdiff(cols, names(x))) x[[k]] <- NA
  x[, c(names(settings), "rep", cols)]
}))
keys <- c("block_size", "a_s", "effect")
rate <- function(p) mean(p < 0.05, na.rm = TRUE)
summary <- do.call(rbind, lapply(split(results, results[, keys], drop = TRUE), function(d)
  data.frame(d[1, keys], blocks = 40 / d$block_size[1],
             naive_g = rate(d$naive_gamma), sub_g = rate(d$sub_gamma), cb_g = rate(d$cb_gamma),
             cbf_g = rate(d$cbf_gamma), naive_m = rate(d$naive_mu), sub_m = rate(d$sub_mu),
             cb_m = rate(d$cb_mu), cbf_m = rate(d$cbf_mu),
             lambda_g = stats::median(d$lambda_gamma, na.rm = TRUE),
             lambda_m = stats::median(d$lambda_mu, na.rm = TRUE),
             cb_na = mean(is.na(d$lambda_gamma) | is.na(d$lambda_mu)),
             errors = sum(!is.na(d$error)))))
summary <- summary[order(summary$effect, summary$block_size, -summary$a_s), ]
cat("S5 corrections:", nrow(results), "replicates in", round(elapsed, 1), "min\n")
options(width = 200)
print(summary, row.names = FALSE, digits = 3)
saveRDS(list(settings = settings, results = results, summary = summary, minutes = elapsed), out)
