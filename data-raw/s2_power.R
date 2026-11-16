# S2 of the DAMOCLES paper (dev/plano-artigo-damocles.md): power of the tests of colonisation
# and local extinction among regimes, with the regimes known (the true ones).
# 1. Trees: sim_shift_tree() (internal, R/sim_shift_tree.R) as in S1: birth-death with
#    speciation 1 and relative extinction 0.5, one heritable shift that doubles speciation,
#    conditioned on `n_tips` living lineages, height one. Trees are drawn again until the
#    shifted (focal) regime holds `fraction` of the tips within a factor 0.8-1.2.
# 2. Communities: simulate_damocles() with the rates of S1 (gamma = 5, mu = 1) in the
#    background and, in the focal regime, gamma, mu or both multiplied by `ratio` (both: the
#    same ratio, so occupancy stays about the same and only turnover changes), `n_comm`
#    communities per tree.
# 3. Tests, likelihood-ratio against constant rates shared by the communities (damocles_lrt()),
#    with the true regimes: gamma and mu free per regime (p_both, df 2), gamma only (p_gamma,
#    df 1) and mu only (p_mu, df 1); the alternatives start from the constant estimates. The
#    estimated ratio of the simulated rate (ratio_hat, from the model that frees that rate, or
#    gamma for "both") is kept for bias.
#    MODE=boot runs instead the parametric bootstrap of damocles_partition_test() (n_rand = 0,
#    N_BOOT simulations under the constant rates, default 49, vary = "both", 24 cores per
#    replicate) on the same trees and communities (same seed); it reports p_boot and the
#    chi-square p of the same statistic (p_boot_chisq).
# Settings (MODE=chisq): grid A, n_comm 5, 10, 14, 20 or 40, n_tips 112 or 250, ratio 1.1, 1.25,
# 1.5, 2 or 3, effect on gamma, mu or both, fraction 0.25; plus a sweep of fraction 0.1 and 0.5 with
# n_tips 112, n_comm 14 or 40, effect on gamma. Type I error (ratio 1) comes from S1.
# MODE=boot: grid B, n_tips 112, n_comm 14 or 40, every ratio, effect on gamma,
# fraction 0.25 (the same settings and seeds as in grid A).
# REPS replicates per setting (default 200 with chisq, 100 with boot), WORKERS processes
# (default 6). SETTINGS (optional, comma-separated setting numbers) restricts the run, for
# pilots; numbers, seeds and file names stay those of the full grid.
# Each replicate is saved in JOBS (default data-raw/s2_jobs, ignored by git; files cXX_rYYYY
# or bXX_rYYYY) with its own seed and skipped on a rerun, so a run cut short resumes (spot
# instances on AWS; data-raw/aws/). When every replicate is done the results go to OUT
# (default data-raw/s2_power.rds or data-raw/s2_boot.rds); the summary gives the rejection
# rate at 5% with a 95% binomial interval (Clopper-Pearson).
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE))
workers <- as.integer(Sys.getenv("WORKERS", "6"))
mode <- match.arg(Sys.getenv("MODE", "chisq"), c("chisq", "boot"))
reps <- as.integer(Sys.getenv("REPS", if (mode == "chisq") "200" else "100"))
n_boot <- as.integer(Sys.getenv("N_BOOT", "49"))
jobs_dir <- Sys.getenv("JOBS", "data-raw/s2_jobs")
out <- Sys.getenv("OUT", c(chisq = "data-raw/s2_power.rds", boot = "data-raw/s2_boot.rds")[[mode]])
RNGkind("L'Ecuyer-CMRG")
dir.create(jobs_dir, showWarnings = FALSE, recursive = TRUE)

# 1-2. Trees and communities -----------------------------------------------------------------
quiet <- function(x) suppressWarnings(x)
focal_tree <- function(n_tips, fraction, max_tries = 1000) {
  for (i in seq_len(max_tries)) {
    sim <- sim_shift_tree(n_tips, 1, 2)
    f <- mean(sim$tip_regime > 0)
    if (f >= 0.8 * fraction && f <= 1.2 * fraction) return(sim)
  }
  stop("no tree with the focal fraction in ", max_tries, " tries", call. = FALSE)
}
simulate <- function(setting) {
  sim <- focal_tree(setting$n_tips, setting$fraction)
  focal <- sim$regimes != "0"
  gamma <- rep(5, length(focal))
  mu <- rep(1, length(focal))
  if (setting$effect %in% c("gamma", "both")) gamma[focal] <- 5 * setting$ratio
  if (setting$effect %in% c("mu", "both")) mu[focal] <- setting$ratio
  # simulate_damocles() requires equal rates on the two edges leaving the root; the focal clade
  # never includes the root, as the shift happens after the first speciation.
  list(sim = sim, comm = simulate_damocles(sim$tree, gamma, mu, n = setting$n_comm))
}

# 3. One replicate ---------------------------------------------------------------------------
p_value <- function(null, alternative) damocles_lrt(null, alternative)$p_value[2]
one_chisq <- function(setting) {
  d <- simulate(setting)
  tree <- d$sim$tree
  regimes <- d$sim$regimes
  null <- quiet(damocles_fit(d$comm, tree, shared = TRUE))
  start <- c(null$gamma[1], null$mu[1])
  fit <- function(vary) quiet(damocles_fit(d$comm, tree, shared = TRUE, regimes = regimes,
                                           vary = vary, start = start))
  both <- fit("both")
  gamma <- fit("gamma")
  mu <- fit("mu")
  ratio_hat <- switch(setting$effect,
                      gamma = gamma$gamma_1[1] / gamma$gamma_0[1],
                      mu = mu$mu_1[1] / mu$mu_0[1],
                      both = both$gamma_1[1] / both$gamma_0[1])
  data.frame(p_both = p_value(null, both), p_gamma = p_value(null, gamma),
             p_mu = p_value(null, mu), ratio_hat = ratio_hat,
             focal_tips = mean(d$sim$tip_regime > 0),
             convergence = max(null$convergence, both$convergence, gamma$convergence,
                               mu$convergence))
}
boot_cores <- 24
one_boot <- function(setting) {
  d <- simulate(setting)
  test <- quiet(damocles_partition_test(d$comm, d$sim$tree, d$sim$regimes, n_rand = 0,
                                        n_boot = n_boot, vary = "both", cores = boot_cores))
  data.frame(p_boot = test$p_boot, p_boot_chisq = test$p_chisq,
             n_boot_ok = length(test$null_boot), focal_tips = mean(d$sim$tip_regime > 0))
}
failed <- function(e, columns) {
  r <- as.data.frame(as.list(stats::setNames(rep(NA_real_, length(columns)), columns)))
  r$error <- conditionMessage(e)
  r
}

# Settings and jobs ----------------------------------------------------------------------------
grid_a <- expand.grid(n_comm = c(5, 10, 14, 20, 40), n_tips = c(112, 250),
                      ratio = c(1.1, 1.25, 1.5, 2, 3), effect = c("gamma", "mu", "both"),
                      fraction = 0.25, stringsAsFactors = FALSE)
sweep <- expand.grid(n_comm = c(14, 40), n_tips = 112, ratio = c(1.1, 1.25, 1.5, 2, 3),
                     effect = "gamma",
                     fraction = c(0.1, 0.5), stringsAsFactors = FALSE)
settings <- rbind(grid_a, sweep)
in_b <- settings$n_tips == 112 & settings$n_comm %in% c(14, 40) & settings$effect == "gamma" &
  settings$fraction == 0.25
run <- if (mode == "chisq") seq_len(nrow(settings)) else which(in_b)
only <- Sys.getenv("SETTINGS", "")
if (nzchar(only)) run <- intersect(run, as.integer(strsplit(only, ",")[[1]]))
jobs <- expand.grid(setting = run, rep = seq_len(reps))
jobs$file <- file.path(jobs_dir, sprintf("%s%03d_r%04d.rds", c(chisq = "c", boot = "b")[[mode]],
                                         jobs$setting, jobs$rep))
# The slowest first: large trees and many communities.
cost <- settings$n_tips[jobs$setting] * settings$n_comm[jobs$setting]
todo <- which(!file.exists(jobs$file))
todo <- todo[order(-cost[todo])]
started <- Sys.time()
cat("S2 (", mode, "): ", length(todo), " of ", nrow(jobs), " replicates to run on ", workers,
    " workers", if (mode == "boot") paste0(" (", boot_cores, " cores per replicate)"), "\n",
    sep = "")
columns <- switch(mode,
  chisq = c("p_both", "p_gamma", "p_mu", "ratio_hat", "focal_tips", "convergence"),
  boot = c("p_boot", "p_boot_chisq", "n_boot_ok", "focal_tips"))
invisible(parallel::mclapply(todo, function(j) {
  s <- settings[jobs$setting[j], ]
  # The same seed in both modes: the bootstrap sees the trees and communities of grid A.
  set.seed(202 + 10000 * jobs$setting[j] + jobs$rep[j])
  t0 <- Sys.time()
  r <- tryCatch(switch(mode, chisq = one_chisq(s), boot = one_boot(s)),
                error = function(e) failed(e, columns))
  if (is.null(r$error)) r$error <- NA_character_
  r$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(cbind(s, rep = jobs$rep[j], r), jobs$file[j])
  NULL
}, mc.cores = if (mode == "chisq") workers else max(1, workers %/% (boot_cores / 2)),
   mc.preschedule = FALSE))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
done <- file.exists(jobs$file)
if (!all(done)) stop(sum(!done), " replicates missing in ", jobs_dir, call. = FALSE)
results <- do.call(rbind, lapply(jobs$file, readRDS))

rate <- function(p) {
  p <- p[!is.na(p)]
  k <- sum(p < 0.05)
  ci <- if (length(p)) stats::binom.test(k, length(p))$conf.int else c(NA, NA)
  sprintf("%.3f [%.3f, %.3f] (n=%d)", k / max(length(p), 1), ci[1], ci[2], length(p))
}
tests <- switch(mode, chisq = c("p_both", "p_gamma", "p_mu"), boot = c("p_boot", "p_boot_chisq"))
keys <- c("n_comm", "n_tips", "ratio", "effect", "fraction")
summary <- do.call(rbind, lapply(split(results, results[, keys], drop = TRUE), function(d) {
  rates <- as.data.frame(lapply(d[tests], rate))
  names(rates) <- sub("^p_", "", tests)
  extra <- data.frame(focal_tips = round(mean(d$focal_tips, na.rm = TRUE), 2),
                      errors = sum(!is.na(d$error)),
                      seconds = round(stats::median(d$seconds), 1))
  if (mode == "chisq") {
    extra$ratio_hat <- round(stats::median(d$ratio_hat, na.rm = TRUE), 2)
    extra$not_converged <- sum(d$convergence != 0, na.rm = TRUE)
  }
  if (mode == "boot") extra$boot_ok <- stats::median(d$n_boot_ok, na.rm = TRUE)
  cbind(d[1, keys], rates, extra)
}))
summary <- summary[order(summary$effect, summary$fraction, summary$n_tips, summary$ratio,
                         summary$n_comm), ]
cat("S2 (", mode, "): ", nrow(results), " replicates in ", round(elapsed, 1), " min\n", sep = "")
print(summary, row.names = FALSE)
saveRDS(list(settings = settings, results = results, summary = summary, minutes = elapsed), out)
