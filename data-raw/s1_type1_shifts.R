# S1 of the DAMOCLES paper (dev/plano-artigo-damocles.md): type I error of the tests of
# colonisation and local extinction among regimes when the tree has shifts in the rate of
# speciation but the communities assemble with constant rates.
# 1. Trees: sim_shift_tree() (internal, R/sim_shift_tree.R; validated in
#    data-raw/validation_shift_tree.R): birth-death with speciation 1 and relative extinction
#    0.5, `n_shifts` heritable shifts that multiply speciation by `ratio` (nested shifts
#    compound), conditioned on `n_tips` living lineages by the general sampling approach, every
#    regime with at least 5% of the tips, height one.
# 2. Communities: simulate_damocles() with gamma = 5 and mu = 1 on every edge (no effect of the
#    regimes; mean occupancy about 0.28, as in the 112 species of Benicio et al. 2021),
#    `n_comm` communities per tree.
# 3. Tests, likelihood-ratio against constant rates shared by the communities (damocles_lrt()):
#    - true regimes, gamma and mu free per regime (vary = "both", df = 2 n_shifts);
#    - true regimes, gamma only (df = n_shifts);
#    - DR per edge (dr_rates()) as covariate of gamma (df = 1);
#    the alternatives start from the estimates of the constant model (one start instead of
#    four; same optimum in a check on a tree of 250 tips, half the time);
#    and, as a separate set of jobs (MODE=partition), damocles_partition_test() with the true
#    regimes (vary = "both", 49 random partitions) in the first PARTITION replicates of each
#    setting, on the same tree and communities as the likelihood-ratio tests of that replicate
#    (same seed); random clades within a factor 2 of the size of the true ones (tolerance = 2)
#    with one shift and 3 with three, as nested regimes often find no placement within 2.
# Settings: n_tips 112 or 250, n_shifts 1 or 3, ratio 2 or 4, n_comm 14 or 40, plus a control
# with ratio 1 (shifts with no change in rate). REPS replicates per setting (default 500),
# PARTITION replicates with the partition test (default 50), WORKERS processes (default 6).
# MODE: "lrt" (default) runs the likelihood-ratio tests, "partition" the partition tests; a
# partition test costs about as much as 100 likelihood-ratio replicates, so the two run apart,
# and its 49 random partitions are fitted on 24 cores (fixed, as the random streams depend on
# the number of cores), with WORKERS %/% 12 replicates at a time (twice as many processes as
# cores, to fill the cores left idle while a replicate fits its constant and observed models
# alone), so that a replicate ends in minutes and survives spot interruptions.
# Each replicate is saved in JOBS (default data-raw/s1_jobs, ignored by git; files sXX_rYYYY
# for the likelihood-ratio tests, pXX_rYYYY for the partition tests) with its own seed, and
# skipped on a rerun, so a run cut short (a spot instance on AWS; data-raw/aws/) resumes.
# When every replicate of the mode is done, the results go to OUT (default
# data-raw/s1_type1_shifts.rds, or data-raw/s1_partition.rds with MODE=partition); the summary
# gives the rejection rate at 5% with a 95% binomial interval (Clopper-Pearson).
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE))
workers <- as.integer(Sys.getenv("WORKERS", Sys.getenv("CORES", "6")))
reps <- as.integer(Sys.getenv("REPS", "500"))
n_partition <- as.integer(Sys.getenv("PARTITION", "50"))
jobs_dir <- Sys.getenv("JOBS", "data-raw/s1_jobs")
mode <- match.arg(Sys.getenv("MODE", "lrt"), c("lrt", "partition"))
out <- Sys.getenv("OUT", if (mode == "lrt") "data-raw/s1_type1_shifts.rds" else
  "data-raw/s1_partition.rds")
RNGkind("L'Ecuyer-CMRG")
dir.create(jobs_dir, showWarnings = FALSE, recursive = TRUE)

# 2-3. One replicate ---------------------------------------------------------------------------
p_value <- function(null, alternative) damocles_lrt(null, alternative)$p_value[2]
quiet <- function(x) suppressWarnings(x)
simulate <- function(setting) {
  sim <- sim_shift_tree(setting$n_tips, setting$n_shifts, setting$ratio)
  list(sim = sim, comm = simulate_damocles(sim$tree, 5, 1, n = setting$n_comm))
}
one_lrt <- function(setting) {
  d <- simulate(setting)
  sim <- d$sim
  tree <- sim$tree
  comm <- d$comm
  null <- quiet(damocles_fit(comm, tree, shared = TRUE))
  start <- c(null$gamma[1], null$mu[1])
  both <- quiet(damocles_fit(comm, tree, shared = TRUE, regimes = sim$regimes, vary = "both",
                             start = start))
  gamma <- quiet(damocles_fit(comm, tree, shared = TRUE, regimes = sim$regimes, vary = "gamma",
                              start = start))
  dr <- quiet(damocles_fit(comm, tree, shared = TRUE, edge_covariate = dr_rates(tree, "edge"),
                           vary = "gamma", start = start))
  data.frame(
    p_both = p_value(null, both), p_gamma = p_value(null, gamma), p_dr = p_value(null, dr),
    shifted_tips = mean(sim$tip_regime > 0), shifted_edges = mean(sim$regimes != "0"),
    convergence = max(null$convergence, both$convergence, gamma$convergence, dr$convergence))
}
partition_cores <- 24
one_partition <- function(setting) {
  d <- simulate(setting)
  test <- quiet(damocles_partition_test(d$comm, d$sim$tree, d$sim$regimes, n_rand = 49,
                                        vary = "both",
                                        tolerance = if (setting$n_shifts > 1) 3 else 2,
                                        cores = partition_cores))
  data.frame(p_partition = test$p_random, p_partition_chisq = test$p_chisq,
             shifted_tips = mean(d$sim$tip_regime > 0))
}
failed <- function(e, columns) {
  r <- as.data.frame(as.list(stats::setNames(rep(NA_real_, length(columns)), columns)))
  r$error <- conditionMessage(e)
  r
}

settings <- expand.grid(n_tips = c(112, 250), n_shifts = c(1, 3), ratio = c(2, 4),
                        n_comm = c(14, 40))
settings <- rbind(settings, data.frame(n_tips = c(112, 250), n_shifts = 1, ratio = 1, n_comm = 14))
jobs <- expand.grid(setting = seq_len(nrow(settings)),
                    rep = seq_len(if (mode == "lrt") reps else n_partition))
jobs$file <- file.path(jobs_dir, sprintf("%s%02d_r%04d.rds", if (mode == "lrt") "s" else "p",
                                         jobs$setting, jobs$rep))
# The slowest first: large trees, many communities and many shifts.
cost <- settings$n_tips[jobs$setting] / 112 * settings$n_comm[jobs$setting] / 14 *
  (1 + settings$n_shifts[jobs$setting])
todo <- which(!file.exists(jobs$file))
todo <- todo[order(-cost[todo])]
started <- Sys.time()
cat("S1 (", mode, "): ", length(todo), " of ", nrow(jobs), " replicates to run on ", workers,
    " workers", if (mode == "partition") paste0(" (", partition_cores, " per replicate)"), "\n",
    sep = "")
columns <- if (mode == "lrt") {
  c("p_both", "p_gamma", "p_dr", "shifted_tips", "shifted_edges", "convergence")
} else {
  c("p_partition", "p_partition_chisq", "shifted_tips")
}
invisible(parallel::mclapply(todo, function(j) {
  s <- settings[jobs$setting[j], ]
  # The same seed in both modes: the partition test sees the trees and communities of the
  # likelihood-ratio tests of the same replicate.
  set.seed(101 + 10000 * jobs$setting[j] + jobs$rep[j])
  r <- tryCatch(if (mode == "lrt") one_lrt(s) else one_partition(s),
                error = function(e) failed(e, columns))
  if (is.null(r$error)) r$error <- NA_character_
  saveRDS(cbind(s, rep = jobs$rep[j], r), jobs$file[j])
  NULL
}, mc.cores = if (mode == "lrt") workers else max(1, workers %/% (partition_cores / 2)),
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
tests <- if (mode == "lrt") c("p_both", "p_gamma", "p_dr") else
  c("p_partition", "p_partition_chisq")
summary <- do.call(rbind, lapply(split(results, results[, names(settings)], drop = TRUE),
  function(d) {
    rates <- as.data.frame(lapply(d[tests], rate))
    names(rates) <- sub("^p_", "", tests)
    extra <- data.frame(shifted_tips = round(mean(d$shifted_tips, na.rm = TRUE), 2),
                        errors = sum(!is.na(d$error)))
    if (mode == "lrt") extra$not_converged <- sum(d$convergence != 0, na.rm = TRUE)
    cbind(d[1, names(settings)], rates, extra)
  }))
summary <- summary[order(summary$ratio, summary$n_tips, summary$n_shifts, summary$n_comm), ]
cat("S1 (", mode, "): ", nrow(results), " replicates in ", round(elapsed, 1), " min\n", sep = "")
print(summary, row.names = FALSE)
saveRDS(list(settings = settings, results = results, summary = summary, minutes = elapsed), out)
