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
#    and, in the first `n_partition` replicates of each setting, damocles_partition_test() with
#    the true regimes (vary = "both", 49 random partitions).
# Settings: n_tips 112 or 250, n_shifts 1 or 3, ratio 2 or 4, n_comm 14 or 40, plus a control
# with ratio 1 (shifts with no change in rate). REPS replicates per setting (default 500),
# PARTITION replicates with the partition test (default 50), CORES workers (default 6).
# Results in data-raw/s1_type1_shifts.rds; the summary gives the rejection rate at 5% with a
# 95% binomial interval (Clopper-Pearson).
suppressMessages(devtools::load_all(".", quiet = TRUE))
cores <- as.integer(Sys.getenv("CORES", "6"))
reps <- as.integer(Sys.getenv("REPS", "500"))
n_partition <- as.integer(Sys.getenv("PARTITION", "50"))
out <- Sys.getenv("OUT", "data-raw/s1_type1_shifts.rds")
RNGkind("L'Ecuyer-CMRG")
set.seed(101)

# 2-3. One replicate ---------------------------------------------------------------------------
p_value <- function(null, alternative) damocles_lrt(null, alternative)$p_value[2]
one <- function(setting, partition) {
  sim <- sim_shift_tree(setting$n_tips, setting$n_shifts, setting$ratio)
  tree <- sim$tree
  comm <- simulate_damocles(tree, 5, 1, n = setting$n_comm)
  quiet <- function(x) suppressWarnings(x)
  null <- quiet(damocles_fit(comm, tree, shared = TRUE))
  start <- c(null$gamma[1], null$mu[1])
  both <- quiet(damocles_fit(comm, tree, shared = TRUE, regimes = sim$regimes, vary = "both",
                             start = start))
  gamma <- quiet(damocles_fit(comm, tree, shared = TRUE, regimes = sim$regimes, vary = "gamma",
                              start = start))
  dr <- quiet(damocles_fit(comm, tree, shared = TRUE, edge_covariate = dr_rates(tree, "edge"),
                           vary = "gamma", start = start))
  result <- data.frame(
    p_both = p_value(null, both), p_gamma = p_value(null, gamma), p_dr = p_value(null, dr),
    p_partition = NA_real_, shifted_tips = mean(sim$tip_regime > 0),
    shifted_edges = mean(sim$regimes != "0"),
    convergence = max(null$convergence, both$convergence, gamma$convergence, dr$convergence))
  if (partition) {
    test <- quiet(damocles_partition_test(comm, tree, sim$regimes, n_rand = 49, vary = "both"))
    result$p_partition <- test$p_random
  }
  result
}

settings <- expand.grid(n_tips = c(112, 250), n_shifts = c(1, 3), ratio = c(2, 4),
                        n_comm = c(14, 40))
settings <- rbind(settings, data.frame(n_tips = c(112, 250), n_shifts = 1, ratio = 1, n_comm = 14))
jobs <- expand.grid(setting = seq_len(nrow(settings)), rep = seq_len(reps))
started <- Sys.time()
results <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  s <- settings[jobs$setting[j], ]
  r <- tryCatch(one(s, jobs$rep[j] <= n_partition), error = function(e) {
    data.frame(p_both = NA, p_gamma = NA, p_dr = NA, p_partition = NA, shifted_tips = NA,
               shifted_edges = NA, convergence = NA, error = conditionMessage(e))
  })
  if (is.null(r$error)) r$error <- NA_character_
  cbind(s, rep = jobs$rep[j], r)
}, mc.cores = cores, mc.preschedule = FALSE)
results <- do.call(rbind, results)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))

rate <- function(p) {
  p <- p[!is.na(p)]
  k <- sum(p < 0.05)
  ci <- if (length(p)) stats::binom.test(k, length(p))$conf.int else c(NA, NA)
  sprintf("%.3f [%.3f, %.3f] (n=%d)", k / max(length(p), 1), ci[1], ci[2], length(p))
}
summary <- do.call(rbind, lapply(split(results, results[, names(settings)], drop = TRUE),
  function(d) {
    cbind(d[1, names(settings)], both = rate(d$p_both), gamma = rate(d$p_gamma),
          dr = rate(d$p_dr), partition = rate(d$p_partition),
          shifted_tips = round(mean(d$shifted_tips, na.rm = TRUE), 2),
          errors = sum(!is.na(d$error)), not_converged = sum(d$convergence != 0, na.rm = TRUE))
  }))
summary <- summary[order(summary$ratio, summary$n_tips, summary$n_shifts, summary$n_comm), ]
cat("S1:", nrow(results), "replicates in", round(elapsed, 1), "min\n")
print(summary, row.names = FALSE)
saveRDS(list(settings = settings, results = results, summary = summary, minutes = elapsed), out)
