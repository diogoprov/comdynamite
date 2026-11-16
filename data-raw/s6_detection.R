# Imperfect detection in the DAMOCLES paper (dev/plano-artigo-damocles.md, variant of S6): how
# false absences in the community data bias the rates and change the size, power and attribution
# of the tests of colonisation (gamma) and local extinction (mu) among regimes.
# 1. Trees and communities as in S2 (data-raw/s2_power.R): sim_shift_tree(n_tips, 1, 2) redrawn
#    until the focal regime holds 0.25 of the tips (factor 0.8-1.2), gamma = 5 and mu = 1 in the
#    background and, in the focal regime, gamma or mu multiplied by `ratio` (ratio 1: no effect,
#    for the size of the tests).
# 2. False absences: each presence becomes an absence with probability q, independently, in
#    0, 0.1, 0.25 and 0.4 (q = 0.2-0.3 is what the comparison of Benicio et al. 2021 with
#    ATLANTIC AMPHIBIANS suggests; data-raw/benicio2021_omission.R). The four levels are applied
#    to the same tree and communities (paired), each with its own mask.
# 3. On each data set the four models of S2 layer C: constant rates, gamma free per regime, mu
#    free, both free (damocles_fit(shared = TRUE), alternatives started from the constant
#    estimates). Kept: the constant estimates of gamma and mu (bias), the per-regime estimates of
#    the model that frees both, the log-likelihoods, the LRT p-values against the constant model
#    and the choice between the gamma and mu models (higher likelihood).
# Settings: n_comm 14 or 40, n_tips 112, ratio 1 (effect none), 1.5 or 3 with effect on gamma or
# mu: 10 settings x 4 levels of q. REPS replicates per setting (default 200), WORKERS processes
# (default 16). Each replicate (all four q) is saved in JOBS (default data-raw/s6_jobs, ignored by
# git) and skipped on a rerun. Results to OUT (default data-raw/s6_detection.rds).
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE))
workers <- as.integer(Sys.getenv("WORKERS", "16"))
reps <- as.integer(Sys.getenv("REPS", "200"))
jobs_dir <- Sys.getenv("JOBS", "data-raw/s6_jobs")
out <- Sys.getenv("OUT", "data-raw/s6_detection.rds")
q_levels <- c(0, 0.1, 0.25, 0.4)
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
simulate <- function(setting) {
  sim <- focal_tree(setting$n_tips, 0.25)
  focal <- sim$regimes != "0"
  gamma <- rep(5, length(focal))
  mu <- rep(1, length(focal))
  if (setting$effect == "gamma") gamma[focal] <- 5 * setting$ratio
  if (setting$effect == "mu") mu[focal] <- setting$ratio
  list(sim = sim, comm = simulate_damocles(sim$tree, gamma, mu, n = setting$n_comm))
}
fit_all <- function(comm, sim) {
  tree <- sim$tree
  regimes <- sim$regimes
  null <- quiet(damocles_fit(comm, tree, shared = TRUE))
  start <- c(null$gamma[1], null$mu[1])
  fit <- function(vary) quiet(damocles_fit(comm, tree, shared = TRUE, regimes = regimes,
                                           vary = vary, start = start))
  both <- fit("both")
  g <- fit("gamma")
  m <- fit("mu")
  ll <- function(x) attr(x, "loglik_total")
  p <- function(alt, df) stats::pchisq(pmax(2 * (ll(alt) - ll(null)), 0), df, lower.tail = FALSE)
  data.frame(gamma_null = null$gamma[1], mu_null = null$mu[1],
             gamma_0 = both$gamma_0[1], gamma_1 = both$gamma_1[1],
             mu_0 = both$mu_0[1], mu_1 = both$mu_1[1],
             ll_null = ll(null), ll_gamma = ll(g), ll_mu = ll(m), ll_both = ll(both),
             p_both = p(both, 2), p_gamma = p(g, 1), p_mu = p(m, 1),
             pick = if (ll(g) >= ll(m)) "gamma" else "mu",
             richness = mean(rowSums(comm)),
             convergence = max(null$convergence, both$convergence, g$convergence, m$convergence))
}
one <- function(setting) {
  d <- simulate(setting)
  do.call(rbind, lapply(q_levels, function(q) {
    comm <- d$comm
    comm[comm > 0 & stats::runif(length(comm)) < q] <- 0L
    cbind(q = q, fit_all(comm, d$sim))
  }))
}

settings <- rbind(
  expand.grid(n_comm = c(14, 40), n_tips = 112, ratio = 1, effect = "none",
              stringsAsFactors = FALSE),
  expand.grid(n_comm = c(14, 40), n_tips = 112, ratio = c(1.5, 3), effect = c("gamma", "mu"),
              stringsAsFactors = FALSE))
only <- Sys.getenv("SETTINGS", "")
run <- seq_len(nrow(settings))
if (nzchar(only)) run <- intersect(run, as.integer(strsplit(only, ",")[[1]]))
jobs <- expand.grid(setting = run, rep = seq_len(reps))
jobs$file <- file.path(jobs_dir, sprintf("d%02d_r%04d.rds", jobs$setting, jobs$rep))
todo <- which(!file.exists(jobs$file))
todo <- todo[order(-settings$n_comm[jobs$setting[todo]])]
cat("S6 detection:", length(todo), "of", nrow(jobs), "replicates on", workers, "workers\n")
started <- Sys.time()
invisible(parallel::mclapply(todo, function(j) {
  s <- settings[jobs$setting[j], ]
  set.seed(606 + 10000 * jobs$setting[j] + jobs$rep[j])
  t0 <- Sys.time()
  r <- tryCatch(one(s), error = function(e) data.frame(q = NA, error = conditionMessage(e)))
  if (is.null(r$error)) r$error <- NA_character_
  r$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(cbind(s, rep = jobs$rep[j], r, row.names = NULL), jobs$file[j])
  NULL
}, mc.cores = workers, mc.preschedule = FALSE))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
done <- file.exists(jobs$file)
if (!all(done)) stop(sum(!done), " replicates missing in ", jobs_dir, call. = FALSE)
results <- do.call(rbind, lapply(jobs$file, function(f) {
  x <- readRDS(f)
  for (k in setdiff(c("gamma_null", "mu_null", "gamma_0", "gamma_1", "mu_0", "mu_1", "ll_null",
                      "ll_gamma", "ll_mu", "ll_both", "p_both", "p_gamma", "p_mu", "pick",
                      "richness", "convergence"), names(x))) x[[k]] <- NA
  x[, c(names(settings), "rep", "q", "gamma_null", "mu_null", "gamma_0", "gamma_1", "mu_0",
        "mu_1", "ll_null", "ll_gamma", "ll_mu", "ll_both", "p_both", "p_gamma", "p_mu", "pick",
        "richness", "convergence", "error", "seconds")]
}))

# Summary: median constant estimates (true background 5 and 1), the per-regime ratio of the
# simulated rate, rejection at 5% of the matched test (gamma -> p_gamma, mu -> p_mu, none -> all
# three) and the share of detected effects (best single-rate model significant) attributed to the
# right rate.
results$p_match <- with(results, ifelse(effect == "gamma", p_gamma,
                                        ifelse(effect == "mu", p_mu, p_both)))
results$best_p <- with(results, stats::pchisq(pmax(2 * (pmax(ll_gamma, ll_mu) - ll_null), 0), 1,
                                              lower.tail = FALSE))
keys <- c("n_comm", "ratio", "effect", "q")
summary <- do.call(rbind, lapply(split(results, results[, keys], drop = TRUE), function(d) {
  det <- d$best_p < 0.05 & !is.na(d$best_p)
  data.frame(d[1, keys],
             gamma_null = stats::median(d$gamma_null, na.rm = TRUE),
             mu_null = stats::median(d$mu_null, na.rm = TRUE),
             gamma_ratio = stats::median(d$gamma_1 / d$gamma_0, na.rm = TRUE),
             mu_ratio = stats::median(d$mu_1 / d$mu_0, na.rm = TRUE),
             p_both = mean(d$p_both < 0.05, na.rm = TRUE),
             p_gamma = mean(d$p_gamma < 0.05, na.rm = TRUE),
             p_mu = mean(d$p_mu < 0.05, na.rm = TRUE),
             detected = mean(det),
             right_rate = if (d$effect[1] == "none" || !any(det)) NA else
               mean(d$pick[det] == d$effect[1]),
             richness = round(mean(d$richness, na.rm = TRUE), 1),
             errors = sum(!is.na(d$error)), not_converged = sum(d$convergence != 0, na.rm = TRUE))
}))
summary <- summary[order(summary$effect, summary$ratio, summary$n_comm, summary$q), ]
cat("S6 detection:", nrow(results), "fits in", round(elapsed, 1), "min\n")
options(width = 200)
print(summary, row.names = FALSE, digits = 3)
saveRDS(list(settings = settings, results = results, summary = summary, minutes = elapsed), out)
