# Recovery of the parameters of simulate_logistic() by a phylogenetic GLMM.
#
# simulate_logistic() follows Appendix C of Peres-Neto et al. (2012), after Ives & Helmus
# (2011): the probability that species j occurs at site t is plogis(a[j] + b[j] * x[t]),
# with a ~ N(intercept_mean, intercept_sd^2) and the slopes b evolving on the tree under
# Brownian motion with Pagel's lambda. Ives & Helmus (2011) analyse such data with a
# binomial phylogenetic GLMM whose slopes along the gradient vary among species with and
# without phylogenetic covariance. This script simulates data with known parameters, fits
# that model with phyr::pglmm() and asks whether the variance components it estimates
# agree with the ones simulated.
#
# Design. A fully balanced tree with 32 species and Grafen branch lengths, as in Appendix
# C, and 40 sites on one gradient. The gradient is evenly spaced and standardised to mean 0
# and variance 1, which pglmm() recommends; the [-15, 15] range of Appendix C with slopes of
# standard deviation 1 gives logits of up to 15 and nearly deterministic presences, and a
# small side check at the end shows what happens then. Scenarios are the three values of
# lambda (0, 0.5, 1) crossed with two site reductions: none (reduction = c(1, 1)) and the
# default U(0.5, 1), which the GLMM does not model. Other parameters are the defaults:
# intercept_mean = 0, intercept_sd = 0.5, slope_sd = 1.
#
# Model. y ~ x + (1 | sp) + (x | sp__), family = "binomial", cov_ranef = list(sp = tree),
# with REML and no observation-level random term. The term (x | sp__) creates two random
# slopes, one with the phylogenetic covariance matrix and one with the identity matrix, so
# the random effects are, in this order, 1|sp, x|sp and x|sp__. The phylogenetic slope is
# tested with pglmm_profile_LRT(), a profile likelihood ratio test conditional on the other
# parameters. In the scenarios with site reduction the model is also fitted with a site
# random intercept, + (1 | site), to see whether that term absorbs the reduction.
#
# Targets. pglmm() standardises the phylogenetic covariance matrix to determinant 1, so its
# phylogenetic variance is multiplied by the diagonal of that matrix (1 / det(C)^(1/n),
# about 7.5 for this tree) to return to the unit-diagonal scale of the simulation. The
# simulator standardises the slopes of each replicate to a sample standard deviation of
# exactly slope_sd, so slope_sd^2 is the realised sample variance, not the variance of the
# process. For a process with covariance s2 * (lambda * C + (1 - lambda) * I) the expected
# sample variance is s2 * n / (n - 1) * (1 - mean of that matrix); the targets for the
# process variances are therefore s2 = slope_sd^2 divided by that factor, split as lambda
# and 1 - lambda. For reference, the same decomposition is also estimated from the true
# slopes with a Gaussian pglmm(b ~ 1 + (1 | sp__)), which shows what 32 species allow when
# the slopes are observed without the binary sampling.
#
# Running. From the package root, with phyr (CRAN, GPL-3; not a dependency of comdynamite)
# installed: Rscript data-raw/validation_pglmm.R. With 100 replicates per scenario and 16
# cores the script takes about 25 minutes (900 binomial fits of 12-16 s each, plus 600
# Gaussian fits of the true slopes); see `timing` in the saved results. Results are
# saved in data-raw/validation_pglmm.rds as a list with the replicate-level data frame, the
# summary tables, the side check and session information.
#
# References (verified against the package documentation and the DOIs):
# Ives, A. R., & Helmus, M. R. (2011). Generalized linear mixed models for phylogenetic
#   analyses of community structure. Ecological Monographs, 81, 511-525.
#   doi:10.1890/10-1264.1
# Li, D., Dinnage, R., Nell, L. A., Helmus, M. R., & Ives, A. R. (2020). phyr: An R package
#   for phylogenetic species-distribution modelling in ecological communities. Methods in
#   Ecology and Evolution, 11, 1455-1463. doi:10.1111/2041-210X.13471
# Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects of spatial
#   contingency and environmental filtering on metacommunity phylogenetics. Ecology, 93,
#   S14-S30. doi:10.1890/11-0494.1

suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(phyr))
library(parallel)

n_rep <- 100
n_sites <- 40
n_sp <- 32
cores <- max(1, min(16, detectCores() - 2))
alpha <- 0.05

tree <- ape::compute.brlen(ape::stree(n_sp, "balanced"), method = "Grafen")
C <- ape::vcv(tree)
phylo_scale <- 1 / det(C)^(1 / n_sp)     # diagonal of the matrix standardised by pglmm()
x <- drop(scale(seq(-1, 1, length.out = n_sites)))

# expected sample variance of a process with covariance lambda * C + (1 - lambda) * I and
# unit process variance, over the n tips
sample_factor <- function(lambda) {
  V <- lambda * C + (1 - lambda) * diag(n_sp)
  n_sp / (n_sp - 1) * (1 - mean(V))
}

# long data frame for pglmm() from a simulation
long_data <- function(sim, x) {
  data.frame(y = as.vector(sim$abundance),
             sp = rep(colnames(sim$abundance), each = length(x)),
             site = rep(rownames(sim$abundance), times = ncol(sim$abundance)),
             x = rep(x, times = ncol(sim$abundance)))
}

# fit the binomial PGLMM; the phylogenetic slope is the third random term (fourth with a
# site intercept, which comes second). Returns NA for a fit that errors. nloptr codes 1 to
# 4 are convergence, 5 is the iteration limit, negative codes are failures.
fit_pglmm <- function(d, site = FALSE) {
  form <- if (site) y ~ x + (1 | sp) + (1 | site) + (x | sp__) else y ~ x + (1 | sp) + (x | sp__)
  out <- tryCatch({
    fit <- pglmm(form, data = d, family = "binomial", cov_ranef = list(sp = tree),
                 add.obs.re = FALSE, REML = TRUE)
    re <- names(fit$random.effects)
    lrt <- tryCatch(pglmm_profile_LRT(fit, re.number = which(re == "x|sp__"))$Pr,
                    error = function(e) NA_real_)
    c(s2_int = unname(fit$s2r[re == "1|sp"]),
      s2_site = if (site) unname(fit$s2r[re == "1|site"]) else NA_real_,
      s2_slope_ind = unname(fit$s2r[re == "x|sp"]),
      s2_slope_phy = unname(fit$s2r[re == "x|sp__"]) * phylo_scale,
      b0 = fit$B[1], b1 = fit$B[2], convcode = fit$convcode, niter = fit$niter,
      p_phylo = lrt, error = 0)
  }, error = function(e) {
    c(s2_int = NA, s2_site = NA, s2_slope_ind = NA, s2_slope_phy = NA, b0 = NA, b1 = NA,
      convcode = NA, niter = NA, p_phylo = NA, error = 1)
  })
  out
}

# decomposition of the true slopes into phylogenetic and independent variance, Gaussian
oracle <- function(b) {
  d <- data.frame(b = b, sp = names(b))
  tryCatch({
    fit <- pglmm(b ~ 1 + (1 | sp__), data = d, family = "gaussian",
                 cov_ranef = list(sp = tree), REML = TRUE)
    # the Gaussian model also has a residual variance, which is independent among species
    # as well, so it is added to the independent component
    re <- names(fit$random.effects)
    c(oracle_ind = unname(fit$s2r[re == "1|sp"]) + fit$s2resid,
      oracle_phy = unname(fit$s2r[re == "1|sp__"]) * phylo_scale)
  }, error = function(e) c(oracle_ind = NA_real_, oracle_phy = NA_real_))
}

one_replicate <- function(lambda, reduction, rep) {
  sim <- simulate_logistic(tree, x, lambda = lambda, reduction = reduction)
  b <- sim$traits$slope[, 1]
  a <- sim$traits$intercept[, 1]
  d <- long_data(sim, x)
  prevalence <- colMeans(sim$abundance)
  t0 <- Sys.time()
  base <- fit_pglmm(d)
  seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  with_site <- if (reduction[1] < 1) fit_pglmm(d, site = TRUE) else NULL
  out <- data.frame(lambda = lambda, reduction = reduction[1], rep = rep,
                    var_a = var(a), var_b = var(b), mean_prevalence = mean(prevalence),
                    n_fixed_species = sum(prevalence %in% c(0, 1)), seconds = seconds,
                    t(base), t(oracle(b)))
  if (!is.null(with_site)) {
    names(with_site) <- paste0("site_", names(with_site))
    out <- cbind(out, t(with_site))
  }
  out
}

scenarios <- expand.grid(lambda = c(0, 0.5, 1), reduction = c(1, 0.5))
set.seed(2011)
RNGkind("L'Ecuyer-CMRG")
t_start <- Sys.time()
replicates <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  lambda <- scenarios$lambda[i]
  reduction <- c(scenarios$reduction[i], 1)
  runs <- mclapply(seq_len(n_rep), function(r) one_replicate(lambda, reduction, r),
                   mc.cores = cores, mc.set.seed = TRUE)
  failed <- vapply(runs, inherits, logical(1), "try-error")
  if (any(failed)) runs <- runs[!failed]
  out <- do.call(rbind, lapply(runs, function(r) {
    # align columns across scenarios with and without the site model
    for (col in c("site_s2_int", "site_s2_site", "site_s2_slope_ind", "site_s2_slope_phy",
                  "site_b0", "site_b1", "site_convcode", "site_niter", "site_p_phylo",
                  "site_error")) if (is.null(r[[col]])) r[[col]] <- NA_real_
    r
  }))
  message(sprintf("lambda = %.1f, reduction = %.1f: %d replicates, %.0f s",
                  lambda, reduction[1], nrow(out),
                  as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  out
}))
timing <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
message(sprintf("main runs: %.1f min on %d cores", timing, cores))

# summaries: targets and estimates per scenario
converged <- function(code) !is.na(code) & code >= 1 & code <= 4
summarise_scenario <- function(d) {
  lambda <- d$lambda[1]
  s2 <- 1 / sample_factor(lambda)
  ok <- converged(d$convcode)
  q <- function(v) sprintf("%.2f (%.2f-%.2f)", median(v, na.rm = TRUE),
                           quantile(v, 0.1, na.rm = TRUE), quantile(v, 0.9, na.rm = TRUE))
  phy_frac <- d$s2_slope_phy / (d$s2_slope_phy + d$s2_slope_ind)
  data.frame(lambda = lambda, reduction = d$reduction[1], n = nrow(d),
             n_error = sum(d$error == 1, na.rm = TRUE), n_not_converged = sum(!ok & d$error == 0),
             target_s2_int = 0.25, realised_s2_int = q(d$var_a), est_s2_int = q(d$s2_int),
             target_s2_slope_phy = round(lambda * s2, 2), est_s2_slope_phy = q(d$s2_slope_phy),
             oracle_s2_slope_phy = q(d$oracle_phy),
             target_s2_slope_ind = round((1 - lambda) * s2, 2), est_s2_slope_ind = q(d$s2_slope_ind),
             oracle_s2_slope_ind = q(d$oracle_ind),
             target_phylo_fraction = lambda, est_phylo_fraction = q(phy_frac),
             oracle_phylo_fraction = q(d$oracle_phy / (d$oracle_phy + d$oracle_ind)),
             est_b0 = q(d$b0), est_b1 = q(d$b1),
             reject_phylo = mean(d$p_phylo <= alpha, na.rm = TRUE),
             mean_prevalence = round(mean(d$mean_prevalence), 2),
             seconds_per_fit = round(mean(d$seconds), 1), stringsAsFactors = FALSE)
}
summary_table <- do.call(rbind, lapply(split(replicates, list(replicates$lambda, replicates$reduction)),
                                       summarise_scenario))
rownames(summary_table) <- NULL

# the site-intercept model in the scenarios with reduction
site_rows <- replicates[replicates$reduction < 1, ]
summary_site <- do.call(rbind, lapply(split(site_rows, site_rows$lambda), function(d) {
  q <- function(v) sprintf("%.2f (%.2f-%.2f)", median(v, na.rm = TRUE),
                           quantile(v, 0.1, na.rm = TRUE), quantile(v, 0.9, na.rm = TRUE))
  data.frame(lambda = d$lambda[1], n_error = sum(d$site_error == 1, na.rm = TRUE),
             n_not_converged = sum(!converged(d$site_convcode) & d$site_error == 0),
             est_s2_int = q(d$site_s2_int), est_s2_site = q(d$site_s2_site),
             est_s2_slope_phy = q(d$site_s2_slope_phy), est_s2_slope_ind = q(d$site_s2_slope_ind),
             est_phylo_fraction = q(d$site_s2_slope_phy / (d$site_s2_slope_phy + d$site_s2_slope_ind)),
             est_b0 = q(d$site_b0), reject_phylo = mean(d$site_p_phylo <= alpha, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
rownames(summary_site) <- NULL

# side check: the gradient of Appendix C, on [-15, 15], with lambda = 1 and no reduction
x_wide <- seq(-15, 15, length.out = n_sites)
set.seed(2012)
t0 <- Sys.time()
wide <- do.call(rbind, mclapply(seq_len(20), function(r) {
  sim <- simulate_logistic(tree, x_wide, lambda = 1, reduction = c(1, 1))
  d <- long_data(sim, x_wide)
  prevalence <- colMeans(sim$abundance)
  fit <- fit_pglmm(d)
  data.frame(rep = r, var_b = var(sim$traits$slope[, 1]),
             n_fixed_species = sum(prevalence %in% c(0, 1)),
             share_deterministic = mean(sim$probability < 0.01 | sim$probability > 0.99),
             t(fit))
}, mc.cores = cores, mc.set.seed = TRUE))
wide_minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

results <- list(
  design = list(n_sp = n_sp, n_sites = n_sites, n_rep = n_rep, tree = "balanced, Grafen",
                gradient = "evenly spaced, standardised", intercept_sd = 0.5, slope_sd = 1,
                formula = "y ~ x + (1 | sp) + (x | sp__)",
                formula_site = "y ~ x + (1 | sp) + (1 | site) + (x | sp__)",
                phylo_scale = phylo_scale,
                sample_factor = sapply(c(0, 0.5, 1), sample_factor)),
  replicates = replicates, summary = summary_table, summary_site = summary_site,
  wide_gradient = wide,
  timing = list(main_minutes = timing, wide_minutes = wide_minutes, cores = cores),
  session = list(phyr = as.character(packageVersion("phyr")),
                 comdynamite = as.character(packageVersion("comdynamite")),
                 R = R.version.string, date = Sys.Date())
)
saveRDS(results, "data-raw/validation_pglmm.rds")
print(summary_table)
print(summary_site)
print(summary(wide[, c("var_b", "n_fixed_species", "share_deterministic", "s2_slope_phy",
                       "s2_slope_ind", "convcode", "p_phylo")]))
