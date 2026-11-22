# rOpenSci statistical standards tested in this file (see the srr tags in R/)
#' @srrstats {BS2.13, BS7.0, BS7.1, BS7.4} Tested in this file, as described in
#'   R/srr-stats-bayesian.R.
NULL

# Extended tests of the analyses of trait variation among individuals: trait-gradient analysis
# (trait_gradient(), its null models, phylogenetic and Bayesian models), trait_flex(),
# nested_variance(), tstats() and tstats_summary(), community_index(), population_shifts() and
# trait_rarefaction(). They check parameter recovery over several seeds (G5.4, G5.6, G5.6b),
# convergence as data grow (G5.7), type I error and power of the permutation and bootstrap
# tests, and, for trait_gradient_bayes(), BS7.0-BS7.2. Run with COMDYNAMITE_EXTENDED_TESTS=true;
# with COMDYNAMITE_REPORT=true the observed numbers are printed as messages.

report_numbers <- function(...) {
  if (identical(Sys.getenv("COMDYNAMITE_REPORT"), "true")) message(...)
}

# sites along a gradient E in [0, 1]; each species occurs where E is near its optimum and has
# one population mean per site around its species mean; with `sorting` the species mean is
# its optimum (traits sort species along the gradient), otherwise independent of it
tg_gradient_data <- function(n_species = 20, n_sites = 20, sorting = FALSE, niche = 0.15,
                             within = 0.3) {
  env <- seq(0, 1, length.out = n_sites)
  optimum <- stats::runif(n_species)
  rows <- lapply(seq_len(n_species), function(i) {
    present <- which(stats::runif(n_sites) < exp(-(env - optimum[i])^2 / (2 * niche^2)))
    if (length(present) < 2) present <- order(abs(env - optimum[i]))[1:2]
    m <- if (sorting) 3 * optimum[i] else stats::rnorm(1)
    data.frame(site = paste0("s", present), species = paste0("sp", i),
               trait = m + stats::rnorm(length(present), 0, within))
  })
  do.call(rbind, rows)
}

# every species at every site with t_ij = a_i + b_i E_j + e_ij; with leave-one-out site means
# and occurrence weights p_j^(-i) = mean_{k != i}(a_k) + mean_{k != i}(b_k) E_j + error, so the
# slope of species i on the site means is b_i / mean_{k != i}(b_k)
tg_slope_data <- function(n_species = 30, n_sites = 20, sd = 0.3) {
  env <- seq(0, 2, length.out = n_sites)
  a <- stats::rnorm(n_species)
  b <- stats::runif(n_species, 0.5, 1.5)
  d <- expand.grid(site = seq_len(n_sites), species = seq_len(n_species))
  d$trait <- a[d$species] + b[d$species] * env[d$site] + stats::rnorm(nrow(d), 0, sd)
  d$site <- paste0("s", d$site)
  d$species <- sprintf("sp%02d", d$species)
  truth <- b / ((sum(b) - b) / (n_species - 1))
  list(data = d, truth = stats::setNames(truth, sprintf("sp%02d", seq_len(n_species))))
}

test_that("trait_gradient() recovers within-species slopes, better with more sites", {
  skip_if_not_extended()
  # G5.4, G5.6, G5.6b: slopes b_i recovered over 5 seeds; G5.7: error shrinks with sites.
  # Tolerance: with the focal species left out, the error of the site means is independent of
  # the response, and its variance (0.09 / 29) is under 1% of the variance of the true site
  # means (var(E) about 0.37), so attenuation is below 1% and the standardised errors
  # (slope - truth) / slope_se follow t distributions with n_sites - 2 df: mean 0 and, pooled
  # over 10 and 40 sites, sd sqrt(mean(c(8 / 6, 38 / 36))) = 1.09. Their mean over n = 300
  # species must be within 3 / sqrt(n) of 0 and their sd within 0.15 of 1.09 (about 3 times the
  # sd of a sample sd of n such values, near 0.05 with the tails of t with 8 df), and the RMSE
  # with 40 sites, near sqrt(10 / 40) = 0.5 of that with 10 sites, must be below 0.75 of it.
  rmse <- c()
  z <- c()
  for (n_sites in c(10, 40)) {
    err <- c()
    for (seed in 1:5) {
      set.seed(seed)
      d <- tg_slope_data(n_sites = n_sites)
      s <- trait_gradient(d$data)$species
      e <- s$slope - d$truth[s$species]
      err <- c(err, e)
      z <- c(z, e / s$slope_se)
    }
    rmse[as.character(n_sites)] <- sqrt(mean(err^2))
  }
  report_numbers("trait_gradient slopes: RMSE 10 sites ", signif(rmse[1], 3), ", 40 sites ",
                 signif(rmse[2], 3), "; z mean ", signif(mean(z), 3), ", z sd ",
                 signif(stats::sd(z), 3))
  expect_lt(abs(mean(z)), 3 / sqrt(length(z)))
  expect_lt(abs(stats::sd(z) - sqrt(mean(c(8 / 6, 38 / 36)))), 0.15)
  expect_lt(rmse[["40"]] / rmse[["10"]], 0.75)
})

test_that("trait_gradient_test() species null: type I error near 5%, high power", {
  skip_if_not_extended()
  # Type I error and power. 200 data sets per scenario, 199 null communities: the two-sided
  # p-value 2 (1 + r) / 200 is at most 0.05 for r <= 4 in either tail, an exact size of 0.05
  # under exchangeability. Under trait means independent of the niches, each statistic must be
  # compatible with 0.05 (exact binomial 99% interval); under sorting (trait mean equal to the
  # niche optimum), the variance of beta must reject in at least 90% of data sets (the
  # Validation section of the help reports 98-100%).
  n_rep <- 200
  run <- function(sorting) {
    t(vapply(seq_len(n_rep), function(r) {
      tg <- trait_gradient(tg_gradient_data(sorting = sorting))
      trait_gradient_test(tg, n_null = 199)$statistics$p <= 0.05
    }, logical(3)))
  }
  set.seed(101)
  null <- run(FALSE)
  set.seed(102)
  sorted <- run(TRUE)
  rates <- rbind(null = colMeans(null), sorting = colMeans(sorted))
  colnames(rates) <- c("cor_beta_alpha", "var_beta", "var_alpha")
  report_numbers("trait_gradient_test rejection rates:\n",
                 paste(utils::capture.output(print(rates)), collapse = "\n"))
  for (k in 1:3) expect_true(rejection_compatible(sum(null[, k]), n_rep))
  expect_gte(mean(sorted[, 2]), 0.9)
})

test_that("trait_gradient_phylo() recovers the Brownian rate of species means", {
  skip_if_not_extended()
  # G5.6, G5.6b: 10 seeds, 30 species on a tree of height 1, species means by Brownian motion
  # with rate 1 (ape::rTraitCont), each species at 12 of 25 sites with population sd 0.3, so
  # the sampling error of a mean (about 0.09) is small and enters the model as known error.
  # The REML estimate of sigma2 has a relative sd of about sqrt(2 / 29) = 0.26, so the mean of
  # 10 estimates must lie within 3 * 0.26 / sqrt(10) = 0.25 of 1. Power of model choice: BM
  # must beat the model without phylogeny by AIC in at least 8 of 10 data sets under Brownian
  # motion, and the reverse in at least 8 of 10 with independent species means.
  fit_one <- function(brownian) {
    tree <- treeshapesim::sim_tree_shape(30)
    tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
    m <- if (brownian) ape::rTraitCont(tree, sigma = 1) else
      stats::setNames(stats::rnorm(30), tree$tip.label)
    rows <- lapply(names(m), function(sp) {
      sites <- sample.int(25, 12)
      data.frame(site = paste0("s", sites), species = sp,
                 trait = m[[sp]] + stats::rnorm(12, 0, 0.3))
    })
    tg <- trait_gradient(do.call(rbind, rows))
    ph <- trait_gradient_phylo(tg, tree, responses = "mean", models = c("white", "BM"),
                               n_boot = 49)
    aic <- stats::setNames(ph$models$aic, ph$models$model)
    c(sigma2 = ph$models$sigma2[ph$models$model == "BM"],
      bm_better = unname(aic["BM"] < aic["white"]))
  }
  set.seed(201)
  bm <- t(replicate(10, fit_one(TRUE)))
  set.seed(202)
  white <- t(replicate(10, fit_one(FALSE)))
  report_numbers("trait_gradient_phylo: mean sigma2 ", signif(mean(bm[, "sigma2"]), 3),
                 " (range ", paste(signif(range(bm[, "sigma2"]), 3), collapse = "-"),
                 "); BM chosen ", sum(bm[, "bm_better"]), "/10 under BM, ",
                 sum(white[, "bm_better"]), "/10 under white noise")
  expect_lt(abs(mean(bm[, "sigma2"]) - 1), 3 * sqrt(2 / 29) / sqrt(10))
  expect_gte(sum(bm[, "bm_better"]), 8)
  expect_lte(sum(white[, "bm_better"]), 2)
})

# species means by Brownian motion on a tree of height 1, niche optima too; populations at
# the sites near the optimum, with sd `within` around the species mean
tgb_data <- function(n_species, n_sites, within = 0.3) {
  tree <- treeshapesim::sim_tree_shape(n_species)
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  optimum <- ape::rTraitCont(tree, sigma = 0.4)
  m <- 2 * ape::rTraitCont(tree, sigma = 1)
  env <- seq(min(optimum) - 0.1, max(optimum) + 0.1, length.out = n_sites)
  rows <- lapply(seq_len(n_species), function(i) {
    present <- which(stats::runif(n_sites) < exp(-(env - optimum[i])^2 / (2 * 0.3^2)))
    if (length(present) < 3) present <- order(abs(env - optimum[i]))[1:3]
    data.frame(site = paste0("s", present), species = tree$tip.label[i],
               trait = m[i] + stats::rnorm(length(present), 0, within))
  })
  list(tree = tree, mean = m, tg = trait_gradient(do.call(rbind, rows)))
}

tgb_fit <- function(d, seed, ...) {
  suppressMessages(trait_gradient_bayes(d$tg, d$tree, chains = 2, iter_warmup = 1000,
                                        iter_sampling = 1000, adapt_delta = 0.99, seed = seed,
                                        refresh = 0, show_messages = FALSE, ...))
}

test_that("trait_gradient_bayes() recovers species means and converges", {
  skip_if_not_extended()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  # BS7.0 (recovery of known parameters) and BS7.2 (the posterior moves away from the prior
  # with informative data), 2 seeds, 12 species, 15 sites, 2 chains of 1000 + 1000 draws.
  # Each fit must have R-hat < 1.01 and no divergent transitions. The 95% intervals of the 24
  # species means must contain the true means at least qbinom(0.01, 24, 0.95) times (the
  # lower 1% quantile of the number covered when coverage is nominal). BS7.2: the median
  # width of these intervals must be below a quarter of the sd of the trait, whereas the
  # prior on the standardised scale (exponential(1) among-species sd) puts the species means
  # within about +-2 sd of the trait.
  covered <- 0
  widths <- c()
  for (seed in 1:2) {
    set.seed(300 + seed)
    d <- tgb_data(12, 15)
    fit <- tgb_fit(d, seed)
    report_numbers("trait_gradient_bayes seed ", seed, ": max R-hat ",
                   signif(fit$diagnostics$max_rhat, 4), ", divergent ",
                   fit$diagnostics$divergent)
    expect_lt(fit$diagnostics$max_rhat, 1.01)
    expect_equal(fit$diagnostics$divergent, 0)
    s <- fit$species
    truth <- d$mean[s$species]
    covered <- covered + sum(s$mean_lower <= truth & truth <= s$mean_upper)
    widths <- c(widths, (s$mean_upper - s$mean_lower) / stats::sd(d$tg$populations$trait))
  }
  report_numbers("trait_gradient_bayes: species means covered ", covered, "/24, median ",
                 "relative width ", signif(stats::median(widths), 3))
  expect_gte(covered, stats::qbinom(0.01, 24, 0.95))
  expect_lt(stats::median(widths), 0.25)
})

test_that("trait_gradient_bayes() follows the prior of lambda when data say little", {
  skip_if_not_extended()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  # BS7.1: with 5 species the phylogenetic share lambda is barely identified, so its posterior
  # stays near the prior. Under beta(2, 8) (median 0.18, 90% interval 0.04-0.43) and beta(8, 2)
  # (median 0.82, 0.57-0.96) the posterior median of lambda of the mean must lie inside the
  # central 90% interval of its prior, and the two medians must differ by more than 0.3 (half
  # the distance between the prior medians). Both fits must have R-hat < 1.01 and no
  # divergent transitions.
  set.seed(311)
  d <- tgb_data(5, 6)
  medians <- c()
  for (shape in list(c(2, 8), c(8, 2))) {
    fit <- tgb_fit(d, 1, lambda_prior = shape)
    expect_lt(fit$diagnostics$max_rhat, 1.01)
    expect_equal(fit$diagnostics$divergent, 0)
    est <- fit$signal$estimate[fit$signal$quantity == "lambda_mean"]
    prior <- stats::qbeta(c(0.05, 0.95), shape[1], shape[2])
    expect_gt(est, prior[1])
    expect_lt(est, prior[2])
    medians <- c(medians, est)
    report_numbers("trait_gradient_bayes prior beta(", shape[1], ", ", shape[2], "): ",
                   "posterior median of lambda_mean ", signif(est, 3), ", max R-hat ",
                   signif(fit$diagnostics$max_rhat, 4), ", divergent ",
                   fit$diagnostics$divergent)
  }
  expect_gt(medians[2] - medians[1], 0.3)
})

# individuals at sites along env in [-1, 1]; five species with means 0, 2, ..., 8, chosen at
# each site with probabilities that shift with env when `turnover` > 0, and an intraspecific
# response `b` to env
flex_data <- function(b = 0, turnover = 0, n_sites = 20, n_ind = 30) {
  env <- seq(-1, 1, length.out = n_sites)
  mu <- c(0, 2, 4, 6, 8)
  rows <- lapply(seq_len(n_sites), function(j) {
    p <- exp(turnover * (mu - 4) / 4 * env[j])
    sp <- sample.int(5, n_ind, replace = TRUE, prob = p)
    data.frame(site = paste0("s", j), species = paste0("sp", sp),
               trait = mu[sp] + b * env[j] + stats::rnorm(n_ind))
  })
  list(ind = do.call(rbind, rows), env = data.frame(site = paste0("s", seq_len(n_sites)),
                                                    env = env))
}

flex_p <- function(d, component) {
  dec <- trait_flex(d$ind, formula = ~env, env = d$env)$decomposition
  dec$p[dec$component == component]
}

test_that("trait_flex() recovers the intraspecific response and its F tests have size 5%", {
  skip_if_not_extended()
  # G5.6, G5.6b: without turnover the species are spread evenly over the gradient, so the
  # difference between the specific and fixed means of a site is b (env_j - mean env) plus
  # sampling noise of variance about 1 / 30, and the intraspecific sum of squares explained
  # by env estimates b^2 SS(env) + noise; over 10 seeds with b = 1 the mean of
  # explained / SS(env) must be within 0.1 of 1 (the noise adds about 1 / (30 SS(env)) = 0.005
  # and the sd of the ratio is about 2 sqrt(var(noise) / SS(env)) = 0.06 per seed).
  # Type I error: without turnover and with b = 0, the F tests of the intraspecific and of the
  # turnover components on env must each reject at a rate compatible with 0.05 over 300 data
  # sets. Power: b = 0.5 with turnover must reject the intraspecific test in at least 90% of
  # 100 data sets.
  set.seed(401)
  ratio <- vapply(1:10, function(seed) {
    d <- flex_data(b = 1)
    dec <- trait_flex(d$ind, formula = ~env, env = d$env)$decomposition
    dec$explained[dec$component == "intraspecific"] / sum((d$env$env - mean(d$env$env))^2)
  }, numeric(1))
  set.seed(402)
  p_null <- replicate(300, {
    d <- flex_data(b = 0, turnover = 0)
    c(flex_p(d, "intraspecific"), flex_p(d, "turnover"))
  })
  set.seed(404)
  p_power <- replicate(100, flex_p(flex_data(b = 0.5, turnover = 2), "intraspecific"))
  report_numbers("trait_flex: mean ratio ", signif(mean(ratio), 3), "; type I intraspecific ",
                 mean(p_null[1, ] <= 0.05), ", turnover ", mean(p_null[2, ] <= 0.05),
                 "; power ", mean(p_power <= 0.05))
  expect_lt(abs(mean(ratio) - 1), 0.1)
  expect_true(rejection_compatible(sum(p_null[1, ] <= 0.05), 300))
  expect_true(rejection_compatible(sum(p_null[2, ] <= 0.05), 300))
  expect_gte(mean(p_power <= 0.05), 0.9)
})

test_that("trait_flex(): the intraspecific F test is conservative under turnover along env", {
  skip_if_not_extended()
  # Documented behaviour (see the details of ?trait_flex), not a size of 5%. The fixed means use
  # species means estimated from the same individuals, so with b = 0 and equal site sizes n
  # the intraspecific site values D are, given the table of counts C (sites x species, with
  # species totals n_i), Gaussian with mean 0 and covariance
  #   V = (I - C diag(1 / n_i) C') / n
  # in units of the variance among individuals: the directions of the columns of C are
  # removed. With turnover along env those columns follow env, so the expected explained sum
  # of squares e' V e (e: centred env scaled to length 1) falls below the expected residual
  # mean square (tr(V) - 1' V 1 / J - e' V e) / (J - 2), and the F test is conservative.
  # Over 300 data sets with turnover = 2:
  # (1) the ratio of the mean explained sum of squares to the mean residual mean square must
  #     agree with that of the expectations above (about 0.6) within 25%: 3 Monte Carlo
  #     standard errors, the relative error of a mean of 300 squared Gaussians being
  #     sqrt(2 / 300) = 0.082 and that of the residual mean square sqrt(2 / (18 * 300)) = 0.019;
  #     the expected ratio must be below 0.75, so that such an estimate stays below 1;
  # (2) the test must not reject more often than the level: at most qbinom(0.99, 300, 0.05)
  #     rejections (the rate was 1.6% of 3000 data sets when the help was written);
  # (3) with the true species means in the fixed mean the cause is removed, and the rate of
  #     the same F test must be compatible with 0.05.
  set.seed(403)
  mu <- c(0, 2, 4, 6, 8)
  sims <- replicate(300, {
    d <- flex_data(b = 0, turnover = 2)
    res <- trait_flex(d$ind, formula = ~env, env = d$env)
    dec <- res$decomposition[res$decomposition$component == "intraspecific", ]
    site <- factor(d$ind$site, d$env$site)
    counts <- unclass(table(site, d$ind$species))
    n_sites <- nrow(counts)
    v <- (diag(n_sites) - counts %*% (t(counts) / colSums(counts)) / 30) / 30
    e <- d$env$env - mean(d$env$env)
    e <- e / sqrt(sum(e^2))
    along <- drop(t(e) %*% v %*% e)
    true_fixed <- tapply(mu[as.integer(sub("sp", "", d$ind$species))], site, mean)
    oracle <- res$sites$specific - as.vector(true_fixed)
    c(p = dec$p, explained = dec$explained, mean_square = dec$residual / (n_sites - 2),
      expected_explained = along,
      expected_mean_square = (sum(diag(v)) - sum(v) / n_sites - along) / (n_sites - 2),
      p_oracle = summary(stats::lm(oracle ~ d$env$env))$coefficients[2, 4])
  })
  m <- rowMeans(sims)
  observed <- unname(m["explained"] / m["mean_square"])
  expected <- unname(m["expected_explained"] / m["expected_mean_square"])
  report_numbers("trait_flex under turnover: type I ", mean(sims["p", ] <= 0.05),
                 "; explained / residual mean square ", signif(observed, 3), ", expected ",
                 signif(expected, 3), "; type I with true species means ",
                 mean(sims["p_oracle", ] <= 0.05))
  expect_lt(expected, 0.75)
  expect_lt(abs(observed / expected - 1), 0.25)
  expect_lte(sum(sims["p", ] <= 0.05), stats::qbinom(0.99, 300, 0.05))
  expect_true(rejection_compatible(sum(sims["p_oracle", ] <= 0.05), 300))
})

test_that("nested_variance() recovers nested variance components, better with more sites", {
  skip_if_not_extended()
  skip_if_not_installed("nlme")
  # G5.6, G5.6b: region / site / individual with variances 2, 1 and 1, 10 regions, 6 sites per
  # region and 10 individuals per site, 20 seeds. REML estimates of variance components are
  # nearly unbiased when the true values are well above 0, so the mean of each estimate over
  # the seeds must lie within 3 Monte Carlo standard errors (sd / sqrt(20)) of the truth.
  # G5.7: the sd of the site variance with 24 sites per region must be below 0.75 of that with
  # 3 sites per region (degrees of freedom 230 against 20 give a ratio near 0.3).
  sim_nested <- function(n_site) {
    d <- expand.grid(ind = 1:10, site = seq_len(n_site), region = 1:10)
    d$site <- paste(d$region, d$site)
    d$trait <- stats::rnorm(10, 0, sqrt(2))[d$region] +
      stats::rnorm(10 * n_site)[match(d$site, unique(d$site))] + stats::rnorm(nrow(d))
    attr(nested_variance(d, "trait", c("region", "site")), "variance")[, "trait"]
  }
  set.seed(501)
  est <- t(replicate(20, sim_nested(6)))
  truth <- c(region = 2, site = 1, within = 1)
  mc_se <- apply(est, 2, stats::sd) / sqrt(20)
  set.seed(502)
  few <- replicate(20, sim_nested(3))["site", ]
  set.seed(503)
  many <- replicate(20, sim_nested(24))["site", ]
  report_numbers("nested_variance: mean estimates ",
                 paste(names(truth), signif(colMeans(est), 3), collapse = ", "),
                 "; sd of site variance 3 sites ", signif(stats::sd(few), 3), ", 24 sites ",
                 signif(stats::sd(many), 3))
  expect_true(all(abs(colMeans(est)[names(truth)] - truth) < 3 * mc_se[names(truth)]))
  expect_lt(stats::sd(many) / stats::sd(few), 0.75)
})

test_that("nested_variance() recovers crossed variance components", {
  skip_if_not_extended()
  skip_if_not_installed("lme4")
  # G5.6, G5.6b: 12 species crossed with 15 sites, 4 individuals per population, variances
  # species 3, site 1, population (site:species) 0.5 and within 1; 15 seeds. As above, the
  # mean of each estimate must lie within 3 Monte Carlo standard errors of the truth.
  sim_crossed <- function() {
    d <- expand.grid(ind = 1:4, species = sprintf("sp%02d", 1:12),
                     site = sprintf("s%02d", 1:15), stringsAsFactors = FALSE)
    pop <- interaction(d$site, d$species, drop = TRUE)
    d$trait <- stats::rnorm(12, 0, sqrt(3))[match(d$species, unique(d$species))] +
      stats::rnorm(15)[match(d$site, unique(d$site))] +
      stats::rnorm(nlevels(pop), 0, sqrt(0.5))[pop] + stats::rnorm(nrow(d))
    v <- attr(suppressWarnings(nested_variance(d, "trait", "site", crossed = "species")),
              "variance")
    stats::setNames(v[, "trait"], rownames(v))
  }
  set.seed(511)
  est <- t(replicate(15, sim_crossed()))
  means <- colMeans(est)
  mc_se <- apply(est, 2, stats::sd) / sqrt(15)
  truth <- c(species = 3, site = 1, within = 1)
  pop_level <- setdiff(colnames(est), names(truth))
  report_numbers("nested_variance crossed: mean estimates ",
                 paste(names(means), signif(means, 3), collapse = ", "))
  expect_length(pop_level, 1)
  expect_true(all(abs(means[names(truth)] - truth) < 3 * mc_se[names(truth)]))
  expect_lt(abs(means[pop_level] - 0.5), 3 * mc_se[pop_level])
})

# eight sites of 40 individuals of six species: under the null every species has the same
# trait distribution and individuals are placed at random; with filtering, species have
# distinct means and each site holds the species whose means are near its environment
tstats_data <- function(filtering = FALSE) {
  if (!filtering) {
    return(data.frame(site = paste0("s", rep(1:8, each = 40)),
                      species = paste0("sp", sample.int(6, 320, replace = TRUE)),
                      trait = stats::rnorm(320)))
  }
  mu <- seq(0, 10, length.out = 6)
  env <- seq(0, 10, length.out = 8)
  rows <- lapply(1:8, function(j) {
    sp <- sample.int(6, 40, replace = TRUE, prob = exp(-(mu - env[j])^2 / (2 * 2^2)))
    data.frame(site = paste0("s", j), species = paste0("sp", sp),
               trait = mu[sp] + stats::rnorm(40, 0, 0.5))
  })
  do.call(rbind, rows)
}

# p-values of site s1 for the three statistics and of the mean tests of tstats_summary() at
# level = 0.9 (one-sided 0.05 in each direction), for data sets from tstats_data()
tstats_rejections <- function(n_rep, filtering) {
  ratios <- c("T_IP.IC", "T_IC.IR", "T_PC.PR")
  out <- t(vapply(seq_len(n_rep), function(r) {
    ts <- tstats(tstats_data(filtering), n_perm = 199)
    s1 <- unlist(ts[ts$site == "s1", paste0("p_", ratios)])
    tr <- tstats_summary(ts, level = 0.9)$traits
    c(s1 <= 0.05, tr$p_mean_low[match(ratios, tr$statistic)] <= 0.05)
  }, logical(6)))
  colnames(out) <- c(paste0("site_", ratios), paste0("mean_", ratios))
  out
}

test_that("tstats() and tstats_summary(): type I error near 5%, high power under filtering", {
  skip_if_not_extended()
  # Type I error and power of the three null models ("less", the default) and of the test of
  # the mean over sites in tstats_summary(). 200 data sets per scenario with 199
  # randomisations: p = (1 + r) / 200 <= 0.05 for r <= 9, an exact size of 0.05. One site per
  # data set (s1) is used so that the replicates are independent. Under the null scenario the
  # individuals of a site are exchangeable with those of the region and among species, which
  # is the null of each randomisation; the rejection rates of the three site tests and of the
  # mean test of T_IP.IC must be compatible with 0.05 (the mean tests of T_IC.IR and T_PC.PR
  # are in the next test). Under filtering, the site and mean tests of T_IP.IC and T_IC.IR
  # must reject in at least 90%.
  set.seed(601)
  null <- tstats_rejections(200, FALSE)
  set.seed(602)
  filt <- tstats_rejections(200, TRUE)
  rates <- rbind(null = colMeans(null), filtering = colMeans(filt))
  report_numbers("tstats rejection rates:\n",
                 paste(utils::capture.output(print(rates)), collapse = "\n"))
  for (k in c("site_T_IP.IC", "site_T_IC.IR", "site_T_PC.PR", "mean_T_IP.IC")) {
    expect_true(rejection_compatible(sum(null[, k]), 200), label = k)
  }
  for (k in c("site_T_IP.IC", "site_T_IC.IR", "mean_T_IP.IC", "mean_T_IC.IR")) {
    expect_gte(mean(filt[, k]), 0.9, label = k)
  }
})

test_that("tstats_summary() mean tests of T_IC.IR and T_PC.PR: type I error near 5%", {
  skip_if_not_extended()
  # Type I error of the mean tests of T_IC.IR and T_PC.PR under the null scenario (600 data
  # sets, compatible with 0.05). The sites of a data set partition the region, so each
  # randomisation must permute the whole region for the mean over sites to have the right
  # null distribution. When each site drew from the region on its own (up to 0.10.0, and in
  # cati), the mean test at 5% rejected in 0/600 data sets for T_IC.IR and in 47/600 = 7.8%
  # (99% CI 5.3-11.1%) for T_PC.PR.
  set.seed(603)
  null <- tstats_rejections(600, FALSE)
  report_numbers("tstats mean tests under the null: T_IC.IR ", mean(null[, "mean_T_IC.IR"]),
                 ", T_PC.PR ", mean(null[, "mean_T_PC.PR"]))
  expect_true(rejection_compatible(sum(null[, "mean_T_IC.IR"]), 600))
  expect_true(rejection_compatible(sum(null[, "mean_T_PC.PR"]), 600))
})

# eight sites of 60 individuals of six species out of 20 with distinct means (0 to 10,
# individuals with sd 0.5). `filtering`: each site holds the six species nearest its
# environment, otherwise six species at random (no filtering of species). `sizes`: individuals
# of the six populations of each site, in random order among its species; a list gives the
# sizes of sites 1-4 and of sites 5-8
species_filter_data <- function(filtering, sizes) {
  mu <- seq(0, 10, length.out = 20)
  env <- seq(0, 10, length.out = 8)
  if (!is.list(sizes)) sizes <- list(sizes, sizes)
  rows <- lapply(1:8, function(j) {
    pool <- if (filtering) order(abs(mu - env[j]))[1:6] else sample.int(20, 6)
    sp <- rep(pool[sample.int(6)], sizes[[1 + (j > 4)]])
    data.frame(site = paste0("s", j), species = paste0("sp", sp),
               trait = mu[sp] + stats::rnorm(length(sp), 0, 0.5))
  })
  do.call(rbind, rows)
}

test_that("tstats(): power of T_PC.PR under filtering of species, both null models", {
  skip_if_not_extended()
  # T_PC.PR compares the variance among the population means of a site with that of the
  # region, so it has power when the species of a site are a narrow part of those of the
  # region and the site has enough populations. (In the scenario of tstats_data() sites hold
  # three to six of six species and the variance among three unweighted means is too variable
  # for a test of one site: 0/200 rejections with either null model.) Here each site holds
  # the six species out of 20 nearest its environment, with uneven populations (30, 15, 8, 4,
  # 2 and 1 individuals): the test of site s1 and the test of the mean must reject in at least
  # 90% of 200 data sets with both null models.
  uneven <- c(30, 15, 8, 4, 2, 1)
  set.seed(604)
  power <- t(vapply(1:200, function(r) {
    d <- species_filter_data(TRUE, uneven)
    unlist(lapply(c("individuals", "populations"), function(null_pc) {
      ts <- tstats(d, n_perm = 199, null_pc = null_pc)
      tr <- tstats_summary(ts, level = 0.9)$traits
      c(ts$p_T_PC.PR[1], tr$p_mean_low[tr$statistic == "T_PC.PR"]) <= 0.05
    }))
  }, logical(4)))
  colnames(power) <- c("site_individuals", "mean_individuals", "site_populations",
                       "mean_populations")
  report_numbers("tstats T_PC.PR power under species filtering: ",
                 paste(colnames(power), colMeans(power), collapse = ", "))
  for (k in colnames(power)) expect_gte(mean(power[, k]), 0.9, label = k)
})

test_that("tstats(): the T_PC.PR null of individuals depends on population sizes", {
  skip_if_not_extended()
  # Documented behaviour (section "Null model of T_PC.PR and the size of populations" of
  # ?tstats). Species are placed at random among sites, so no site is filtered, but sites
  # 1-4 have six populations of 10 individuals and sites 5-8 populations of 45, 10, 2, 1, 1
  # and 1. With null_pc = "individuals" (regional.pop of Taudiere & Violle 2016) the null
  # population means of large populations stay near the regional mean and those of small
  # ones spread, so the test for high values at s1 rejects far more often than 5% (84.5% in
  # 200 data sets when the test was written; at least 50% is required). With null_pc =
  # "populations" the sizes play no part: the test for high values at s1 and that for low
  # values at s8 must not reject more often than 5% (lower limit of the 99% interval of the
  # rate at or below 0.05; the test is conservative, because populations of one species at
  # several sites are permuted as independent values).
  sizes <- list(rep(10, 6), c(45, 10, 2, 1, 1, 1))
  set.seed(605)
  rejected <- t(vapply(1:200, function(r) {
    d <- species_filter_data(FALSE, sizes)
    high <- function(null_pc) {
      tstats(d, n_perm = 199, alternative = "greater", null_pc = null_pc)$p_T_PC.PR[1]
    }
    low <- tstats(d, n_perm = 199, null_pc = "populations")$p_T_PC.PR[8]
    c(high("individuals"), high("populations"), low) <= 0.05
  }, logical(3)))
  colnames(rejected) <- c("high_individuals", "high_populations", "low_populations")
  report_numbers("tstats T_PC.PR without filtering, uneven sizes: ",
                 paste(colnames(rejected), colMeans(rejected), collapse = ", "))
  expect_gte(mean(rejected[, "high_individuals"]), 0.5)
  for (k in c("high_populations", "low_populations")) {
    lower <- stats::binom.test(sum(rejected[, k]), 200, conf.level = 0.99)$conf.int[1]
    expect_lte(lower, 0.05, label = k)
  }
})

test_that("tstats() with null_ip = 'species': type I error near 5%, power", {
  skip_if_not_extended()
  # Type I error of the null model of T_IP.IC that exchanges population means within species,
  # two-sided, 200 data sets, 199 randomisations. Five species with the same mean are at every
  # site with 6 individuals around a population mean that is the species mean plus an
  # independent deviation (sd 1), so population means are exchangeable within species. The
  # help documents the test as conservative when populations differ in size; with equal sizes
  # the rate must be compatible with 0.05. Power: under local adaptation (all populations of a
  # site share its deviation, sd 1), the test must reject in at least 80% of 50 data sets.
  sim <- function(adapted) {
    d <- expand.grid(ind = 1:6, species = 1:5, site = 1:10)
    env <- stats::rnorm(10)
    dev <- if (adapted) env[d$site] else stats::rnorm(50)[(d$site - 1) * 5 + d$species]
    d$trait <- dev + stats::rnorm(nrow(d), 0, 0.3)
    d$site <- paste0("s", d$site)
    d$species <- paste0("sp", d$species)
    ts <- tstats(d, n_perm = 199, alternative = "two.sided", null_ip = "species")
    ts$p_T_IP.IC[ts$site == "s1"] <= 0.05
  }
  set.seed(611)
  null <- replicate(200, sim(FALSE))
  set.seed(612)
  adapted <- replicate(50, sim(TRUE))
  report_numbers("tstats null_ip = species: type I ", mean(null), ", power ", mean(adapted))
  expect_true(rejection_compatible(sum(null), 200))
  expect_gte(mean(adapted), 0.8)
})

test_that("community_index(): type I error near 5%, high power under filtering", {
  skip_if_not_extended()
  # Type I error and power, two-sided (the default), 200 data sets, 199 randomisations:
  # p = 2 (1 + r) / 200 <= 0.05 for r <= 4 in either tail, an exact size of 0.05. Site s1 of
  # each data set. "regional.ind" with the variance of the site, and "local" with the variance
  # of the population means of the site, under the null scenario of tstats_data(), whose
  # individuals are exchangeable within sites and with the region; rates must be compatible
  # with 0.05. Under filtering the variance of s1 against "regional.ind" must reject in at
  # least 90%.
  pop_var <- function(x, species) stats::var(tapply(x, species, mean))
  p_s1 <- function(d, index, null) {
    res <- community_index(d, index, null = null, n_perm = 199)
    res$p[res$site == "s1"] <= 0.05
  }
  set.seed(701)
  regional <- replicate(200, p_s1(tstats_data(), stats::var, "regional.ind"))
  set.seed(702)
  local <- replicate(200, p_s1(tstats_data(), pop_var, "local"))
  set.seed(703)
  power <- replicate(100, p_s1(tstats_data(TRUE), stats::var, "regional.ind"))
  report_numbers("community_index: type I regional.ind ", mean(regional), ", local ",
                 mean(local), "; power ", mean(power))
  expect_true(rejection_compatible(sum(regional), 200))
  expect_true(rejection_compatible(sum(local), 200))
  expect_gte(mean(power), 0.9)
})

# individuals of a simulated metacommunity (25 sites, 20 species), with known local adaptation
shifts_data <- function(n = 6, local_adaptation = 0.3, displacement = 0) {
  tree <- treeshapesim::sim_tree_shape(20)
  sim <- simulate_metacommunity(landscape_implicit(25), species_traits(tree),
                                competition_matrix(tree, "none"), n_steps = 100,
                                seed_steps = 1:5)
  sample_individuals(sim, n = n, local_adaptation = local_adaptation, intraspecific_sd = 1,
                     displacement = displacement, displacement_width = 5)
}

test_that("population_shifts() recovers local adaptation, better with more individuals", {
  skip_if_not_extended()
  # G5.6, G5.6b: h = 0.3 from sample_individuals() over 10 seeds (no displacement, width fixed
  # at 5 trait units); the mean estimate must lie within 3 Monte Carlo standard errors of 0.3.
  # G5.7: the error of a population mean falls with the number of individuals, so the RMSE of
  # h over 25 seeds with 20 individuals per population must be below 0.75 of that with 2. If
  # only the sampling error of the means mattered the ratio would be sqrt(2 / 20) = 0.32; it
  # is larger because the displacement c, estimated jointly, carries an error that shrinks
  # more slowly (in pilots the correlation of the errors of h and c was 0.75 with n = 20).
  h_of <- function(n) population_shifts(shifts_data(n = n), width = 5)$local_adaptation
  set.seed(801)
  h6 <- replicate(10, h_of(6))
  set.seed(802)
  h2 <- replicate(25, h_of(2))
  set.seed(803)
  h20 <- replicate(25, h_of(20))
  rmse <- function(h) sqrt(mean((h - 0.3)^2))
  report_numbers("population_shifts: mean h ", signif(mean(h6), 4), " (sd ",
                 signif(stats::sd(h6), 3), "); RMSE n = 2 ", signif(rmse(h2), 3), ", n = 20 ",
                 signif(rmse(h20), 3))
  expect_lt(abs(mean(h6) - 0.3), 3 * stats::sd(h6) / sqrt(10))
  expect_lt(rmse(h20) / rmse(h2), 0.75)
})

test_that("population_shifts() bootstrap test of displacement: type I error and power", {
  skip_if_not_extended()
  # Type I error and power of p_displacement, 120 data sets from sample_individuals() without
  # displacement (h = 0.3), 99 bootstrap simulations each, width fixed at 5. The help reports
  # 2-3% rejections at 5% (conservative); the rate must be compatible with 0.05 or below it
  # (not inflated: the lower end of the 99% interval at most 0.05). Power: displacement = 2
  # (displacement_width = 5) must be detected in at least 80% of 40 data sets.
  p_of <- function(c) {
    population_shifts(shifts_data(displacement = c), width = 5, nboot = 99)$p_displacement
  }
  set.seed(811)
  p_null <- replicate(120, p_of(0))
  set.seed(812)
  p_alt <- replicate(40, p_of(2))
  rejections <- sum(p_null <= 0.05)
  report_numbers("population_shifts p_displacement: type I ", signif(rejections / 120, 3),
                 ", power ", mean(p_alt <= 0.05))
  lower <- stats::binom.test(rejections, 120, conf.level = 0.99)$conf.int[1]
  expect_lte(lower, 0.05)
  expect_gte(mean(p_alt <= 0.05), 0.8)
})

test_that("trait_rarefaction() converges to the full-data statistic as effort grows", {
  skip_if_not_extended()
  # G5.7: with 30 individuals per population, subsamples of n individuals vary less as n
  # grows, and at n = 30 every subsample is the full data. Over 3 seeds the width of the
  # interval of each T-statistic must shrink from n = 2 to n = 20, and at n = 30 the interval
  # must have width 0 and equal the statistic of the full data (to 1e-10).
  for (seed in 1:3) {
    set.seed(900 + seed)
    d <- expand.grid(ind = 1:30, species = paste0("sp", 1:4), site = paste0("s", 1:6),
                     stringsAsFactors = FALSE)
    d$trait <- match(d$species, unique(d$species)) +
      stats::rnorm(24)[match(paste(d$site, d$species), unique(paste(d$site, d$species)))] +
      stats::rnorm(nrow(d))
    rare <- trait_rarefaction(d, n = c(2, 20, 30), reps = 30)
    width <- rare$high - rare$low
    for (s in unique(rare$statistic)) {
      w <- width[rare$statistic == s]
      expect_lt(w[rare$effort[rare$statistic == s] == 20],
                w[rare$effort[rare$statistic == s] == 2])
    }
    full <- tstats(d, n_perm = 0)
    full <- colMeans(full[c("T_IP.IC", "T_IC.IR", "T_PC.PR")], na.rm = TRUE)
    at30 <- rare[rare$effort == 30, ]
    expect_equal(at30$high - at30$low, rep(0, nrow(at30)), tolerance = 1e-10)
    expect_equal(at30$mean, unname(full[at30$statistic]), tolerance = 1e-10)
  }
})

test_that("trait_gradient_bayes() can run without printing anything", {
  skip_if_not_extended()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "CmdStan is not installed")
  # BS2.13: with refresh = 0, show_messages = FALSE and show_exceptions = FALSE passed through
  # `...` to $sample(), a fit prints nothing to the console and emits no messages. The only
  # condition left is the package's own warning when the diagnostics fail, which is allowed
  # here because short chains may trigger it. A first fit compiles the model and shows the
  # one-time citation message, so that the second fit checks the sampler's own output only.
  set.seed(401)
  d <- tgb_data(8, 10)
  quiet <- function() {
    trait_gradient_bayes(d$tg, d$tree, chains = 1, iter_warmup = 300, iter_sampling = 300,
                         seed = 1, refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  }
  suppressWarnings(suppressMessages(quiet()))
  messages <- character()
  warnings <- character()
  printed <- utils::capture.output(fit <- withCallingHandlers(
    quiet(),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ))
  expect_length(printed, 0)
  expect_length(messages, 0)
  expect_length(grep("has not converged", warnings, invert = TRUE, value = TRUE), 0)
  expect_s3_class(fit, "mc_trait_gradient_bayes")
})
