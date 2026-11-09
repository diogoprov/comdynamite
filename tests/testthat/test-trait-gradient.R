toy <- data.frame(site = c("A", "A", "B", "B", "B", "C"),
                  species = c("s1", "s2", "s1", "s3", "s2", "s3"),
                  trait = c(1, 3, 2, 6, 4, 8), abundance = c(1, 1, 2, 1, 1, 1))

tg_sim <- function(seed = 1, n_tips = 20, n_sites = 15, lambda = 1) {
  set.seed(seed)
  tree <- treeshapesim::sim_tree_shape(n_tips)
  sim <- simulate_metacommunity(landscape_implicit(n_sites), species_traits(tree, lambda = lambda),
                                competition_matrix(tree, "none"), n_steps = 100,
                                seed_steps = 1:5)
  list(tree = tree, ind = sample_individuals(sim, n = 3, local_adaptation = 0.5,
                                             intraspecific_sd = 1))
}

test_that("beta and alpha match the values computed by hand", {
  loo <- trait_gradient(toy)
  s <- loo$species
  expect_equal(s$beta, c(4, 2.5, 3))
  expect_equal(s$mean, c(1.5, 3.5, 7))
  expect_equal(s$alpha, s$mean - s$beta)
  expect_equal(s$n_alone, c(0L, 0L, 1L))
  expect_equal(s$breadth, c(2, 3, 0))
  expect_equal(loo$sites$mean, c(2, 4, 8))
  original <- trait_gradient(toy, leave_one_out = FALSE)$species
  expect_equal(original$beta, c(3, 3, 6))
  weighted <- trait_gradient(toy, abundance = "abundance", leave_one_out = FALSE)
  expect_equal(weighted$sites$mean, c(2, 3.5, 8))
  expect_equal(weighted$species$beta[1], 3)
  expect_equal(weighted$species$mean[1], 5 / 3)
  mixed <- trait_gradient(toy, abundance = "abundance", species_weights = "occurrence",
                          leave_one_out = FALSE)
  expect_equal(mixed$species$beta[1], (2 + 3.5) / 2)
})

test_that("individuals are averaged and unmeasured populations get the species mean", {
  ind <- rbind(toy, data.frame(site = "A", species = "s1", trait = 3, abundance = 1))
  s <- trait_gradient(ind)$species
  expect_equal(s$mean[1], (2 + 2) / 2)
  gap <- toy
  gap$trait[3] <- NA
  tg <- trait_gradient(gap)
  expect_equal(tg$populations$trait[3], 1)
  expect_false(tg$populations$measured[3])
  # slopes use measured values only and need min_slope sites
  expect_true(all(is.na(tg$species$slope)))
  line <- data.frame(site = rep(letters[1:4], each = 2), species = rep(c("x", "y"), 4),
                     trait = c(1, 0, 2, 1, 3, 2, 4, 3))
  fit <- trait_gradient(line, min_slope = 3)$species
  expect_equal(fit$slope, c(1, 1))
  expect_equal(fit$slope_se, c(0, 0))
  expect_true(is.na(trait_gradient(line, min_slope = 5)$species$slope[1]))
})

test_that("the input is checked", {
  expect_error(trait_gradient(as.matrix(toy)), "data frame")
  expect_error(trait_gradient(toy, trait = "z"), "no column z")
  expect_error(trait_gradient(transform(toy, trait = as.character(trait))), "numeric")
  expect_error(trait_gradient(toy, site_weights = "abundance"), "need the column")
  bad <- rbind(toy, data.frame(site = "A", species = "s1", trait = 3, abundance = 9))
  expect_error(trait_gradient(bad, abundance = "abundance"), "one value per population")
  expect_error(trait_gradient(transform(toy, abundance = -1), abundance = "abundance"),
               "0 or more")
  expect_error(trait_gradient(transform(toy, trait = NA_real_)), "no values")
  none <- toy
  none$trait[none$species == "s2"] <- NA
  expect_error(trait_gradient(none), "without any trait value: s2")
  expect_error(trait_gradient(toy[toy$site == "B", ]), "two sites")
  expect_error(trait_gradient(toy, leave_one_out = NA))
  expect_error(trait_gradient(toy, min_slope = 1), "`min_slope`")
  zero <- transform(toy, abundance = c(0, 1, 2, 1, 1, 1))
  expect_equal(trait_gradient(zero, abundance = "abundance")$species$n_sites, c(1L, 2L, 2L))
})

test_that("the original analysis is the leave_one_out = FALSE case", {
  d <- tg_sim()
  loo <- trait_gradient(d$ind, abundance = "abundance")
  orig <- trait_gradient(d$ind, abundance = "abundance", leave_one_out = FALSE)
  expect_equal(loo$species$mean, orig$species$mean)
  # with the focal species in the site means, the abundance-weighted mean of beta over the
  # populations equals that of the species means
  p <- orig$populations
  w <- p$weight
  b <- orig$species$beta[match(p$species, orig$species$species)]
  m <- orig$species$mean[match(p$species, orig$species$species)]
  expect_equal(sum(w * b), sum(w * m))
  expect_output(print(loo), "without the focal species")
  expect_output(print(orig), "with the focal species")
  expect_identical(as.data.frame(loo), loo$species)
})

test_that("trait_gradient_weights() compares the four choices", {
  d <- tg_sim()
  w <- trait_gradient_weights(d$ind, abundance = "abundance")
  expect_named(w, c("beta", "alpha", "correlation"))
  expect_equal(ncol(w$beta), 5)
  expect_equal(w$correlation$beta[1], 1)
  expect_equal(w$beta[["abundance/abundance"]],
               trait_gradient(d$ind, abundance = "abundance")$species$beta)
  expect_error(trait_gradient_weights(d$ind), "needs the column")
})

test_that("the curveball keeps the richness of sites and the frequency of species", {
  set.seed(2)
  d <- tg_sim()
  x <- tg_records(trait_gradient(d$ind, abundance = "abundance"))
  moved <- curveball_sites(x$site, x$species, 500)
  expect_equal(as.vector(table(moved)[sort(unique(x$site))]),
               as.vector(table(x$site)[sort(unique(x$site))]))
  expect_false(anyDuplicated(paste(moved, x$species)) > 0)
  expect_false(identical(moved, x$site))
  expect_identical(curveball_sites(x$site, x$species, 0), x$site)
  expect_identical(curveball_sites(rep("a", 3), c("p", "q", "r"), 5), rep("a", 3))
  # sites with the same species cannot trade
  expect_identical(curveball_sites(c("a", "b"), c("p", "p"), 5), c("a", "b"))
})

test_that("trait_gradient_test() returns the statistics and species tests", {
  d <- tg_sim()
  tg <- trait_gradient(d$ind, abundance = "abundance")
  set.seed(3)
  res <- trait_gradient_test(tg, n_null = 19, burn_in = 50, thin = 5)
  expect_s3_class(res, "mc_trait_gradient_test")
  expect_equal(res$statistics$statistic, c("cor_beta_alpha", "var_beta", "var_alpha"))
  expect_equal(nrow(res$species), nrow(tg$species))
  expect_equal(dim(res$null), c(19, 3 + 2 * nrow(tg$species)))
  expect_true(all(res$statistics$p > 0 & res$statistics$p <= 1))
  expect_output(print(res), "curveball")
  perm <- trait_gradient_test(tg, n_null = 9, null = "permute")
  expect_equal(perm$settings$null, "permute")
  expect_error(trait_gradient_test(list()), "come from trait_gradient")
  expect_error(trait_gradient_test(tg, n_null = 0), "`n_null`")
})

test_that("the GLS with known error matches nlme and reduces to BM", {
  skip_if_not_installed("nlme")
  d <- tg_sim(n_tips = 25, n_sites = 20)
  tg <- trait_gradient(d$ind, abundance = "abundance")
  set.seed(4)
  ph <- trait_gradient_phylo(tg, d$tree, n_boot = 29)
  dat <- ph$data
  parts <- tg_tree_parts(ph$tree, dat$species)
  bm <- gls_known_error(dat$beta, rep(0, nrow(dat)), parts, "BM")
  ref <- nlme::gls(beta ~ 1, data = dat, method = "REML",
                   correlation = ape::corBrownian(1, ph$tree, form = ~species))
  expect_equal(bm$loglik, as.numeric(stats::logLik(ref)), tolerance = 1e-6)
  expect_equal(bm$mu, unname(stats::coef(ref)), tolerance = 1e-6)
  white <- gls_known_error(dat$beta, rep(0, nrow(dat)), parts, "white")
  expect_equal(white$mu, mean(dat$beta), tolerance = 1e-6)
  expect_equal(white$sigma2, stats::var(dat$beta), tolerance = 1e-4)
  # error makes the phylogenetic variance smaller
  noisy <- gls_known_error(dat$beta, rep(stats::var(dat$beta) / 2, nrow(dat)), parts, "white")
  expect_lt(noisy$sigma2, white$sigma2)
  ou <- gls_known_error(dat$beta, rep(0, nrow(dat)), parts, "OU")
  expect_gte(ou$loglik, bm$loglik - 0.01)
})

test_that("trait_gradient_phylo() fits the models and checks its input", {
  d <- tg_sim(n_tips = 25, n_sites = 20)
  tg <- trait_gradient(d$ind, abundance = "abundance")
  set.seed(5)
  ph <- trait_gradient_phylo(tg, d$tree, responses = c("beta", "alpha"), n_boot = 19)
  expect_s3_class(ph, "mc_trait_gradient_phylo")
  expect_equal(nrow(ph$models), 6)
  expect_equal(as.vector(tapply(ph$models$delta_aic, ph$models$response, min)), c(0, 0))
  expect_true(all(ph$data$se_beta >= 0))
  expect_true(all(ph$data$n_sites >= 3))
  expect_setequal(ph$tree$tip.label, ph$data$species)
  expect_output(print(ph), "species")
  expect_error(trait_gradient_phylo(list(), d$tree), "come from trait_gradient")
  expect_error(trait_gradient_phylo(tg, "x"), "`tree`")
  expect_error(trait_gradient_phylo(tg, d$tree, min_sites = 1000), "fewer than 4 species")
  se <- tg_boot_se(trait_gradient(toy), 9)
  expect_true(is.na(se["s3", "beta"]))
})

test_that("plots of the trait-gradient analysis", {
  skip_if_not_installed("ggplot2")
  d <- tg_sim()
  tg <- trait_gradient(d$ind, abundance = "abundance")
  expect_ggplot(plot(tg))
  expect_ggplot(plot(trait_gradient(d$ind, leave_one_out = FALSE)))
  expect_ggplot(plot(tg, type = "beta_alpha", species_labels = TRUE))
  expect_ggplot(plot(trait_gradient_test(tg, n_null = 9, burn_in = 10, thin = 2)))
  expect_error(plot(tg, colour = c("red", "blue")))
})
