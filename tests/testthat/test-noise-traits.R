# Quick checks of sensitivity to noise (G5.9a) and to the seed (G5.9b) of the analyses of trait
# variation among individuals. G5.9a: a relative jitter of 1e-10 on the trait values (and on
# abundances), with the same seed, leaves the estimates equal up to 1e-6 (relative) and
# permutation p-values equal up to one step of the permutation distribution, 1 / (n_perm + 1),
# since a null value tied with the observed one can fall on either side. G5.9b: with different
# seeds, a one-sided permutation p-value has Monte Carlo standard error about
# sqrt(p (1 - p) / n_perm), and a two-sided one, twice a one-sided tail, about
# sqrt(p (2 - p) / n_perm); p-values must lie within 4 such errors (plus one step) of their
# mean, a bound exceeded with probability below 1e-4 per p-value under a normal approximation.

traits_noise_data <- function() {
  set.seed(41)
  d <- expand.grid(ind = 1:5, species = paste0("sp", 1:6), site = paste0("s", 1:8),
                   stringsAsFactors = FALSE)
  d <- d[stats::runif(nrow(d)) < 0.7, ]
  env <- stats::setNames(seq(-1, 1, length.out = 8), paste0("s", 1:8))
  pop <- paste(d$site, d$species)
  d$trait <- match(d$species, paste0("sp", 1:6)) + 0.5 * env[d$site] +
    stats::rnorm(length(unique(pop)), 0, 0.3)[match(pop, unique(pop))] + stats::rnorm(nrow(d))
  d$abundance <- stats::rpois(length(unique(pop)), 10)[match(pop, unique(pop))] + 1
  d$region <- ifelse(d$site %in% paste0("s", 1:4), "north", "south")
  list(ind = d, env = data.frame(site = names(env), env = unname(env)))
}

# jitter of each value of `columns`; abundances (one value per population) get one factor per
# population. The seed is set here, so call it before setting the seed of an analysis
traits_jitter <- function(d, columns = "trait") {
  set.seed(42)
  pop <- paste(d$site, d$species)
  for (k in columns) {
    noise <- stats::rnorm(nrow(d))
    if (k == "abundance") noise <- noise[match(pop, pop)]
    d[[k]] <- d[[k]] * (1 + 1e-10 * noise)
  }
  d
}

expect_traits_seed_stable <- function(p, n_perm, two_sided = FALSE) {
  p <- as.matrix(p)
  centre <- rowMeans(p)
  var_p <- if (two_sided) centre * (2 - centre) else centre * (1 - centre)
  bound <- 4 * sqrt(pmax(var_p, 1 / n_perm) / n_perm) + 1 / (n_perm + 1)
  expect_true(all(abs(p - centre) <= bound, na.rm = TRUE))
}

test_that("trait_gradient() and trait_gradient_weights() ignore trivial noise", {
  # G5.9a: jitter on the trait values and on the abundances
  d <- traits_noise_data()$ind
  j <- traits_jitter(d, c("trait", "abundance"))
  a <- trait_gradient(d, abundance = "abundance")$species
  b <- trait_gradient(j, abundance = "abundance")$species
  num <- vapply(a, is.numeric, logical(1))
  expect_equal(b[num], a[num], tolerance = 1e-6)
  wa <- trait_gradient_weights(d, abundance = "abundance")
  wb <- trait_gradient_weights(j, abundance = "abundance")
  expect_equal(wb$beta[-1], wa$beta[-1], tolerance = 1e-6)
  expect_equal(wb$alpha[-1], wa$alpha[-1], tolerance = 1e-6)
})

test_that("trait_gradient_test() ignores trivial noise and is stable across seeds", {
  # G5.9a with the same seed; G5.9b over 3 seeds, two-sided p-values with 299 null communities
  d <- traits_noise_data()$ind
  tg <- trait_gradient(d)
  tg_j <- trait_gradient(traits_jitter(d))
  set.seed(5)
  a <- trait_gradient_test(tg, n_null = 99)$statistics
  set.seed(5)
  b <- trait_gradient_test(tg_j, n_null = 99)$statistics
  expect_equal(b$observed, a$observed, tolerance = 1e-6)
  expect_true(all(abs(b$p - a$p) <= 2 / 100 + 1e-12)) # two-sided: one step is 2 / 100
  p <- vapply(1:3, function(seed) {
    set.seed(seed)
    trait_gradient_test(tg, n_null = 299)$statistics$p
  }, numeric(3))
  expect_traits_seed_stable(p, 299, two_sided = TRUE)
})

test_that("trait_flex() and nested_variance() ignore trivial noise", {
  # G5.9a: sums of squares, p-values and variance components
  d <- traits_noise_data()
  j <- traits_jitter(d$ind)
  a <- trait_flex(d$ind, formula = ~env, env = d$env)$decomposition
  b <- trait_flex(j, formula = ~env, env = d$env)$decomposition
  num <- vapply(a, is.numeric, logical(1))
  expect_equal(b[num], a[num], tolerance = 1e-6)
  skip_if_not_installed("nlme")
  va <- nested_variance(d$ind, "trait", c("region", "site"))
  vb <- nested_variance(j, "trait", c("region", "site"))
  expect_equal(vb$trait, va$trait, tolerance = 1e-5)
})

test_that("tstats() and community_index() ignore trivial noise and are stable across seeds", {
  # G5.9a with the same seed; G5.9b over 3 seeds with 199 randomisations ("less" for tstats(),
  # two-sided for community_index())
  d <- traits_noise_data()$ind
  j <- traits_jitter(d)
  ratios <- c("T_IP.IC", "T_IC.IR", "T_PC.PR")
  set.seed(5)
  a <- tstats(d, n_perm = 99)
  set.seed(5)
  b <- tstats(j, n_perm = 99)
  expect_equal(b[ratios], a[ratios], tolerance = 1e-6)
  pa <- as.matrix(a[paste0("p_", ratios)])
  pb <- as.matrix(b[paste0("p_", ratios)])
  expect_true(all(abs(pb - pa) <= 1 / 100 + 1e-12, na.rm = TRUE))
  set.seed(6)
  ca <- community_index(d, stats::var, null = "regional.ind", n_perm = 99)
  set.seed(6)
  cb <- community_index(j, stats::var, null = "regional.ind", n_perm = 99)
  expect_equal(cb$observed, ca$observed, tolerance = 1e-6)
  expect_true(all(abs(cb$p - ca$p) <= 2 / 100 + 1e-12))
  p_t <- vapply(1:3, function(seed) {
    set.seed(seed)
    as.vector(as.matrix(tstats(d, n_perm = 199)[paste0("p_", ratios)]))
  }, numeric(3 * 8))
  expect_traits_seed_stable(p_t, 199)
  p_c <- vapply(1:3, function(seed) {
    set.seed(seed)
    community_index(d, stats::var, null = "regional.ind", n_perm = 199)$p
  }, numeric(8))
  expect_traits_seed_stable(p_c, 199, two_sided = TRUE)
})

test_that("population_shifts() ignores trivial noise", {
  # G5.9a: with the width fixed and the same seed (the split of each population in halves is
  # random), the estimates do not change
  d <- traits_noise_data()
  ind <- merge(d$ind, d$env, by = "site")
  ind_j <- traits_jitter(ind)
  set.seed(7)
  a <- population_shifts(ind, width = 1)
  set.seed(7)
  b <- population_shifts(ind_j, width = 1)
  expect_equal(b$local_adaptation, a$local_adaptation, tolerance = 1e-6)
  expect_equal(b$displacement, a$displacement, tolerance = 1e-6)
})
