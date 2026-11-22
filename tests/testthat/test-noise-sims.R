# Cheap checks of sensitivity to noise (G5.9a) and to the seed (G5.9b) for the simulators and
# summary_stats(); they run on CRAN and CI.

# multiplies by 1 + 1e-10 relative noise, so zeros stay zeros
sims_jitter <- function(x) x * (1 + stats::runif(length(x), -1e-10, 1e-10))

test_that("summary_stats is insensitive to trivial noise in its numeric inputs", {
  # G5.9a: abundances, trait, environment and coordinates jittered by 1e-10 relative; with the
  # same seed for the null model of ses_mpd and ses_mntd, every statistic changes by less
  # than 1e-6 relative (the jitter propagated through ratios and correlations is near 1e-10)
  set.seed(40)
  tree <- ape::rcoal(8)
  comm <- matrix(stats::rpois(20 * 8, 3), 20, 8, dimnames = list(NULL, tree$tip.label))
  series <- array(stats::rpois(20 * 8 * 4, 3), c(20, 8, 4),
                  dimnames = list(NULL, tree$tip.label, NULL))
  trait <- stats::setNames(stats::rnorm(8), tree$tip.label)
  env <- stats::runif(20)
  coords <- cbind(stats::runif(20), stats::runif(20))
  # jitter drawn before the seed is set, so the null model gets the same random numbers
  jittered <- list(comm = sims_jitter(comm), trait = sims_jitter(trait), env = sims_jitter(env),
                   coords = sims_jitter(coords))
  set.seed(1)
  base <- summary_stats(comm, tree, trait, env, coords, n_null = 19)
  set.seed(1)
  noisy <- summary_stats(jittered$comm, tree, jittered$trait, jittered$env, jittered$coords,
                         n_null = 19)
  expect_equal(noisy, base, tolerance = 1e-6)
  expect_equal(summary_stats(sims_jitter(series)), summary_stats(series), tolerance = 1e-6)
})

test_that("the null-model statistics of summary_stats vary with the seed within MC error", {
  # G5.9b: ses_mpd with 99 shuffles under ten seeds. Its Monte Carlo error comes from the
  # estimated null mean and standard deviation of each site, about 1 / sqrt(99) = 0.1 per
  # site before averaging over 15 sites; the ten values must have a standard deviation below
  # 0.2 (development: 0.06) and a mean within 0.1 of the value with 999 shuffles
  # (development: 0.03)
  set.seed(41)
  tree <- ape::rcoal(12)
  comm <- matrix(stats::rbinom(15 * 12, 1, 0.5), 15, 12, dimnames = list(NULL, tree$tip.label))
  values <- vapply(1:10, function(seed) {
    set.seed(seed)
    summary_stats(comm, tree, n_null = 99)[["ses_mpd"]]
  }, numeric(1))
  set.seed(99)
  reference <- summary_stats(comm, tree, n_null = 999)[["ses_mpd"]]
  expect_lt(stats::sd(values), 0.2)
  expect_lt(abs(mean(values) - reference), 0.1)
})

test_that("simulators are insensitive to trivial noise in traits and growth rates", {
  # G5.9a: expected abundances of simulate_static change by about 1e-10 relative; Poisson
  # counts with the same seed are identical, as a jitter of 1e-10 moves a quantile boundary
  # past the uniform draw with negligible probability
  set.seed(42)
  traits <- species_traits(6, optima_range = c(10, 90), rmax_range = c(2, 20))
  env <- stats::runif(30, 0, 100)
  noisy <- traits
  noisy[] <- lapply(traits, sims_jitter)
  noisy_env <- sims_jitter(env)
  expect_equal(simulate_static(noisy, noisy_env, sampling = "expected")$abundance,
               simulate_static(traits, env, sampling = "expected")$abundance,
               tolerance = 1e-8)
  set.seed(1)
  a <- simulate_static(traits, env)$abundance
  set.seed(1)
  b <- simulate_static(noisy, noisy_env)$abundance
  expect_identical(b, a)
  # mc_dynamics: same seed, growth rates and competition jittered
  r <- matrix(stats::runif(10 * 3, 2, 5), 10, 3)
  alpha <- matrix(0.01, 3, 3) + diag(0.02, 3)
  noisy_r <- sims_jitter(r)
  noisy_alpha <- sims_jitter(alpha)
  set.seed(2)
  a <- mc_dynamics(r, alpha, dispersal = 0.05, n_steps = 30)$abundance
  set.seed(2)
  b <- mc_dynamics(noisy_r, noisy_alpha, dispersal = 0.05, n_steps = 30)$abundance
  expect_identical(b, a)
})
