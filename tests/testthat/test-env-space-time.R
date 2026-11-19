test_that("env_space_time has the shape, burn-in and ranks of Thompson et al. (2020)", {
  set.seed(1)
  xy <- cbind(runif(20, 1, 100), runif(20, 1, 100))
  e <- env_space_time(xy, n_steps = 300, burn_in = 50, n_features = 500)
  expect_equal(dim(e), c(20, 300))
  expect_true(all(e > 0 & e <= 1))
  expect_equal(max(e), 1)
  expect_true(all(e[, 1:51] == e[, 1]))
  expect_false(all(e[, 52] == e[, 51]))
  # ranks: uniform over all sites and steps
  expect_lt(abs(mean(e) - 0.5), 0.05)
  raw <- env_space_time(xy, n_steps = 300, burn_in = 50, ranks = FALSE, n_features = 500)
  expect_equal(range(raw), c(0, 1))
  expect_gt(diff(range(raw[, 51])), 0.6)
})

test_that("env_space_time accepts a landscape and checks its input", {
  set.seed(2)
  land <- landscape_grid(side = 4)
  e <- env_space_time(land, n_steps = 20, burn_in = 5, min_range = 0, n_features = 200)
  expect_equal(dim(e), c(16, 20))
  expect_error(env_space_time(landscape_implicit(5), 10), "no coordinates")
  expect_error(env_space_time(matrix(1:3), 10), "two columns")
  expect_error(env_space_time(cbind(1:3, 1:3), 10, burn_in = 10))
  # one site has no range at the first step, so `min_range` can never be met
  expect_error(env_space_time(cbind(1, 1), 3, burn_in = 1, n_features = 5), "lower it")
})

test_that("the spectral field has the exponential covariance", {
  set.seed(3)
  # site-steps on the scale of the covariance: distances 0.6 in space and 0.1 in time
  p <- rbind(c(0, 0, 0), c(0.6, 0, 0), c(0, 0, 0.1), c(0.3, 0.4, 1.2))
  draws <- replicate(4000, comdynamite:::exp_field(p, 200))
  v <- tcrossprod(draws) / ncol(draws)
  expect_equal(diag(v), rep(1, 4), tolerance = 0.08)
  expect_equal(v[1, 2:4], exp(-c(0.6, 0.1, 1.3)), tolerance = 0.1)
})
