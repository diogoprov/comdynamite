test_that("env_seasonal() adds a cycle with amplitudes per site and yearly anomalies", {
  e <- env_seasonal(c(10, 20), n_steps = 40, period = 20, amplitude = c(1, 4))
  expect_equal(dim(e), c(2, 40))
  expect_equal(e[, 5], c(10, 20) + c(1, 4) * sin(2 * pi * 5 / 20))
  expect_equal(rowMeans(e), c(10, 20))
  expect_equal(e[, 1:20], e[, 21:40])
  set.seed(1)
  y <- env_seasonal(c(10, 20), n_steps = 60, period = 20, amplitude = 0, year_sd = 2)
  # one anomaly per year, shared by the sites
  expect_equal(length(unique(round(y[1, ] - 10, 10))), 3)
  expect_equal(y[1, ] - 10, y[2, ] - 20)
  base <- matrix(1:6, 2, 3)
  expect_equal(env_seasonal(base, period = 4, amplitude = 1, phase = 1)[, 1], c(1, 2))
  expect_error(env_seasonal(1:3, period = 5, amplitude = 1), "n_steps")
  expect_error(env_seasonal(base, n_steps = 5, period = 5, amplitude = 1), "differs")
  expect_error(env_seasonal(c(1, NA), n_steps = 5, period = 5, amplitude = 1), "numeric")
})

test_that("river_floods() gives washout probabilities by season, year and river size", {
  net <- data.frame(downstream = c(NA, 1, 2, 2), length = c(0, 1, 1, 3))
  land <- landscape_river(net)
  f <- river_floods(land, n_steps = 30, period = 10, flood_steps = 2, start = 3,
                    probability = 0.6)
  expect_equal(dim(f), c(4, 30))
  expect_equal(which(colSums(f) > 0), c(3, 4, 13, 14, 23, 24))
  expect_equal(attr(f, "flood_years"), 1:3)
  size <- land$network$upstream_length
  expect_equal(f[, 3], 0.6 * size / max(size))
  flat <- river_floods(land, n_steps = 10, period = 10, size_effect = 0, probability = 0.2)
  expect_equal(flat[, 1], rep(0.2, 4))
  head <- river_floods(land, n_steps = 10, period = 10, size_effect = -1, size = "magnitude",
                       probability = 0.5)
  expect_true(all(head[, 1] <= 1))
  expect_gt(head[3, 1], head[1, 1])
  set.seed(2)
  some <- river_floods(land, n_steps = 1000, period = 10, year_prob = 0.3)
  expect_equal(length(attr(some, "flood_years")) / 100, 0.3, tolerance = 0.25)
  expect_error(river_floods(landscape_grid(3), 10, 5), "landscape_river")
  expect_error(river_floods(land, 10, 5, flood_steps = 6), "at most")
  expect_error(river_floods(land, 10, 5, size = c(1, 2)), "size")
})

test_that("simulate_metacommunity() takes disturbance as a matrix of sites x steps", {
  set.seed(3)
  land <- landscape_river(river_network(12))
  traits <- data.frame(optimum = 50, sigma = 40, r_max = 3, row.names = "sp1")
  alpha <- matrix(0.05, dimnames = list("sp1", "sp1"))
  wipe <- matrix(0, 12, 30)
  wipe[1:6, 30] <- 1
  sim <- simulate_metacommunity(land, traits, alpha, n_steps = 30, disturbance = wipe)
  expect_true(all(sim$abundance[1:6, ] == 0))
  expect_true(any(sim$abundance[7:12, ] > 0))
  # a matrix of zeros gives the same run as no disturbance
  set.seed(4)
  a <- simulate_metacommunity(land, traits, alpha, n_steps = 20)
  set.seed(4)
  b <- simulate_metacommunity(land, traits, alpha, n_steps = 20,
                              disturbance = matrix(0, 12, 20))
  expect_identical(a$abundance, b$abundance)
  expect_error(simulate_metacommunity(land, traits, alpha, n_steps = 20,
                                      disturbance = matrix(0, 12, 5)), "matrix")
  expect_error(simulate_metacommunity(land, traits, alpha, n_steps = 20,
                                      disturbance = matrix(2, 12, 20)), "probabilit")
})
