test_that("env_pulse changes only the disturbed sites during the pulse", {
  e <- env_pulse(c(0.6, 0.7, 0.8), n_steps = 10, sites = c(1, 3), start = 4, duration = 3,
                 value = 0.1)
  expect_equal(dim(e), c(3L, 10L))
  expect_equal(e[2, ], rep(0.7, 10))
  expect_equal(e[1, ], c(rep(0.6, 3), rep(0.1, 3), rep(0.6, 4)))
  expect_equal(e[3, 4:6], rep(0.1, 3))
  # permanent change, one value per site, logical sites
  p <- env_pulse(c(0.6, 0.7, 0.8), 5, sites = c(TRUE, FALSE, TRUE), start = 2,
                 value = c(0.2, 0.3))
  expect_equal(p[, 5], c(0.2, 0.7, 0.3))
  expect_equal(p[, 1], c(0.6, 0.7, 0.8))
})

test_that("env_pulse keeps the course of an environment that varies in time", {
  set.seed(1)
  et <- env_over_time(c(10, 20), 20, rho = 0.5, sd = 1)
  e <- env_pulse(et, sites = 2, start = 5, duration = 5, value = 0)
  expect_equal(e[1, ], et[1, ])
  expect_equal(e[2, -(5:9)], et[2, -(5:9)])
  expect_true(all(e[2, 5:9] == 0))
})

test_that("env_pulse rejects wrong sites and steps", {
  expect_error(env_pulse(c(1, 2), 5, sites = 3, start = 1, value = 0), "indices")
  expect_error(env_pulse(c(1, 2), 5, sites = TRUE, start = 1, value = 0), "one value per site")
  expect_error(env_pulse(matrix(1, 2, 5), 6, sites = 1, start = 1, value = 0), "n_steps")
  expect_error(env_pulse(c(1, 2), 5, sites = 1, start = 6, value = 0), "after the last")
  expect_error(env_pulse(c(1, 2), 5, sites = 1, start = 2.5, value = 0), "whole number")
  expect_error(env_pulse(c(1, 2), 5, sites = 1, start = 1, duration = 0.5, value = 0),
               "whole number")
  expect_error(env_pulse(c(1, 2), sites = 1, start = 1, value = 0), "n_steps")
  expect_error(env_pulse(c(1, 2), 5, sites = 1, start = 1, value = NA_real_), "value")
  expect_error(env_pulse(c(1, NA), 5, sites = 1, start = 1, value = 0), "NA")
})

test_that("a permanent change reaches the last step", {
  e <- env_pulse(c(1, 2), 5, sites = 1, start = 2, value = 0)
  expect_equal(e[1, ], c(1, 0, 0, 0, 0))
})

test_that("a pulse below every optimum empties disturbed sites, which are recolonised", {
  set.seed(3)
  ls <- landscape_grid(side = 4, disp_scale = 10)
  traits <- data.frame(optimum = seq(0.55, 0.95, length.out = 6), sigma = 0.25 * sqrt(2),
                       r_max = 5, row.names = paste0("sp", 1:6))
  env <- stats::runif(16, 0.5, 1)
  et <- env_pulse(env, 160, sites = 1:8, start = 60, duration = 50, value = 0.1)
  sim <- simulate_metacommunity(ls, traits, competition_matrix(6, "random", inter_max = 0.5),
                                dispersal = 0.05, n_steps = 160, env_time = et,
                                seed_steps = numeric(0), seed_mean = 1,
                                immigration = "multinomial", record_steps = c(59, 109, 160))
  richness <- apply(sim$history > 0, c(1, 3), sum)
  expect_true(all(richness[1:8, "step109"] < richness[1:8, "step59"]))
  expect_gt(mean(richness[1:8, "step160"]), mean(richness[1:8, "step109"]))
})

test_that("env_pulse without disturbed sites returns the environment unchanged", {
  expect_equal(env_pulse(c(1, 2), 5, sites = integer(0), start = 1, value = 0),
               matrix(c(1, 2), 2, 5))
  expect_equal(env_pulse(c(1, 2), 5, sites = c(FALSE, FALSE), start = 1, value = 0),
               matrix(c(1, 2), 2, 5))
})
