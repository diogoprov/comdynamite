test_that("mean abundance matches the Gaussian response", {
  set.seed(1)
  traits <- data.frame(optimum = c(20, 80), sigma = c(10, 5), r_max = c(30, 10),
                       row.names = c("a", "b"))
  env <- c(20, 30, 80)
  sims <- replicate(2000, simulate_static(traits, env)$abundance)
  expected <- cbind(30 * exp(-c(0, 10, 60)^2 / 200), 10 * exp(-c(60, 50, 0)^2 / 50))
  expect_equal(unname(apply(sims, 1:2, mean)), expected, tolerance = 0.02)
})

test_that("output has sites by species, with names and class", {
  set.seed(2)
  traits <- species_traits(ape::rcoal(5), rmax_range = c(0, 30))
  sim <- simulate_static(traits, runif(7, 0, 100))
  expect_s3_class(sim, "mc_sim")
  expect_equal(dim(sim$abundance), c(7L, 5L))
  expect_equal(colnames(sim$abundance), rownames(traits))
  expect_equal(rownames(sim$abundance), paste0("site", 1:7))
  expect_output(print(sim), "7 sites x 5 species")
})

test_that("incidence gives zeros and ones", {
  set.seed(3)
  traits <- species_traits(5, rmax_range = c(5, 30))
  sim <- simulate_static(traits, seq(0, 100, 10), incidence = TRUE)
  expect_true(all(sim$abundance %in% 0:1))
  expect_gt(sum(sim$abundance), 0)
})

test_that("several gradients multiply the responses", {
  set.seed(4)
  a <- data.frame(optimum = c(20, 80), sigma = c(10, 5), r_max = c(30, 10),
                  row.names = c("a", "b"))
  b <- data.frame(optimum = c(0, 10), sigma = c(20, 20), r_max = c(1, 0.5),
                  row.names = c("a", "b"))
  env <- data.frame(e = c(20, 30, 80), s = c(0, 10, 20))
  sims <- replicate(2000, simulate_static(list(a, b), env)$abundance)
  expected <- growth_rates(env$e, a) * growth_rates(env$s, b)
  expect_equal(unname(apply(sims, 1:2, mean)), unname(expected), tolerance = 0.03)
  one <- simulate_static(list(a), list(env$e))
  expect_equal(dim(one$abundance), c(3L, 2L))
})

test_that("gradients and traits must match", {
  a <- data.frame(optimum = 1, sigma = 1, r_max = 1, row.names = "a")
  b <- data.frame(optimum = 1, sigma = 1, r_max = 1, row.names = "b")
  expect_error(simulate_static(list(a, a), list(1:3)))
  expect_error(simulate_static(list(a, a), list(1:3, 1:2)))
  expect_error(simulate_static(list(a, b), list(1:3, 1:3)))
  expect_error(simulate_static(a, "x"))
})
