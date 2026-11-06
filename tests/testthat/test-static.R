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

test_that("sampling = 'expected' returns the expected abundances", {
  traits <- data.frame(optimum = 20, sigma = 10, r_max = 0.8, row.names = "a")
  sim <- simulate_static(traits, c(20, 30), sampling = "expected")
  expect_equal(unname(sim$abundance[, 1]), 0.8 * exp(-c(0, 10)^2 / 200))
  expect_error(simulate_static(traits, 20, incidence = TRUE, sampling = "expected"))
})

test_that("logistic presences follow their probabilities", {
  set.seed(5)
  tr <- ape::compute.brlen(ape::stree(8, "balanced"), method = "Grafen")
  g <- logistic_gradients(1)[c("x", "y")]
  sim <- simulate_logistic(tr, g, slope_sd = 0.3)
  expect_s3_class(sim, "mc_sim")
  expect_equal(dim(sim$abundance), c(31L, 8L))
  expect_true(all(sim$abundance %in% 0:1))
  expect_true(all(sim$probability >= 0 & sim$probability <= 1))
  expect_equal(dim(sim$traits$slope), c(8L, 2L))
  # no slopes, no intercepts, no reduction: every probability is 0.25
  flat <- simulate_logistic(tr, g, intercept_sd = 0, slope_sd = 0, reduction = c(1, 1))
  expect_equal(unname(flat$probability), matrix(0.25, 31, 8))
  many <- replicate(300, simulate_logistic(tr, g, intercept_sd = 0, slope_sd = 0,
                                           reduction = c(1, 1))$abundance)
  expect_equal(mean(many), 0.25, tolerance = 0.03)
  normal <- replicate(100, simulate_logistic(tr, 1:10, intercept_mean = -30, slope_sd = 0,
                                             intercept_sd = 0, draw = "normal")$abundance)
  expect_equal(mean(normal), 0.5, tolerance = 0.03)
})

test_that("slopes carry phylogenetic signal and centring makes no difference", {
  set.seed(6)
  tr <- ape::compute.brlen(ape::stree(16, "balanced"), method = "Grafen")
  b <- replicate(2000, simulate_logistic(tr, 1:5)$traits$slope[, 1])
  r <- cor(t(b))
  sisters <- r[cbind(c(1, 3, 5), c(2, 4, 6))]
  expect_gt(mean(sisters), 0.5)
  expect_equal(mean(apply(b, 1, sd)), 1, tolerance = 0.1)
  set.seed(7); p1 <- simulate_logistic(tr, 1:5)$probability
  set.seed(7); p2 <- simulate_logistic(tr, 101:105)$probability
  expect_equal(p1, p2)
})

test_that("the seven scenarios have the expected gradients", {
  set.seed(8)
  crossings <- function(v) sum(diff(sign(v)) != 0)
  s <- lapply(1:7, logistic_gradients)
  expect_true(all(vapply(s, nrow, 1L) == 31))
  expect_equal(s[[1]]$x, -15:15)
  expect_equal(sort(s[[1]]$y), -15:15)
  expect_gt(cor(s[[2]]$x, s[[2]]$y), 0.95)
  expect_identical(s[[1]]$x_test, s[[1]]$x)
  expect_false(isTRUE(all.equal(s[[4]]$x_test, s[[4]]$x)))
  expect_lt(max(abs(s[[4]]$x_test)), 5)
  expect_equal(crossings(s[[5]]$y - mean(s[[5]]$y)), 1)
  expect_equal(crossings(s[[6]]$y - mean(s[[6]]$y)), 3)
  expect_equal(crossings(s[[7]]$x - mean(s[[7]]$x)), 1)
  expect_equal(range(s[[7]]$y), c(-15, 15))
  expect_error(logistic_gradients(8))
})

test_that("transect MEMs are orthonormal", {
  m <- transect_mem(20)
  expect_equal(crossprod(m), diag(ncol(m)), tolerance = 1e-8)
  expect_lt(ncol(m), 20)
})

test_that("binomial sampling draws presences with the expected probability", {
  set.seed(8)
  traits <- data.frame(optimum = 50, sigma = 10, r_max = 0.8, row.names = "sp1")
  env <- rep(c(50, 60), each = 5000)
  sim <- simulate_static(traits, env, sampling = "binomial")
  expect_setequal(sim$abundance, c(0, 1))
  expect_equal(unname(c(tapply(sim$abundance[, 1], env, mean))),
               0.8 * exp(-c(0, 10)^2 / (2 * 10^2)), tolerance = 0.03)
  traits$r_max <- 2
  expect_error(simulate_static(traits, env, sampling = "binomial"), "at most 1")
})
