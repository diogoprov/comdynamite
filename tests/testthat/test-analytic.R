# Results that can be derived by hand from the model. Equilibria come from
# setting the expectation of r N / (1 + alpha N) equal to N.

test_that("a species alone converges to (r - 1) / alpha and dies out if r < 1", {
  set.seed(11)
  out <- mc_dynamics(matrix(3, 60, 1), matrix(0.02), dispersal = 0, n_steps = 300)
  expect_equal(mean(out$abundance), (3 - 1) / 0.02, tolerance = 0.02)
  lost <- mc_dynamics(matrix(0.8, 20, 1), matrix(0.02), dispersal = 0, n_steps = 300)
  expect_equal(sum(lost$abundance), 0)
})

test_that("coexisting species reach the solution of alpha N = r - 1", {
  set.seed(12)
  alpha <- matrix(c(1, 0.3, 0.2,
                    0.1, 1, 0.4,
                    0.3, 0.2, 1), 3, byrow = TRUE) * 0.02
  r <- c(4, 5, 6)
  expected <- solve(alpha, r - 1)
  stopifnot(all(expected > 0))
  out <- mc_dynamics(matrix(r, 50, 3, byrow = TRUE), alpha, dispersal = 0, n_steps = 400)
  expect_equal(unname(colMeans(out$abundance)), expected, tolerance = 0.02)
})

test_that("a species that cannot invade is excluded", {
  set.seed(13)
  alpha <- matrix(c(1, 0.5,
                    1.5, 1), 2, byrow = TRUE) * 0.02
  # species 2 cannot invade: at N1* its growth is r / (1 + alpha21 N1*) < 1
  out <- mc_dynamics(matrix(5, 30, 2), alpha, dispersal = 0, n_steps = 600)
  expect_equal(sum(out$abundance[, 2]), 0)
  expect_equal(mean(out$abundance[, 1]), 4 / 0.02, tolerance = 0.02)
})

test_that("with equal competition, total abundance matches a single species", {
  set.seed(14)
  out <- mc_dynamics(matrix(5, 40, 6), matrix(0.05, 6, 6), dispersal = 0, n_steps = 300)
  expect_equal(mean(rowSums(out$abundance)), 4 / 0.05, tolerance = 0.03)
})

test_that("a source sustains a sink at the mass-effect equilibrium", {
  # r = 0 in the sink: N_source = ((1 - d) r - 1) / alpha, N_sink = d N_source / (1 - d)
  set.seed(15)
  d <- 0.2
  source_n <- ((1 - d) * 5 - 1) / 0.05
  sink_n <- d * source_n / (1 - d)
  steps <- 200:1200
  implicit <- mc_dynamics(matrix(c(5, 0), 2, 1), matrix(0.05), dispersal = d,
                          n_steps = 1200, record_steps = steps)
  expect_equal(rowMeans(implicit$history[, 1, ]), c(source_n, sink_n),
               tolerance = 0.03, ignore_attr = TRUE)

  pairs <- kronecker(diag(30), matrix(c(0, 1, 1, 0), 2))
  explicit <- mc_dynamics(matrix(rep(c(5, 0), 30), 60, 1), matrix(0.05), dispersal = d,
                          n_steps = 300, disp_matrix = pairs)
  expect_equal(mean(explicit$abundance[c(TRUE, FALSE), ]), source_n, tolerance = 0.03)
  expect_equal(mean(explicit$abundance[c(FALSE, TRUE), ]), sink_n, tolerance = 0.03)
})

test_that("one step of dispersal has the expected flow through the kernel", {
  set.seed(16)
  d <- 0.3
  r <- c(1, 2, 3)
  kernel <- matrix(c(0, 0.5, 0.5,
                     0.9, 0, 0.1,
                     0.2, 0.8, 0), 3, byrow = TRUE)
  start <- 5000
  expected <- (1 - d) * r * start + d * drop(crossprod(kernel, r * start))
  sims <- replicate(200, mc_dynamics(matrix(r, 3, 1), matrix(0), dispersal = d, n_steps = 1,
                                     seed_steps = integer(0), seed_mean = start,
                                     disp_matrix = kernel)$abundance[, 1])
  expect_equal(unname(rowMeans(sims)), expected, tolerance = 0.005)
})

test_that("mean static abundance equals the Gaussian response", {
  set.seed(17)
  traits <- data.frame(optimum = 50, sigma = 10, r_max = 20, row.names = "sp1")
  env <- rep(c(50, 60, 70), each = 3000)
  sim <- simulate_static(traits, env)
  means <- tapply(sim$abundance[, 1], env, mean)
  expect_equal(unname(c(means)), 20 * exp(-c(0, 10, 20)^2 / (2 * 10^2)), tolerance = 0.02)
})

test_that("isolated populations survive extirpation geometrically", {
  # r = 5 and alpha = 0.05 keep populations near 80, so only extirpation removes them
  set.seed(20)
  q <- 0.01
  steps <- 50
  sites <- 5000
  out <- mc_dynamics(matrix(5, sites, 1), matrix(0.05), dispersal = 0, n_steps = steps,
                     seed_steps = integer(0), seed_mean = 20, extirpation = q)
  expected <- (1 - q)^steps
  expect_lt(abs(mean(out$abundance > 0) - expected), 4 * sqrt(expected * (1 - expected) / sites))
})
