test_that("a species alone converges to N* = (r - 1) / alpha", {
  set.seed(2)
  out <- mc_dynamics(matrix(5, 50, 1), matrix(0.05), dispersal = 0, n_steps = 300)
  expect_equal(mean(out$abundance), 80, tolerance = 0.05)
})

test_that("uniform dispersal moves individuals without creating or losing them", {
  r <- matrix(1, 20, 3)
  a <- matrix(0, 3, 3)
  set.seed(9)
  n1 <- mc_dynamics(r, a, dispersal = 0, n_steps = 1, seed_steps = integer(0), seed_mean = 3)
  set.seed(9)
  n2 <- mc_dynamics(r, a, dispersal = 0.5, n_steps = 1, seed_steps = integer(0), seed_mean = 3)
  expect_equal(colSums(n1$abundance), colSums(n2$abundance))
  expect_false(isTRUE(all.equal(n1$abundance, n2$abundance)))
})

test_that("recorded history matches the final state and keeps extinct species at zero", {
  set.seed(3)
  ls <- landscape_implicit(15)
  tr <- ape::rcoal(12)
  sim <- simulate_metacommunity(ls, species_traits(tr), competition_matrix(tr, "random", 0.5),
                                n_steps = 150, record_steps = c(50, 150))
  expect_s3_class(sim, "mc_sim")
  expect_equal(dim(sim$history), c(15, 12, 2))
  expect_equal(sim$history[, , 2], sim$abundance, ignore_attr = TRUE)
  expect_identical(colnames(sim$abundance), tr$tip.label)
  expect_true(all(sim$abundance >= 0))
})

test_that("time-varying environment and growth-rate functions run", {
  set.seed(4)
  ls <- landscape_implicit(10)
  tr <- species_traits(8)
  et <- env_over_time(ls$env, 100, rho = 0.9, sd = 5)
  expect_equal(dim(et), c(10, 100))
  sim <- simulate_metacommunity(ls, tr, competition_matrix(8), n_steps = 100, env_time = et)
  expect_equal(dim(sim$abundance), c(10, 8))
  expect_error(simulate_metacommunity(ls, tr, competition_matrix(8), n_steps = 100,
                                      env_time = et[, 1:50]))
})

test_that("spatially explicit simulation runs on a grid", {
  set.seed(5)
  ls <- landscape_grid(side = 5)
  tr <- ape::rcoal(10)
  sim <- simulate_metacommunity(ls, species_traits(tr), competition_matrix(tr),
                                dispersal = 0.3, n_steps = 80)
  expect_equal(dim(sim$abundance), c(25, 10))
  expect_output(print(sim), "25 sites x 10 species")
})
