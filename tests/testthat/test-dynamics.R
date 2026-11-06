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

test_that("print.mc_sim reports the simulation and returns it invisibly", {
  set.seed(1)
  tr <- ape::rcoal(4)
  sim <- simulate_metacommunity(landscape_implicit(3), species_traits(tr),
                                competition_matrix(tr, "none"), dispersal = 0.01,
                                n_steps = 3, record_steps = 1:3)
  expect_output(out <- print(sim), "3 sites x 4 species")
  expect_output(print(sim), "recorded steps: 3")
  expect_identical(out, sim)
})

test_that("extirpation acts on each species and leaves other runs unchanged", {
  r <- matrix(5, 30, 2)
  a <- diag(0.05, 2)
  set.seed(1)
  out <- mc_dynamics(r, a, dispersal = 0, n_steps = 5, seed_steps = integer(0), seed_mean = 20,
                     extirpation = c(0, 1))
  expect_true(all(out$abundance[, 1] > 0))
  expect_true(all(out$abundance[, 2] == 0))

  set.seed(2)
  none <- mc_dynamics(r, a, dispersal = 0.1, n_steps = 20)
  set.seed(2)
  zero <- mc_dynamics(r, a, dispersal = 0.1, n_steps = 20, extirpation = 0)
  expect_identical(none, zero)

  expect_error(mc_dynamics(r, a, dispersal = 0, n_steps = 5, extirpation = 2))
})

test_that("extirpation reaches the metacommunity and is kept with the result", {
  set.seed(3)
  tree <- ape::rcoal(5)
  land <- landscape_implicit(40)
  traits <- species_traits(tree)
  alpha <- competition_matrix(tree, "none")
  with_loss <- simulate_metacommunity(land, traits, alpha, dispersal = 0, n_steps = 150,
                                      extirpation = 0.2)
  without <- simulate_metacommunity(land, traits, alpha, dispersal = 0, n_steps = 150)
  expect_equal(with_loss$extirpation, 0.2)
  expect_lt(sum(with_loss$abundance > 0), sum(without$abundance > 0))
})
