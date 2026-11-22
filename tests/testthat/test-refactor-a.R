# Paths of the simulation engines that no other test runs, and the helpers that the engines
# were split into.

test_that("one trait evolves along one axis of `env` that changes through time", {
  set.seed(1)
  ls <- landscape_implicit(6)
  traits <- data.frame(optimum = c(40, 60), sigma = 20, r_max = 5, row.names = c("sp1", "sp2"))
  a <- matrix(0.05, 2, 2, dimnames = list(rownames(traits), rownames(traits)))
  env <- matrix(ls$env, ncol = 1, dimnames = list(NULL, "temp"))
  env_time <- list(temp = matrix(rep(seq(20, 80, length.out = 30), each = 6), 6, 30))
  set.seed(4)
  sim <- simulate_metacommunity(ls, traits, a, n_steps = 30, env = env, env_time = env_time,
                                evolution = evolution_control(genetic_var = 4, ne_ratio = 0.5),
                                record_steps = c(10, 30))
  expect_s3_class(sim, "mc_sim")
  # the environment of the evolving trait is that of the axis at each step
  expect_identical(sim$evolution$environment, env_time$temp[, 30])
  expect_identical(unname(sim$evolution$environment_history[, "step10"]), env_time$temp[, 10])
  expect_identical(dim(sim$evolution$genetic), c(6L, 2L))
  expect_named(sim$env, "temp")
  # and the fixed one without `env_time`
  set.seed(4)
  fixed <- simulate_metacommunity(ls, traits, a, n_steps = 20, env = env,
                                  evolution = evolution_control(genetic_var = 4, ne_ratio = 0.5))
  expect_identical(fixed$evolution$environment, as.numeric(ls$env))
  expect_identical(fixed$evolution$reference, mean(ls$env))
})

test_that("records of a step follow the buddings of that step with several traits", {
  set.seed(2)
  grid <- landscape_grid(4, disp_scale = 0.5)
  env <- cbind(svl = grid$env, head = seq(-3, 3, length.out = 16))
  anc <- list(svl = data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1"),
              head = data.frame(optimum = 0, sigma = 3, r_max = 1, row.names = "sp1"))
  set.seed(5)
  sim <- simulate_evolving_metacommunity(
    grid, anc, speciation_control("fission", rate = 0.2),
    evolution = evolution_control(matrix(c(10, 2, 2, 1), 2), env_var = diag(2)),
    env = env, dispersal = 0.02, n_steps = 40, record_steps = 1:40
  )
  k <- ncol(sim$abundance)
  expect_gt(nrow(sim$events), 0)
  expect_identical(dim(sim$history), c(16L, k, 40L))
  expect_identical(dim(sim$evolution$history), c(16L, k, 2L, 40L))
  # the last record is the state at the end, after the buddings of the last step
  expect_identical(unname(sim$history[, , 40]), unname(sim$abundance))
  expect_identical(sim$evolution$history[, , , 40], sim$evolution$genetic)
})

test_that("a run that starts without individuals has an empty series", {
  set.seed(2)
  grid <- landscape_grid(4, disp_scale = 0.5)
  one <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
  set.seed(6)
  sim <- simulate_evolving_metacommunity(
    grid, one, speciation_control("point", rate = 0.001), n_steps = 5,
    initial = matrix(0, 16, 1, dimnames = list(NULL, "sp1"))
  )
  expect_s3_class(sim, "mc_evolving")
  expect_identical(nrow(sim$series), 0L)
  expect_true(all(c("generation", "lineage", "occupancy", "abundance") %in% names(sim$series)))
  expect_identical(sim$richness, rep(0L, 5))
  expect_identical(sum(sim$abundance), 0L)
  expect_identical(nrow(sim$events), 0L)
  expect_null(sim$tree)
})

test_that("helpers of the checks of the engines return what the engines use", {
  expect_identical(detection_by_species(0.5, c("a", "b")), c(a = 0.5, b = 0.5))
  expect_identical(detection_by_species(c(b = 0.2, a = 0.1), c("a", "b")), c(a = 0.1, b = 0.2))
  expect_error(detection_by_species(c(0.1, 0.2), c("a", "b")), "named by species")
  expect_identical(evolution_each(2, 3), c(2, 2, 2))
  expect_error(evolution_each(c(1, 2), 3), "one value or one per")
  expect_true(is_covariance(diag(2)))
  expect_false(is_covariance(1))
  expect_true(is_variances(c(0, 1)))
  expect_false(is_variances(-1))
  expect_identical(mc_lottery_microsites(5, 3), c(5L, 5L, 5L))
  expect_error(mc_lottery_microsites(0.5, 3), "needs `microsites`")
  expect_identical(mc_sensitivity(1, 2), c(1, 1))
  expect_error(mc_sensitivity(2, 2), "disturbance_sensitivity")
  steps <- matrix(1, 4, 2, dimnames = list(NULL, c("fish", "frog")))
  expect_identical(kernel_group_columns(c("frog", "fish"), steps), c(2L, 1L))
  expect_identical(kernel_group_columns(c(1, 2), steps), c(1L, 2L))
  expect_error(kernel_group_columns(c(1, 3), steps), "a column of `kernel_step`")
  expect_identical(dispersal_time_matrix(matrix(0.1, 3, 1), 2, 2), matrix(0.1, 2, 2))
  expect_error(dispersal_time_matrix(matrix(2, 3, 1), 2, 2), "emigration probabilities")
})

test_that("arrays of the state grow by species", {
  x <- matrix(1:4, 2)
  expect_identical(grow_species_array(x, NA_integer_, 1L), cbind(x, NA_integer_))
  arr <- array(1:8, c(2, 2, 2))
  grown <- grow_species_array(arr, 0L, 2)
  expect_identical(dim(grown), c(2L, 4L, 2L))
  expect_identical(grown[, 1:2, ], arr)
  expect_identical(grown[, 3:4, ], array(0L, c(2, 2, 2)))
  expect_null(grow_species_history(NULL, 0L, "sp2"))
  hist <- array(1L, c(2, 1, 3), dimnames = list(c("s1", "s2"), "sp1", paste0("step", 1:3)))
  grown <- grow_species_history(hist, 0L, c("sp2", "sp3"))
  expect_identical(dimnames(grown)[[2]], c("sp1", "sp2", "sp3"))
  expect_identical(grown[, "sp1", ], hist[, "sp1", ])
  expect_identical(unname(grown[, "sp3", ]), matrix(0L, 2, 3))
})
