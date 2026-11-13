start <- function(n, m) matrix(n, m, 1, dimnames = list(NULL, "sp1"))
evo_traits <- function(optimum = 0, sigma = 1e6, r_max = 5) {
  data.frame(optimum = optimum, sigma = sigma, r_max = r_max, row.names = "sp1")
}

test_that("without genetic variance, evolution leaves abundances as they were", {
  set.seed(1)
  ls <- landscape_implicit(15)
  tr <- ape::rcoal(6)
  traits <- species_traits(tr)
  a <- competition_matrix(tr, "phylo")
  set.seed(2)
  s0 <- simulate_metacommunity(ls, traits, a, n_steps = 150, extirpation = 0.01)
  for (nu in c(Inf, 1)) {
    set.seed(2)
    s1 <- simulate_metacommunity(ls, traits, a, n_steps = 150, extirpation = 0.01,
                                 evolution = evolution_control(0, ne_ratio = nu))
    expect_identical(s1$abundance, s0$abundance)
    present <- s1$abundance > 0
    expect_equal(s1$evolution$genetic[present],
                 matrix(traits$optimum, 15, 6, byrow = TRUE)[present])
    expect_true(all(is.na(s1$evolution$genetic[!present])))
  }
})

test_that("drift alone accumulates G times the sum of 1 / Ne (Lande 1976)", {
  set.seed(3)
  m <- 1000
  ev <- evolution_control(genetic_var = 1, ne_ratio = 0.5, founder = FALSE,
                          env = rep(0, m), traits = evo_traits())
  steps <- 40
  out <- mc_dynamics(matrix(5, m, 1), matrix(0.05), dispersal = 0, n_steps = steps,
                     seed_steps = numeric(0), initial = start(80L, m),
                     record_steps = seq_len(steps), evolution = ev)
  n <- out$history[, 1, ]
  ok <- apply(n > 0, 1, all)
  expected <- mean(rowSums(1 / (0.5 * n[ok, ])))
  expect_equal(var(out$evolution$genetic[ok, 1]), expected, tolerance = 0.15)
})

test_that("two populations reach the divergence of migration-selection balance", {
  set.seed(4)
  env <- c(0, 10)
  a <- 4 / (10^2 + 4)
  p <- 0.05
  ev <- evolution_control(genetic_var = 4, ne_ratio = Inf, env = env,
                          traits = evo_traits(optimum = 5, sigma = 10))
  out <- mc_dynamics(matrix(5, 2, 1), matrix(0.001), dispersal = p, n_steps = 1500,
                     seed_steps = numeric(0), initial = start(4000L, 2), evolution = ev)
  d <- unname(diff(out$evolution$genetic[, 1]))
  expected <- (1 - 2 * p) * a * diff(env) / (1 - (1 - 2 * p) * (1 - a))
  expect_equal(d, expected, tolerance = 0.05)
})

test_that("full plasticity around the optimum removes selection on the genetic mean", {
  set.seed(5)
  ev <- evolution_control(genetic_var = 4, ne_ratio = Inf, plasticity = 1, env = c(0, 10),
                          traits = evo_traits(optimum = 5, sigma = 3))
  out <- mc_dynamics(matrix(5, 2, 1), matrix(0.01), dispersal = 0.01, n_steps = 200,
                     seed_steps = numeric(0), initial = start(100L, 2), evolution = ev)
  expect_equal(out$evolution$genetic[, 1], c(site1 = 5, site2 = 5))
  expect_equal(unname(out$evolution$plastic[, 1]), c(-5, 5))
})

test_that("founder sampling makes new populations differ from their source", {
  m <- 200
  run <- function(founder) {
    set.seed(6)
    ev <- evolution_control(genetic_var = 1, ne_ratio = Inf, founder = founder,
                            env = rep(0, m), traits = evo_traits())
    mc_dynamics(matrix(5, m, 1), matrix(0.05), dispersal = 0.001, n_steps = 50,
                seed_steps = 1, seed_mean = 0.5, extirpation = 0.05, evolution = ev)
  }
  without <- run(FALSE)$evolution$genetic[, 1]
  with <- run(TRUE)$evolution$genetic[, 1]
  expect_equal(unique(stats::na.omit(without)), 0)
  expect_gt(stats::sd(with, na.rm = TRUE), 0.05)
})

test_that("steps per generation slow the response", {
  run <- function(per) {
    set.seed(7)
    ev <- evolution_control(genetic_var = 4, ne_ratio = Inf, steps_per_generation = per,
                            founder = FALSE,
                            env = 10, traits = evo_traits(optimum = 0, sigma = 10))
    mc_dynamics(matrix(5, 1, 1), matrix(0.01), dispersal = 0, n_steps = 10,
                seed_steps = numeric(0), initial = start(100L, 1),
                evolution = ev)$evolution$genetic[1, 1]
  }
  a <- 4 / 104
  expect_equal(run(1), 10 * (1 - (1 - a)^10))
  expect_equal(run(2), 10 * (1 - (1 - a / 2)^10))
})

test_that("evolution is recorded and checked", {
  set.seed(8)
  ls <- landscape_grid(4)
  tr <- ape::rcoal(4)
  traits <- species_traits(tr)
  a <- competition_matrix(tr, "phylo")
  sim <- simulate_metacommunity(ls, traits, a, n_steps = 60, record_steps = c(30, 60),
                                evolution = evolution_control(1))
  expect_equal(dim(sim$evolution$history), c(16, 4, 2))
  expect_equal(sim$evolution$history[, , 2], sim$evolution$genetic)
  expect_s3_class(sim$evolution$settings, "mc_evolution")
  expect_output(print(evolution_control(1)), "mc_evolution")
  expect_error(evolution_control(), "genetic_var")
  expect_error(evolution_control(1, steps_per_generation = 0.5))
  expect_error(simulate_metacommunity(ls, traits, a, n_steps = 5, evolution = list(1)),
               "evolution_control")
  # several axes: one trait per axis, with G = diag(1)
  two <- simulate_metacommunity(ls, list(traits, traits), a, n_steps = 5,
                                env = cbind(ls$env, ls$env), evolution = evolution_control(1))
  expect_equal(dim(two$evolution$genetic), c(16, 4, 2))
  expect_error(mc_dynamics(matrix(5, 2, 1), matrix(0.05), 0, n_steps = 5,
                           evolution = evolution_control(1)), "env")
})

test_that("founders from divergent sources carry the variance among them", {
  em <- matrix(c(1, 1, 0))
  z <- matrix(c(0, 10, NA))
  mean_in <- immigrant_mean(em, z, NULL, z)
  between <- immigrant_mean(em, z^2, NULL, z^2) - mean_in^2
  expect_equal(mean_in[3], 5)
  expect_equal(between[3], 25)
})

test_that("evolution follows env_time and checks its settings", {
  set.seed(9)
  ls <- landscape_implicit(5)
  tr <- ape::rcoal(3)
  traits <- species_traits(tr)
  a <- competition_matrix(tr, "phylo")
  env_time <- matrix(rep(seq(0, 50, length.out = 40), each = 5), 5, 40)
  sim <- simulate_metacommunity(ls, traits, a, n_steps = 40, env_time = env_time,
                                evolution = evolution_control(1, plasticity = 1))
  present <- sim$abundance > 0
  expect_equal(sim$evolution$plastic[present],
               (50 - mean(env_time[, 1])) * rep(1, sum(present)))
  ev <- function(...) evolution_control(1, ...)
  r <- matrix(5, 2, 1)
  expect_error(mc_dynamics(r, matrix(0.05), 0, n_steps = 2, evolution = list(1)),
               "evolution_control")
  expect_error(mc_dynamics(r, matrix(0.05), 0, n_steps = 2,
                           evolution = ev(env = 1:3, traits = evo_traits())), "one value per site")
  expect_error(mc_dynamics(r, matrix(0.05), 0, n_steps = 2,
                           evolution = ev(env = 1:2, traits = rbind(evo_traits(), evo_traits()))),
               "one row per species")
  expect_error(mc_dynamics(r, matrix(0.05), 0, n_steps = 2,
                           evolution = evolution_control(c(1, 2), env = 1:2,
                                                         traits = evo_traits())),
               "one value or one per species")
})

test_that("simulate_metacommunity() checks alpha and dispersal", {
  ls <- landscape_implicit(3)
  tr <- ape::rcoal(2)
  traits <- species_traits(tr)
  expect_error(simulate_metacommunity(ls, traits, diag(3), n_steps = 2), "one row and one column")
  a <- diag(2)
  dimnames(a) <- list(c("x", "y"), c("x", "y"))
  expect_error(simulate_metacommunity(ls, traits, a, n_steps = 2), "names of species")
  expect_error(mc_dynamics(matrix(5, 1, 1), matrix(0.05), dispersal = 0.1, n_steps = 2),
               "at least two sites")
})
