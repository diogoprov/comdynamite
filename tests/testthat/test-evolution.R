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

test_that("immigrant_selection filters maladapted immigrants source by source", {
  xy <- expand.grid(x = 0:4, y = 0:3)
  env <- rep(c(40, 60), 10)
  ls <- landscape_points(xy, env, disp_scale = 1)
  tr <- data.frame(optimum = c(40, 60), sigma = 5, r_max = 5, row.names = c("sp1", "sp2"))
  init <- cbind(sp1 = (env == 40) * 100L, sp2 = (env == 60) * 100L)
  run <- function(sel) {
    set.seed(4)
    simulate_metacommunity(ls, tr, matrix(0.04, 2, 2), dispersal = 0.2, n_steps = 30,
                           seed_steps = integer(0), initial = init,
                           evolution = evolution_control(genetic_var = 1, ne_ratio = Inf,
                                                         immigrant_selection = sel))
  }
  plain <- run(FALSE)
  filtered <- run(TRUE)
  # maladapted immigrants are mostly filtered: fewer individuals of each species in the
  # patch type of the other, and their populations closer to the local environment
  wrong <- function(s) sum(s$abundance[env == 60, "sp1"]) + sum(s$abundance[env == 40, "sp2"])
  expect_lt(wrong(filtered), wrong(plain))
  # the helper weights each source by its fitness in the destination
  disp <- matrix(c(0, 1, 1, 0), 2)
  sel <- immigrants_selected(matrix(c(10, 0), 2, 1), matrix(c(40, 0), 2, 1),
                             matrix(c(40, 0), 2, 1), c(40, 60), 25, disp)
  expect_equal(sel$lambda[2, 1], 10 * exp(-400 / 50))
  expect_equal(sel$mean[2, 1], 40)
  # a species without emigrants sends no one and keeps its own mean
  none <- immigrants_selected(matrix(0, 2, 1), matrix(c(40, 45), 2, 1), matrix(c(40, 45), 2, 1),
                              c(40, 60), 25, disp)
  expect_equal(none$lambda, matrix(0, 2, 1))
  expect_equal(none$mean, matrix(c(40, 45), 2, 1))
  expect_error(evolution_control(1, immigrant_selection = NA))
  expect_error(simulate_metacommunity(landscape_implicit(5), tr, matrix(0.04, 2, 2),
                                      n_steps = 5,
                                      evolution = evolution_control(1, immigrant_selection = TRUE)),
               "spatially explicit")
})

test_that("mutation_var gives each population a genetic variance with mutation-drift dynamics", {
  ls <- landscape_points(expand.grid(x = 1:3, y = 1:2), rep(50, 6), disp_scale = 1)
  tr <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
  a <- matrix(0.01, dimnames = list("sp1", "sp1"))
  init <- cbind(sp1 = rep(400L, 6))
  run <- function(vm, g = 10, n_steps = 300, nu = 0.5) {
    set.seed(5)
    simulate_metacommunity(ls, tr, a, dispersal = 0, n_steps = n_steps,
                           seed_steps = integer(0), initial = init,
                           evolution = evolution_control(genetic_var = g, ne_ratio = nu,
                                                         mutation_var = vm))
  }
  # equilibrium 2 Ne Vm: Ne = 0.5 x N (N about 395), so Vm = 0.025 keeps G near 10
  eq <- run(0.025, n_steps = 1000)
  expect_equal(mean(eq$evolution$genetic_variance), 0.025 * mean(eq$abundance), tolerance = 0.1)
  # no mutation: the variance decays by about 1 / (2 Ne) per generation
  dec <- run(0, n_steps = 100)
  expected <- 10 * (1 - 1 / mean(dec$abundance))^100
  expect_equal(mean(dec$evolution$genetic_variance), expected, tolerance = 0.03)
  # constant G has no genetic_variance
  expect_null(run(NULL)$evolution$genetic_variance)
  expect_error(evolution_control(1, mutation_var = -1), "mutation_var")
  expect_error(evolution_control(1, ne_ratio = Inf, mutation_var = 0.1), "finite")
  # founders: a population founded by few immigrants starts with less variance
  xy <- expand.grid(x = 0:3, y = 0:1)
  ls2 <- landscape_points(xy, rep(50, 8), disp_scale = 0.5)
  init2 <- cbind(sp1 = c(400L, rep(0L, 7)))
  set.seed(6)
  sim <- simulate_metacommunity(ls2, tr, a, dispersal = 0.01, n_steps = 40,
                                seed_steps = integer(0), initial = init2,
                                evolution = evolution_control(genetic_var = 10, ne_ratio = 0.5,
                                                              mutation_var = 0.025,
                                                              immigrant_selection = TRUE))
  v <- sim$evolution$genetic_variance
  expect_true(any(v[-1] < v[1], na.rm = TRUE))
  # growth with one variance per population matches the vector version
  p <- matrix(c(40, 55), 2, 1)
  expect_equal(evolution_growth_cells(p, c(50, 50), 20, 5, matrix(3, 2, 1)),
               evolution_growth(p, c(50, 50), 20, 5, 3))
})

test_that("mutation_var works with seeding and with gene flow without immigrant selection", {
  ls <- landscape_points(expand.grid(x = 1:4, y = 1:3), rep(c(40, 60), 6), disp_scale = 1)
  tr <- data.frame(optimum = c(45, 55), sigma = 20, r_max = 5, row.names = c("sp1", "sp2"))
  set.seed(7)
  sim <- simulate_metacommunity(ls, tr, matrix(0.01, 2, 2), dispersal = 0.05, n_steps = 150,
                                evolution = evolution_control(genetic_var = 10, ne_ratio = 0.5,
                                                              mutation_var = 0.025))
  v <- sim$evolution$genetic_variance
  present <- sim$abundance > 0
  expect_true(any(present))
  expect_true(all(is.finite(v[present]) & v[present] > 0))
  expect_true(all(is.na(v[!present])))
})
