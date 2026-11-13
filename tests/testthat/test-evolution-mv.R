mv_setup <- function(m = 25) {
  ls <- landscape_grid(sqrt(m), disp_scale = 0.5)
  list(ls = ls, env = cbind(svl = ls$env, head = 0),
       tr = list(svl = data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1"),
                 head = data.frame(optimum = 0, sigma = Inf, r_max = 1, row.names = "sp1")),
       alpha = matrix(0.05, dimnames = list("sp1", "sp1")))
}

test_that("G = 0 gives the abundances without evolution", {
  s <- mv_setup()
  set.seed(3)
  a <- simulate_metacommunity(s$ls, s$tr, s$alpha, env = s$env, dispersal = 0.01, n_steps = 50)
  set.seed(3)
  b <- simulate_metacommunity(s$ls, s$tr, s$alpha, env = s$env, dispersal = 0.01, n_steps = 50,
                              evolution = evolution_control(matrix(0, 2, 2)))
  expect_identical(a$abundance, b$abundance)
  expect_equal(dim(b$evolution$genetic), c(25, 1, 2))
})

test_that("drift makes the covariance among isolated populations G times sum(1 / Ne)", {
  set.seed(2)
  m <- 400
  ls <- landscape_grid(20, disp_scale = 0.5)
  tr <- data.frame(optimum = 0, sigma = Inf, r_max = 2, row.names = "sp1")
  G <- matrix(c(1, 0.6, 0.6, 0.5), 2)
  sim <- simulate_metacommunity(ls, list(a = tr, b = tr),
                                matrix(0.02, dimnames = list("sp1", "sp1")),
                                env = cbind(a = ls$env, b = 0), dispersal = 0, n_steps = 60,
                                evolution = evolution_control(G, founder = FALSE),
                                seed_steps = 0, record_steps = 1:60,
                                initial = matrix(50, m, 1, dimnames = list(NULL, "sp1")))
  expected <- G * mean(rowSums(1 / sim$history[, 1, -1]))
  observed <- stats::cov(sim$evolution$genetic[, 1, ])
  expect_equal(unname(observed), expected, tolerance = 0.15)
})

test_that("selection moves the unselected trait by its genetic covariance", {
  s <- mv_setup()
  G <- matrix(c(10, 6, 6, 8), 2)
  sim <- simulate_metacommunity(s$ls, s$tr, s$alpha, env = s$env, dispersal = 0, n_steps = 30,
                                evolution = evolution_control(G, ne_ratio = Inf,
                                                              founder = FALSE),
                                seed_steps = 0, initial = matrix(20, 25, 1,
                                                                 dimnames = list(NULL, "sp1")))
  z <- sim$evolution$genetic[, 1, ]
  # Delta z = G[, 1] (omega + P11)^-1 d: the two traits move in the ratio G21 / G11
  expect_equal(unname(z[, "head"]), unname(z[, "svl"] - 50) * 6 / 10, tolerance = 1e-8)
  expect_gt(stats::cor(z[, "svl"], s$ls$env, use = "complete.obs"), 0.9)
  # one step of the deterministic response
  one <- simulate_metacommunity(s$ls, s$tr, s$alpha, env = s$env, dispersal = 0, n_steps = 1,
                                evolution = evolution_control(G, ne_ratio = Inf,
                                                              founder = FALSE),
                                seed_steps = 0, initial = matrix(20, 25, 1,
                                                                 dimnames = list(NULL, "sp1")))
  d <- s$env[, "svl"] - 50
  expect_equal(unname(one$evolution$genetic[, 1, "svl"]), 50 + 10 / (20^2 + 10) * d)
})

test_that("multivariate settings are checked", {
  s <- mv_setup()
  run <- function(ev, traits = s$tr, env = s$env) {
    simulate_metacommunity(s$ls, traits, s$alpha, env = env, dispersal = 0, n_steps = 2,
                           evolution = ev)
  }
  expect_error(run(evolution_control(matrix(c(1, 2, 2, 1), 2))), "positive semidefinite")
  expect_error(run(evolution_control(matrix(c(1, 0.5, 0, 1), 2))), "symmetric")
  expect_error(run(evolution_control(diag(3))), "symmetric 2 x 2")
  expect_error(run(evolution_control(list(diag(2), diag(2)))), "one matrix per species")
  expect_error(run(evolution_control(diag(2), plasticity = matrix(0, 2, 2))), "plasticity")
  expect_error(simulate_metacommunity(s$ls, s$tr$svl, s$alpha, dispersal = 0, n_steps = 2,
                                      evolution = evolution_control(diag(2))),
               "several axes")
  expect_error(evolution_control(diag(2), env_var = -1))
  ev <- evolution_control(list(sp1 = diag(2)), env_var = diag(2), plasticity = c(0.5, 0))
  expect_output(print(ev), "1 matrices")
  expect_output(print(evolution_control(diag(2))), "2 x 2 matrix")
  sim <- run(ev)
  expect_equal(unname(sim$evolution$plasticity[1, ]), c(0.5, 0))
  sp <- speciation_control("point", rate = 0.01)
  expect_error(simulate_evolving_metacommunity(s$ls, s$tr$svl, sp,
                                               evolution = evolution_control(diag(2))),
               "several axes")
})

test_that("sample_individuals() draws several traits from N(mean, G + E)", {
  set.seed(5)
  s <- mv_setup()
  G <- matrix(c(10, 6, 6, 8), 2)
  E <- diag(2, 2)
  sim <- simulate_metacommunity(s$ls, s$tr, s$alpha / 20, env = s$env, dispersal = 0.01,
                                n_steps = 40, record_steps = c(20, 40),
                                evolution = evolution_control(G, env_var = E, ne_ratio = 0.5))
  ind <- sample_individuals(sim, n = 1e4, source = "evolution")
  expect_true(all(c("env_svl", "mean_svl", "svl", "env_head", "mean_head", "head") %in%
                    names(ind)))
  expect_identical(attr(ind, "traits"), c("svl", "head"))
  dev <- as.matrix(ind[, c("svl", "head")] - ind[, c("mean_svl", "mean_head")])
  expect_equal(unname(stats::cov(dev)), G + E, tolerance = 0.08)
  ind2 <- sample_individuals(sim, n = 1e4, source = "evolution", intraspecific_sd = c(1, 0))
  dev2 <- as.matrix(ind2[, c("svl", "head")] - ind2[, c("mean_svl", "mean_head")])
  expect_equal(stats::var(dev2[, 1]), 13, tolerance = 0.08)
  dates <- sample_individuals(sim, n = 3, source = "evolution", step = c(20, 40),
                              session_sd = 1)
  expect_identical(sort(unique(dates$step)), c(20, 40))
  at20 <- sample_individuals(sim, n = 3, source = "evolution", step = 20)
  expect_equal(unique(at20$mean_svl[at20$site == "site1"]),
               unname(sim$evolution$history["site1", "sp1", "svl", "step20"]))
  expect_error(sample_individuals(sim, source = "evolution", step = 10), "not recorded")
  expect_error(sample_individuals(sim, source = "evolution", env = s$env[, 1]), "leave `env`")
  expect_error(sample_individuals(sim, source = "evolution", population_sd = 1), "stay 0")
  expect_error(sample_individuals(sim, source = "evolution", intraspecific_sd = c(1, 2, 3)),
               "one per trait")
  frac <- sample_individuals(sim, fraction = 0.5, source = "evolution", cover = 1)
  expect_true(nrow(frac) > 0)
})

test_that("genetic_matrices() builds the three families", {
  set.seed(1)
  p <- genetic_matrices(3, n_traits = 3, correlation = 0.5, heritability = 0.4)
  expect_named(p$genetic_var, c("sp1", "sp2", "sp3"))
  expect_equal(p$genetic_var$sp1 / p$env_var$sp1, matrix(0.4 / 0.6, 3, 3))
  sh <- genetic_matrices(c("a", "b"), n_traits = 3, family = "shared", correlation = 0.3,
                         divergence = 0.5)
  g <- sh$genetic_var$a
  e <- sh$env_var$a
  expect_equal(g %*% e, e %*% g, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(g / e, matrix(1, 3, 3))))
  P <- matrix(c(2, 1, 1, 3), 2, dimnames = list(c("x", "y"), c("x", "y")))
  ind <- genetic_matrices(2, family = "independent", P = P, heritability = 0.5,
                          divergence = 0.3)
  expect_equal(sum(diag(ind$genetic_var$sp1)), 2.5)
  expect_equal(sum(diag(ind$env_var$sp2)), 2.5)
  expect_identical(dimnames(ind$genetic_var$sp1), dimnames(P))
  expect_true(all(eigen(ind$genetic_var$sp1)$values > -1e-10))
  expect_false(isTRUE(all.equal(ind$genetic_var$sp1 %*% ind$env_var$sp1,
                                ind$env_var$sp1 %*% ind$genetic_var$sp1)))
  expect_error(genetic_matrices(2, heritability = 2))
  expect_error(genetic_matrices(2, P = matrix(c(1, 1, 1, 1), 2)), "positive definite")
  s <- mv_setup()
  m <- genetic_matrices("sp1", correlation = 0.5)
  sim <- simulate_metacommunity(s$ls, s$tr, s$alpha, env = s$env, dispersal = 0, n_steps = 3,
                                evolution = evolution_control(m$genetic_var,
                                                              env_var = m$env_var))
  expect_equal(unname(sim$evolution$genetic_var$sp1), unname(m$genetic_var$sp1))
})

test_that("several species, migrants and plasticity follow the same engine", {
  set.seed(11)
  ls <- landscape_grid(4, disp_scale = 0.5)
  sp <- c("sp1", "sp2", "sp3")
  tr <- list(a = data.frame(optimum = c(40, 50, 60), sigma = c(15, 20, 25), r_max = 4,
                            row.names = sp),
             b = data.frame(optimum = 0, sigma = c(Inf, 5, 10), r_max = c(1, 1.5, 1),
                            row.names = sp))
  env <- cbind(a = ls$env, b = seq(-2, 2, length.out = 16))
  alpha <- matrix(0.03, 3, 3, dimnames = list(sp, sp))
  run <- function(ev, seed = 4, ...) {
    set.seed(seed)
    simulate_metacommunity(ls, tr, alpha, env = env, dispersal = 0.05, extirpation = 0.02,
                           n_steps = 40, evolution = ev, record_steps = c(20, 40), ...)
  }
  base <- run(NULL)
  zero <- run(evolution_control(matrix(0, 2, 2)))
  expect_identical(base$abundance, zero$abundance)
  expect_identical(base$history, zero$history)
  expect_true(all(is.na(zero$evolution$genetic[, , 1]) == (zero$abundance == 0)))
  # G as a named list out of order is matched by species
  Gs <- list(sp3 = diag(c(3, 1)), sp1 = diag(2), sp2 = matrix(c(2, 1, 1, 2), 2))
  sim <- run(evolution_control(Gs, env_var = diag(2), plasticity = c(0.3, 0)))
  expect_equal(unname(sim$evolution$genetic_var$sp3), diag(c(3, 1)))
  plastic <- sim$evolution$plastic[, , "a"]
  ok <- !is.na(plastic[, 1])
  expect_equal(unname(plastic[ok, 1]), (env[ok, "a"] - mean(env[, "a"])) * 0.3)
  expect_true(all(sim$evolution$plastic[, , "b"] == 0, na.rm = TRUE))
  # implicit landscape: immigrants from all other sites
  imp <- landscape_implicit(6)
  t1 <- list(a = tr$a[1, , drop = FALSE], b = tr$b[1, , drop = FALSE])
  set.seed(2)
  s_imp <- simulate_metacommunity(imp, t1, alpha[1, 1, drop = FALSE],
                                  env = cbind(a = imp$env * 10 + 40, b = 0),
                                  dispersal = 0.2, n_steps = 20,
                                  evolution = evolution_control(diag(c(4, 1))))
  expect_equal(dim(s_imp$evolution$genetic), c(6, 1, 2))
  expect_error(run(evolution_control(list(sp1 = diag(2), sp2 = diag(2), spx = diag(2)))),
               "lacks matrices for")
  expect_error(run(evolution_control(diag(2), plasticity = c(0.1, 0.2, 0.3))), "plasticity")
  bm <- matrix(c(0.1, 0.2, 0.3, 0, 0, 0), 3, dimnames = list(c("sp3", "sp1", "sp2"), NULL))
  sim_b <- run(evolution_control(diag(2), plasticity = bm), seed = 1)
  expect_equal(unname(sim_b$evolution$plasticity[, 1]), c(0.2, 0.3, 0.1))
})

test_that("the engine refuses inconsistent multivariate settings", {
  s <- mv_setup()
  expect_error(mc_dynamics(matrix(5, 2, 1), matrix(0.05), 0, n_steps = 5,
                           evolution = c(evolution_control(diag(2)), list(env = 1:2,
                                         traits = data.frame(optimum = 1, sigma = 1,
                                                             r_max = 1)))),
               "evolution_control")
  ev <- evolution_control(diag(2))
  ev$env <- 1:2
  ev$traits <- data.frame(optimum = 1, sigma = 1, r_max = 1)
  expect_error(mc_dynamics(matrix(5, 2, 1), matrix(0.05), 0, n_steps = 5, evolution = ev),
               "several axes")
  bad <- s$tr
  names(bad) <- c("site", "head")
  env <- s$env
  colnames(env) <- c("site", "head")
  expect_error(simulate_metacommunity(s$ls, bad, s$alpha, env = env, dispersal = 0,
                                      n_steps = 2, evolution = ev), "distinct names")
  zero <- s$tr
  zero$svl$sigma <- 0
  expect_error(simulate_metacommunity(s$ls, zero, s$alpha, env = s$env, dispersal = 0,
                                      n_steps = 2, evolution = evolution_control(0)),
               "sigma")
})

test_that("the multivariate state checks traits and env, and axes without selection have no load", {
  s <- mv_setup()
  ev <- evolution_control(diag(2))
  ev$traits <- s$tr
  ev$env <- s$env
  expect_error(evolution_state_mv(ev, 25, 2, NULL), "one row per species")
  ev$env <- s$env[1:5, ]
  expect_error(evolution_state_mv(ev, 25, 1, NULL), "one row per site")
  free <- s$tr
  free$svl$sigma <- Inf
  set.seed(1)
  sim <- simulate_metacommunity(s$ls, free, s$alpha, env = s$env, dispersal = 0.01,
                                n_steps = 10, evolution = evolution_control(diag(2)))
  m <- maladaptation(sim)
  expect_true(all(m$distance == 0 & m$load == 0))
})

test_that("env_time changes every axis of a niche, with and without evolution", {
  set.seed(6)
  s <- mv_setup()
  tr <- s$tr
  tr$head$sigma <- 5
  n_steps <- 30
  warm <- list(svl = sapply(1:n_steps, function(t) s$env[, "svl"] + t),
               head = sapply(1:n_steps, function(t) s$env[, "head"] - 0.1 * t))
  run <- function(env_time, ev = NULL, seed = 1) {
    set.seed(seed)
    simulate_metacommunity(s$ls, tr, s$alpha, env = s$env, env_time = env_time,
                           dispersal = 0.02, n_steps = n_steps, record_steps = c(10, 30),
                           evolution = ev)
  }
  # constant env_time equals no env_time
  same <- lapply(s$env |> as.data.frame(), function(x) matrix(x, length(x), n_steps))
  expect_identical(run(same)$abundance, run(NULL)$abundance)
  # as an array, in a different order of the list
  arr <- array(c(warm$svl, warm$head), c(25, n_steps, 2))
  expect_identical(run(warm[c("head", "svl")])$abundance, run(arr)$abundance)
  ev <- evolution_control(diag(c(4, 1)), env_var = diag(2))
  sim <- run(warm, ev)
  eh <- sim$evolution$environment_history
  expect_equal(unname(eh[, "svl", "step10"]), s$env[, "svl"] + 10)
  expect_equal(unname(eh[, "head", "step30"]), s$env[, "head"] - 3)
  expect_equal(unname(sim$evolution$environment[, "svl"]), s$env[, "svl"] + 30)
  m <- maladaptation(sim, step = 10)
  expect_equal(m$env_svl, unname(s$env[match(m$site, rownames(sim$abundance)), "svl"] + 10))
  # G = 0 with env_time gives the abundances without evolution
  expect_identical(run(warm, evolution_control(matrix(0, 2, 2)))$abundance,
                   run(warm)$abundance)
  expect_error(run(list(warm$svl)), "one matrix")
  expect_error(run(list(svl = warm$svl, head = warm$head[, 1:5])), "columns")
})

test_that("env_time as a matrix works with a single column of env and one evolving trait", {
  s <- mv_setup()
  n_steps <- 20
  warm <- sapply(1:n_steps, function(t) s$ls$env + t)
  run <- function(...) {
    set.seed(8)
    simulate_metacommunity(s$ls, s$tr$svl, s$alpha, dispersal = 0.01, n_steps = n_steps,
                           evolution = evolution_control(4), record_steps = c(10, 20), ...)
  }
  a <- run(env = data.frame(svl = s$ls$env), env_time = warm)
  b <- run(env_time = warm)
  expect_equal(a$evolution$genetic, b$evolution$genetic)
  expect_equal(unname(a$evolution$environment_history[, "step10"]), warm[, 10])
})

test_that("immigrant_selection is refused with several traits and in evolving simulations", {
  s <- mv_setup()
  ev <- evolution_control(genetic_var = diag(2), immigrant_selection = TRUE)
  expect_error(simulate_metacommunity(s$ls, s$tr, s$alpha, env = s$env, n_steps = 5,
                                      evolution = ev), "several traits")
  expect_error(simulate_evolving_metacommunity(s$ls, s$tr$svl,
                                               speciation_control("geographic", Inf),
                                               evolution = evolution_control(
                                                 10, immigrant_selection = TRUE),
                                               n_steps = 5), "not yet available")
  expect_error(simulate_metacommunity(s$ls, s$tr, s$alpha, env = s$env, n_steps = 5,
                                      evolution = evolution_control(diag(2), ne_ratio = 0.5,
                                                                    mutation_var = 0.1)),
               "several traits")
})
