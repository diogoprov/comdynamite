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
               "single trait")
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
