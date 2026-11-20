test_that("competition_matrix(type = 'trait') is the overlap of Gaussian niches", {
  tr <- data.frame(optimum = c(10, 10, 40), sigma = c(5, 5, 10), r_max = 5,
                   row.names = c("a", "b", "c"))
  a <- competition_matrix(type = "trait", traits = tr, inter_max = 1, scale = 1)
  expect_equal(dimnames(a), list(c("a", "b", "c"), c("a", "b", "c")))
  expect_equal(unname(diag(a)), c(1, 1, 1))
  expect_equal(a["a", "b"], 1)
  expect_equal(a["a", "c"], exp(-30^2 / (2 * (25 + 100))))
  expect_equal(a, t(a))
  a2 <- competition_matrix(type = "trait", traits = tr, inter_max = 0.5)
  expect_equal(a2["a", "c"], 0.05 * 0.5 * exp(-30^2 / 250))
  expect_equal(a2["a", "a"], 0.05)
  # several axes: product of the overlaps
  tr2 <- tr
  tr2$optimum <- c(0, 20, 0)
  aa <- competition_matrix(type = "trait", traits = list(tr, tr2), inter_max = 1, scale = 1)
  expect_equal(aa["a", "b"], exp(-400 / (2 * 50)))
  expect_equal(aa["a", "c"], a["a", "c"])
  expect_error(competition_matrix(type = "trait"), "needs `traits`")
  expect_error(competition_matrix(3, "random", traits = tr), "only with")
  expect_error(competition_matrix(type = "trait", traits = data.frame(x = 1)), "optimum")
})

test_that("the same kernel as niche_overlap in the speciation engine", {
  tr <- data.frame(optimum = c(20, 35), sigma = c(8, 12), r_max = 5,
                   row.names = c("a", "b"))
  a <- competition_matrix(type = "trait", traits = tr, inter_max = 1, scale = 1)
  w <- 8^2 + 12^2
  expect_equal(a["a", "b"], exp(-(15)^2 / (2 * w)))
})

test_that("niche_correlation gives a multivariate Gaussian niche", {
  tr <- list(t = data.frame(optimum = c(40, 60), sigma = c(10, 20), r_max = c(5, 4),
                            row.names = c("a", "b")),
             o = data.frame(optimum = c(40, 60), sigma = c(10, 5), r_max = 1,
                            row.names = c("a", "b")))
  env <- cbind(t = c(30, 50, 30), o = c(30, 50, 50))
  base <- growth_rates(env, tr)
  expect_equal(growth_rates(env, niche_correlation(tr, diag(2))), base)
  rho <- 0.6
  r <- growth_rates(env, niche_correlation(tr, matrix(c(1, rho, rho, 1), 2)))
  d <- c(-1, -1)
  expect_equal(unname(r[1, "a"]), 5 * exp(-drop(d %*% solve(matrix(c(1, rho, rho, 1), 2)) %*% d) / 2))
  # same-direction departures are favoured, opposite ones penalised
  expect_gt(r[1, "a"], base[1, "a"])
  expect_lt(r[3, "a"], base[3, "a"])
  # one matrix per species
  rr <- growth_rates(env, niche_correlation(tr, list(b = diag(2), a = matrix(c(1, rho, rho, 1), 2))))
  expect_equal(rr[, "a"], r[, "a"])
  expect_equal(rr[, "b"], base[, "b"])
  expect_error(niche_correlation(tr, matrix(c(1, 2, 2, 1), 2)), "positive definite")
  expect_error(niche_correlation(tr, matrix(c(1, 0.2, 0.3, 1), 2)), "symmetric")
  expect_error(niche_correlation(tr[1], diag(2)), "at least two")
  expect_error(niche_correlation(tr, list(a = diag(2))), "one matrix per species")
})

test_that("landscape_grid gives several correlated axes, and simulations use them", {
  set.seed(4)
  one <- landscape_grid(side = 12)
  set.seed(4)
  three <- landscape_grid(side = 12, n_axes = 3, axis_cor = 0.7)
  expect_equal(three$env, one$env)
  expect_equal(dim(three$env_axes), c(144, 3))
  expect_equal(three$env_axes[, 1], one$env)
  expect_true(all(three$env_axes >= 0 & three$env_axes <= 100))
  expect_null(one$env_axes)
  set.seed(5)
  r <- replicate(20, {
    l <- landscape_grid(side = 12, n_axes = 2, axis_cor = 0.8)
    cor(l$env_axes)[1, 2]
  })
  expect_gt(mean(r), 0.5)
  expect_error(landscape_grid(side = 4, n_axes = 3, axis_cor = -0.9), "between")
  expect_error(landscape_grid(side = 4, n_axes = 2, axis_cor = matrix(1, 3, 3)), "positive")

  set.seed(6)
  l <- landscape_grid(side = 5, n_axes = 2, axis_cor = 0.5)
  tr <- list(species_traits(3), species_traits(3, rmax_range = c(1, 1)))
  a <- competition_matrix(3)
  set.seed(7)
  s1 <- simulate_metacommunity(l, tr, a, dispersal = 0.1, n_steps = 30)
  set.seed(7)
  s2 <- simulate_metacommunity(l, tr, a, dispersal = 0.1, n_steps = 30, env = l$env_axes)
  expect_equal(s1$abundance, s2$abundance)
  set.seed(7)
  s3 <- simulate_metacommunity(l, niche_correlation(tr, matrix(c(1, 0.9, 0.9, 1), 2)), a,
                               dispersal = 0.1, n_steps = 30)
  expect_false(isTRUE(all.equal(s1$abundance, s3$abundance)))
  expect_error(simulate_metacommunity(l, niche_correlation(tr, diag(2)), a, dispersal = 0.1,
                                      n_steps = 5, evolution = evolution_control(10)),
               "correlated axes")
})

test_that("the remaining checks of the new options", {
  # unnamed species and mismatched sets of traits
  tr <- data.frame(optimum = c(10, 40), sigma = c(5, 10))
  expect_equal(rownames(competition_matrix(type = "trait", traits = tr)), c("sp1", "sp2"))
  a <- data.frame(optimum = 1:2, sigma = 1, row.names = c("a", "b"))
  b <- data.frame(optimum = 1:2, sigma = 1, row.names = c("a", "c"))
  expect_error(competition_matrix(type = "trait", traits = list(a, b)), "different species")
  # correlation matrices
  tr2 <- list(species_traits(2), species_traits(2, rmax_range = c(1, 1)))
  expect_error(niche_correlation(tr2, diag(3)), "2 x 2")
  set.seed(8)
  l <- landscape_grid(side = 4, n_axes = 2, axis_cor = matrix(c(1, -0.5, -0.5, 1), 2))
  expect_equal(dim(l$env_axes), c(16, 2))
  # speciation: axes from the landscape, and no correlated niche yet
  set.seed(9)
  l2 <- landscape_grid(side = 3, n_axes = 2)
  one <- list(species_traits(1), species_traits(1, rmax_range = c(1, 1)))
  ev <- evolution_control(genetic_var = diag(2) * 5, env_var = diag(2))
  sim <- simulate_evolving_metacommunity(l2, one, dispersal = 0.1, n_steps = 3,
                                         speciation = speciation_control(threshold = 50),
                                         evolution = ev)
  expect_s3_class(sim, "mc_evolving")
  expect_error(simulate_evolving_metacommunity(l2, niche_correlation(one, diag(2)),
                                               dispersal = 0.1, n_steps = 3,
                                               speciation = speciation_control(threshold = 50),
                                               evolution = ev),
               "correlated axes")
})
