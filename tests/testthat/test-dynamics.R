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

test_that("multinomial immigration conserves emigrants and follows the kernel", {
  set.seed(9)
  ls <- landscape_grid(side = 4)
  em <- matrix(rpois(16 * 3, 3), 16, 3)
  imm <- place_emigrants(em, ls$dispersal)
  expect_equal(colSums(imm), colSums(em))
  means <- Reduce(`+`, lapply(1:2000, function(i) place_emigrants(em, ls$dispersal))) / 2000
  expect_equal(means, crossprod(ls$dispersal, em), tolerance = 0.05)
  r <- matrix(3, 16, 3)
  set.seed(10)
  a <- mc_dynamics(r, diag(0.05, 3), 0.1, n_steps = 20, disp_matrix = ls$dispersal)
  set.seed(10)
  b <- mc_dynamics(r, diag(0.05, 3), 0.1, n_steps = 20, disp_matrix = ls$dispersal,
                   immigration = "poisson")
  expect_identical(a$abundance, b$abundance)
  sim <- simulate_metacommunity(ls, species_traits(3), diag(0.05, 3), dispersal = 0.2,
                                n_steps = 30, immigration = "multinomial")
  expect_s3_class(sim, "mc_sim")
  expect_error(mc_dynamics(r, diag(0.05, 3), 0.1, immigration = "other"), "should be one of")
})

test_that("a niche with several axes multiplies the responses to each axis", {
  set.seed(11)
  t1 <- species_traits(4)
  t2 <- species_traits(4, rmax_range = c(1, 1))
  env <- cbind(e1 = c(10, 50, 90), e2 = c(30, 30, 70))
  g <- growth_rates(env, list(t1, t2))
  expect_equal(g, growth_rates(env[, 1], t1) * growth_rates(env[, 2], t2), ignore_attr = TRUE)
  expect_equal(growth_rates(as.data.frame(env), list(t1, t2)), g)
  expect_equal(growth_rates(env[, 1, drop = FALSE], t1), growth_rates(env[, 1], t1),
               ignore_attr = TRUE)
  expect_error(growth_rates(env, list(t1)), "2 axes")
  expect_error(growth_rates(env, t1), "need a list of traits")
  t3 <- t2
  rownames(t3)[1] <- "other"
  expect_error(growth_rates(env, list(t1, t3)), "different species")

  ls <- landscape_implicit(5)
  alpha <- competition_matrix(4)
  set.seed(1)
  a <- simulate_metacommunity(ls, list(t1, t2), alpha, n_steps = 20,
                              env = cbind(ls$env, rep(50, 5)))
  set.seed(1)
  b <- mc_dynamics(growth_rates(cbind(ls$env, rep(50, 5)), list(t1, t2)), alpha,
                   dispersal = 0.01, n_steps = 20)
  expect_equal(unname(a$abundance), unname(b$abundance))
  expect_length(a$traits, 2)
  expect_equal(ncol(as_mc_data(a)$env), 2)
  # one axis given as env reproduces the landscape
  set.seed(2)
  c1 <- simulate_metacommunity(ls, t1, alpha, n_steps = 20)
  set.seed(2)
  c2 <- simulate_metacommunity(ls, t1, alpha, n_steps = 20, env = matrix(ls$env))
  expect_equal(c1$abundance, c2$abundance)
  expect_error(simulate_metacommunity(ls, list(t1, t2), alpha, n_steps = 5), "needs `env`")
  expect_error(simulate_metacommunity(ls, list(t1, t2), alpha, n_steps = 5,
                                      env = cbind(1:4, 1:4)), "one row per site")
  expect_error(simulate_metacommunity(ls, list(t1, t2), alpha, n_steps = 5,
                                      env = cbind(ls$env, ls$env),
                                      env_time = matrix(0, 5, 5)), "env_time")
})

test_that("a simulation with several niche axes keeps the tree and works downstream", {
  set.seed(12)
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  t1 <- species_traits(tree)
  t2 <- species_traits(tree, rmax_range = c(1, 1))
  ls <- landscape_implicit(6)
  alpha <- competition_matrix(tree, "none")
  sim <- simulate_metacommunity(ls, list(t1, t2), alpha, n_steps = 30,
                                env = cbind(ls$env, rev(ls$env)))
  expect_identical(sim$tree, attr(t1, "tree"))
  expect_named(sim$env, c("axis1", "axis2"))
  named <- simulate_metacommunity(ls, list(temp = t1, rain = t2), alpha, n_steps = 5,
                                  env = cbind(ls$env, rev(ls$env)))
  expect_named(named$env, c("temp", "rain"))
  d <- as_mc_data(sim)
  expect_identical(d$tree, sim$tree)
  expect_named(d$traits, c("optimum_axis1", "optimum_axis2"))
  expect_equal(d$traits$optimum_axis2, t2$optimum)
  # summary_stats() uses the first axis, as if given its trait and environment
  set.seed(3)
  s <- summary_stats(sim, n_null = 5)
  set.seed(3)
  s1 <- summary_stats(sim, trait = stats::setNames(t1$optimum, rownames(t1)),
                      env = sim$env$axis1, n_null = 5)
  expect_equal(s, s1)
  expect_error(sample_individuals(sim), "one environmental axis")
})

test_that("summary_stats() and sample_individuals() use the env that set growth", {
  set.seed(13)
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  t1 <- species_traits(tree)
  ls <- landscape_implicit(6)
  e <- rev(ls$env)
  sim <- simulate_metacommunity(ls, t1, competition_matrix(tree, "none"), n_steps = 30,
                                env = matrix(e))
  set.seed(4)
  s <- summary_stats(sim, n_null = 5)
  set.seed(4)
  expect_equal(s, summary_stats(sim, env = e, n_null = 5))
  ind <- sample_individuals(sim, n = 2)
  expect_equal(ind$env, e[match(ind$site, rownames(sim$abundance))])
})

test_that("disturbance empties whole sites and is validated", {
  r <- matrix(5, 30, 2)
  a <- matrix(0.05, 2, 2)
  # every site disturbed at the last step: nothing left
  all_gone <- mc_dynamics(r, a, dispersal = 0, n_steps = 20, seed_steps = integer(0),
                          initial = matrix(10L, 30, 2, dimnames = list(NULL, c("sp1", "sp2"))), disturbance = 1)
  expect_true(all(all_gone$abundance == 0))
  # without dispersal, a disturbed site stays empty
  set.seed(2)
  out <- mc_dynamics(r, a, dispersal = 0, n_steps = 50, seed_steps = integer(0),
                     initial = matrix(10L, 30, 2, dimnames = list(NULL, c("sp1", "sp2"))), disturbance = 0.02)
  empty <- rowSums(out$abundance) == 0
  expect_true(any(empty) && !all(empty))
  # no draws with disturbance = 0
  set.seed(3)
  a0 <- mc_dynamics(r, a, dispersal = 0.05, n_steps = 30)
  set.seed(3)
  b0 <- mc_dynamics(r, a, dispersal = 0.05, n_steps = 30, disturbance = 0)
  expect_identical(a0$abundance, b0$abundance)
  expect_error(mc_dynamics(r, a, dispersal = 0, disturbance = c(0.1, 0.2)), "disturbance")
  expect_error(mc_dynamics(r, a, dispersal = 0, disturbance = 2), "disturbance")
  ls <- landscape_implicit(10)
  tr <- data.frame(optimum = c(30, 70), sigma = 20, r_max = 5,
                   row.names = c("sp1", "sp2"))
  sim <- simulate_metacommunity(ls, tr, matrix(0.05, 2, 2), n_steps = 20, disturbance = 0.1)
  expect_equal(sim$disturbance, 0.1)
})

test_that("pool_from opens a site to a species at a given step", {
  ls <- landscape_matrix(matrix(1, 3, 3), env = c(40, 40, 60))
  tr <- data.frame(optimum = c(40, 40), sigma = 35, r_max = 5, row.names = c("sp1", "sp2"))
  init <- cbind(sp1 = c(500L, 0L, 0L), sp2 = c(0L, 500L, 0L))
  from <- cbind(sp1 = c(1, 1, 10), sp2 = c(1, 1, 30))
  set.seed(8)
  sim <- simulate_metacommunity(ls, tr, matrix(0.004, 2, 2), dispersal = 0.05, n_steps = 40,
                                seed_steps = integer(0), initial = init, pool_from = from,
                                record_steps = 1:40)
  expect_true(all(sim$history[3, "sp1", 1:9] == 0))
  expect_true(any(sim$history[3, "sp1", 10:15] > 0))
  expect_true(all(sim$history[3, "sp2", 1:29] == 0))
  expect_true(any(sim$history[3, "sp2", 30:40] > 0))
  expect_error(simulate_metacommunity(ls, tr, matrix(0.004, 2, 2), n_steps = 5,
                                      initial = init, pool_from = from * 0),
               "pool_from")
  bad <- from
  bad[1, 1] <- 5
  expect_error(simulate_metacommunity(ls, tr, matrix(0.004, 2, 2), n_steps = 5,
                                      initial = init, pool_from = bad), "opens after")
})

test_that("a site closed by pool_from is not a sink for emigrants", {
  tr <- data.frame(optimum = 40, sigma = 35, r_max = 5, row.names = "sp1")
  init <- cbind(sp1 = c(1000L, 0L, 0L))
  closed <- cbind(sp1 = c(1, 100, 100))
  for (ls in list(landscape_matrix(matrix(1, 3, 3), env = c(40, 40, 60)),
                  landscape_implicit(3, env_range = c(39, 41)))) {
    for (how in c("poisson", "multinomial")) {
      set.seed(3)
      sim <- simulate_metacommunity(ls, tr, matrix(0.004, 1, 1), dispersal = 0.5, n_steps = 20,
                                    seed_steps = integer(0), initial = init, pool_from = closed,
                                    immigration = how, record_steps = 1:20)
      # no emigration: the population stays near (5 - 1) / 0.004 = 1000
      expect_true(all(sim$history[1, 1, ] > 900))
      expect_true(all(sim$history[2:3, 1, ] == 0))
    }
  }
  # one site open: half of the emigrants go there, none are lost to the closed one
  set.seed(4)
  sim <- simulate_metacommunity(landscape_matrix(matrix(1, 3, 3), env = c(40, 40, 60)), tr,
                                matrix(0.004, 1, 1), dispersal = 0.5, n_steps = 30,
                                seed_steps = integer(0), initial = init,
                                pool_from = cbind(sp1 = c(1, 1, 100)))
  expect_true(all(sim$abundance[1:2, 1] > 800))
  # the same on a spatially implicit landscape
  set.seed(5)
  sim <- simulate_metacommunity(landscape_implicit(3, env_range = c(39, 41)), tr,
                                matrix(0.004, 1, 1), dispersal = 0.5, n_steps = 30,
                                seed_steps = integer(0), initial = init,
                                pool_from = cbind(sp1 = c(1, 1, 100)))
  expect_true(all(sim$abundance[1:2, 1] > 800))
  expect_equal(sim$abundance[3, 1], 0L)
})

test_that("lottery recruitment fills the microsites and emigration does not empty the source", {
  tr <- data.frame(optimum = 40, sigma = 35, r_max = 5, row.names = "sp1")
  ls <- landscape_matrix(matrix(1, 2, 2), env = c(40, 40))
  init <- cbind(sp1 = c(300L, 0L))
  set.seed(1)
  sim <- simulate_metacommunity(ls, tr, matrix(0.004, 1, 1), dispersal = 0.5, n_steps = 30,
                                seed_steps = integer(0), initial = init,
                                pool = cbind(sp1 = c(1, 0)), recruitment = "lottery",
                                microsites = 500, record_steps = 1:30)
  # half of the offspring are lost to the closed site, but they still outnumber the vacancies
  expect_true(all(sim$history[1, 1, 2:30] == 500))
  expect_equal(sim$recruitment, "lottery")
  # Beverton-Holt with the same loss stays well below (5 * 0.5 - 1) / 0.004 = 375
  set.seed(1)
  bh <- simulate_metacommunity(ls, tr, matrix(0.004, 1, 1), dispersal = 0.5, n_steps = 30,
                               seed_steps = integer(0), initial = init,
                               pool = cbind(sp1 = c(1, 0)))
  expect_lt(bh$abundance[1, 1], 450)
  # two equivalent species share the microsites; overlapping generations keep survivors
  tr2 <- data.frame(optimum = c(40, 40), sigma = 35, r_max = 5, row.names = c("a", "b"))
  set.seed(2)
  sim <- simulate_metacommunity(landscape_matrix(matrix(1, 3, 3), env = c(40, 40, 60)), tr2,
                                matrix(0.004, 2, 2), dispersal = c(0.05, 0.2), n_steps = 50,
                                seed_steps = integer(0), initial = cbind(a = c(100L, 0L, 0L),
                                                                         b = c(0L, 100L, 0L)),
                                recruitment = "lottery", microsites = c(200, 300, 400),
                                mortality = 0.2, immigration = "multinomial")
  expect_equal(unname(rowSums(sim$abundance)), c(200, 300, 400))
  expect_error(simulate_metacommunity(ls, tr, matrix(0.004, 1, 1), recruitment = "lottery"),
               "microsites")
  expect_error(simulate_metacommunity(ls, tr, matrix(0.004, 1, 1), recruitment = "lottery",
                                      microsites = 10, mortality = 0), "mortality")
})

test_that("lottery_draw fills vacancies without exceeding the offspring", {
  set.seed(3)
  pool <- matrix(c(5L, 0L, 30L, 2L, 0L, 10L, 1L, 0L, 50L), 3, 3)
  out <- lottery_draw(pool, c(4, 3, 100))
  expect_true(all(out <= pool))
  expect_equal(rowSums(out), c(4, 0, 90))
  expect_equal(out[3, ], pool[3, ])
})

test_that("lottery keeps the microsites, counts immigrant recruits and checks its settings", {
  tr <- data.frame(optimum = c(40, 40), sigma = 35, r_max = 5, row.names = c("a", "b"))
  ls <- landscape_points(cbind(1:12, 0), env = rep(40, 12), disp_scale = 2)
  set.seed(5)
  sim <- simulate_metacommunity(ls, tr, matrix(0.004, 2, 2), dispersal = 0.5, n_steps = 40,
                                seed_mean = 3, recruitment = "lottery", microsites = 10,
                                mortality = 0.1, record_steps = 1:40, disturbance = 0.05,
                                extirpation = 0.01,
                                evolution = evolution_control(genetic_var = 4, ne_ratio = 0.5))
  expect_true(all(apply(sim$history, c(1, 3), sum) <= 10))
  expect_true(all(sim$migrant_load <= 40, na.rm = TRUE))
  expect_true(all(sim$propagule <= 10, na.rm = TRUE))
  expect_true(is.integer(sim$abundance))
  expect_equal(sim$microsites, 10)
  expect_error(simulate_metacommunity(ls, tr, matrix(0.004, 2, 2), recruitment = "lottery",
                                      microsites = 5, seed_steps = integer(0),
                                      initial = cbind(a = rep(4L, 12), b = rep(4L, 12))),
               "more individuals")
  expect_error(simulate_metacommunity(ls, tr, matrix(0.004, 2, 2), recruitment = "lottery",
                                      microsites = 5,
                                      evolution = evolution_control(genetic_var = 4,
                                                                    steps_per_generation = 2)),
               "steps_per_generation")
  # pool_from with a lottery: the closed site stays empty, then fills
  set.seed(6)
  sim <- simulate_metacommunity(landscape_matrix(matrix(1, 3, 3), env = c(40, 40, 60)), tr,
                                matrix(0.004, 2, 2), dispersal = 0.1, n_steps = 30,
                                seed_steps = integer(0),
                                initial = cbind(a = c(50L, 0L, 0L), b = c(0L, 50L, 0L)),
                                pool_from = cbind(a = c(1, 1, 10), b = c(1, 1, 10)),
                                recruitment = "lottery", microsites = 50, record_steps = 1:30,
                                evolution = evolution_control(genetic_var = 4, ne_ratio = 0.5,
                                                              founder = FALSE))
  expect_true(all(sim$history[3, , 1:9] == 0))
  expect_equal(sum(sim$abundance[3, ]), 50)
})
