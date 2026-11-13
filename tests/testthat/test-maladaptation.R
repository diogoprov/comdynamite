mal_setup <- function() {
  land <- landscape_grid(6, disp_scale = 0.5)
  list(land = land,
       tr = data.frame(optimum = c(45, 55), sigma = c(15, 20), r_max = 4,
                       row.names = c("sp1", "sp2")),
       alpha = matrix(0.05, 2, 2, dimnames = list(c("sp1", "sp2"), c("sp1", "sp2"))))
}

test_that("maladaptation() gives the lag load of every population", {
  set.seed(1)
  s <- mal_setup()
  fixed <- simulate_metacommunity(s$land, s$tr, s$alpha, dispersal = 0.05, n_steps = 60,
                                  evolution = evolution_control(0, env_var = 2),
                                  record_steps = c(30, 60))
  m <- maladaptation(fixed)
  expect_named(m, c("site", "species", "abundance", "env", "phenotype", "deviation",
                    "distance", "load"))
  expect_true(all(m$abundance > 0))
  theta <- s$tr[m$species, "optimum"]
  width <- s$tr[m$species, "sigma"]^2 + 2
  expect_equal(m$phenotype, theta)
  expect_equal(m$load, 1 - exp(-(m$env - theta)^2 / (2 * width)))
  expect_equal(nrow(m), sum(fixed$abundance > 0))
  both <- maladaptation(fixed, step = c(30, 60))
  expect_identical(sort(unique(both$step)), c(30, 60))
  expect_equal(both[both$step == 60, -1], m, ignore_attr = TRUE)
  set.seed(1)
  evolving <- simulate_metacommunity(s$land, s$tr, s$alpha, dispersal = 0.05, n_steps = 200,
                                     evolution = evolution_control(4, env_var = 2))
  expect_lt(mean(maladaptation(evolving)$load), mean(m$load))
  # a niche of infinite breadth carries no load
  wide <- s$tr
  wide$sigma <- Inf
  w <- simulate_metacommunity(s$land, wide, s$alpha, dispersal = 0.05, n_steps = 10,
                              evolution = evolution_control(1))
  expect_true(all(maladaptation(w)$load == 0))
  expect_error(maladaptation(simulate_metacommunity(s$land, s$tr, s$alpha, n_steps = 5)),
               "with `evolution`")
  expect_error(maladaptation(fixed, step = 10), "not recorded")
})

test_that("maladaptation() uses the Mahalanobis distance with several traits", {
  set.seed(2)
  land <- landscape_grid(5, disp_scale = 0.5)
  tr <- list(svl = data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1"),
             head = data.frame(optimum = 0, sigma = 4, r_max = 1, row.names = "sp1"))
  G <- matrix(c(10, 6, 6, 8), 2)
  E <- diag(2, 2)
  sim <- simulate_metacommunity(land, tr, matrix(0.05, dimnames = list("sp1", "sp1")),
                                env = cbind(svl = land$env, head = seq(-3, 3, length.out = 25)),
                                dispersal = 0.05, n_steps = 50, record_steps = c(25, 50),
                                evolution = evolution_control(G, env_var = E))
  m <- maladaptation(sim)
  expect_true(all(c("env_svl", "phenotype_svl", "deviation_svl", "env_head",
                    "phenotype_head", "deviation_head", "distance", "load") %in% names(m)))
  d <- as.matrix(m[, c("deviation_svl", "deviation_head")])
  A <- solve(diag(c(400, 16)) + G + E)
  expect_equal(m$distance, sqrt(rowSums((d %*% A) * d)))
  expect_equal(m$load, 1 - exp(-m$distance^2 / 2))
  expect_equal(m$deviation_svl, m$env_svl - m$phenotype_svl)
  both <- maladaptation(sim, step = c(25, 50))
  expect_equal(both[both$step == 50, -1], m, ignore_attr = TRUE)
  tr$head$sigma <- Inf
  sim2 <- simulate_metacommunity(land, tr, matrix(0.05, dimnames = list("sp1", "sp1")),
                                 env = cbind(svl = land$env, head = 0), dispersal = 0.05,
                                 n_steps = 20, evolution = evolution_control(G))
  m2 <- maladaptation(sim2)
  expect_equal(m2$distance, abs(m2$deviation_svl) / sqrt(400 + 10))
})

test_that("maladaptation() takes evolving metacommunities", {
  set.seed(3)
  land <- landscape_grid(5, disp_scale = 0.5)
  anc <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
  sim <- simulate_evolving_metacommunity(land, anc, speciation_control("point", rate = 1e-3),
                                         evolution = evolution_control(2), n_steps = 60)
  m <- maladaptation(sim)
  expect_true(nrow(m) > 0)
  expect_true(all(m$load >= 0 & m$load < 1))
})

test_that("maladaptation() reads recorded steps, plasticity and a changing environment", {
  set.seed(4)
  s <- mal_setup()
  env_time <- sapply(1:40, function(t) s$land$env + 0.5 * t)
  sim <- simulate_metacommunity(s$land, s$tr, s$alpha, dispersal = 0.05, n_steps = 40,
                                env_time = env_time, record_steps = c(20, 30),
                                evolution = evolution_control(3, env_var = 1,
                                                              plasticity = 0.4))
  m20 <- maladaptation(sim, step = 20)
  ev <- sim$evolution
  z <- ev$history[, , "step20"]
  i <- match(m20$site, rownames(z))
  j <- match(m20$species, colnames(z))
  expected <- z[cbind(i, j)] + 0.4 * (env_time[i, 20] - ev$reference)
  expect_equal(m20$env, env_time[i, 20])
  expect_equal(m20$phenotype, expected)
  width <- s$tr[m20$species, "sigma"]^2 + 3 + 1
  expect_equal(m20$load, 1 - exp(-(env_time[i, 20] - expected)^2 / (2 * width)))
  # the last step falls back to the final state when it was not recorded
  expect_equal(maladaptation(sim, step = 40), maladaptation(sim))
  expect_error(maladaptation(sim, step = numeric(0)), "recorded steps")
  expect_error(maladaptation(sim, step = NA), "recorded steps")
})

test_that("a niche of infinite breadth exerts no selection", {
  r <- evolution_growth(matrix(c(10, 20), 1), matrix(c(0, 50), 1), sigma = c(Inf, 5),
                        r_max = c(3, 2), v = 4)
  expect_identical(r[1, 1], 3)
  expect_equal(r[1, 2], 2 * sqrt(25 / 29) * exp(-30^2 / 58))
  s <- mal_setup()
  wide <- s$tr
  wide$sigma <- Inf
  sim <- simulate_metacommunity(s$land, wide, s$alpha, dispersal = 0, n_steps = 30,
                                evolution = evolution_control(5, ne_ratio = Inf,
                                                              founder = FALSE))
  z <- sim$evolution$genetic
  expect_equal(unname(z[, "sp1"][!is.na(z[, "sp1"])]), rep(45, sum(!is.na(z[, "sp1"]))))
})
