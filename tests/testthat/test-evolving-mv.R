ev_mv_setup <- function(side = 5) {
  ls <- landscape_grid(side, disp_scale = 0.5)
  list(ls = ls, env = cbind(svl = ls$env, head = seq(-3, 3, length.out = side^2)),
       anc = list(svl = data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1"),
                  head = data.frame(optimum = 0, sigma = 3, r_max = 1, row.names = "sp1")))
}

test_that("several traits evolve and speciate in evolving metacommunities", {
  s <- ev_mv_setup(6)
  G <- matrix(c(10, 2, 2, 1), 2)
  set.seed(2)
  sim <- simulate_evolving_metacommunity(s$ls, s$anc, speciation_control("fission", rate = 0.02),
                                         evolution = evolution_control(G, env_var = diag(2)),
                                         env = s$env, dispersal = 0.005, n_steps = 150,
                                         record_steps = c(75, 150))
  k <- ncol(sim$abundance)
  expect_gt(k, 1)
  expect_equal(dim(sim$genetic), c(36, ncol(sim$lineages), 2))
  expect_named(sim$traits, c("svl", "head"))
  expect_identical(rownames(sim$traits$svl), colnames(sim$abundance))
  expect_true(all(c("genetic_svl", "genetic_head") %in% names(sim$series)))
  ev <- sim$evolution
  expect_equal(dim(ev$genetic), c(36, k, 2))
  expect_equal(dim(ev$history), c(36, k, 2, 2))
  expect_equal(dim(ev$environment_history), c(36, 2, 2))
  expect_length(ev$genetic_var, k)
  expect_equal(unname(ev$genetic_var[[1]]), G)
  expect_true(all(is.na(ev$genetic[, , 1]) == (sim$abundance == 0)))
  # the series ends with the genetic means of the lineages, weighted by abundance
  last <- sim$series[sim$series$generation == 150 & sim$series$abundance > 0, ]
  l1 <- last$lineage[1]
  w <- sim$lineages[, l1]
  expect_equal(last$genetic_svl[1], sum(w * sim$genetic[, l1, "svl"], na.rm = TRUE) / sum(w))
  mc <- as_metacommunity(sim)
  expect_s3_class(mc, "mc_sim")
  expect_identical(names(mc$env), c("svl", "head"))
  ind <- sample_individuals(sim, n = 2, source = "evolution", step = c(75, 150))
  expect_true(all(c("svl", "head", "step") %in% names(ind)))
  m <- maladaptation(sim, step = 75)
  expect_true(all(m$load >= 0 & m$load < 1))
  expect_output(print(speciation_control("point", rate = 1e-4, jump_var = diag(2))),
                "2 x 2 matrix")
})

test_that("jumps follow jump_var and species without genetic variance keep their optimum", {
  s <- ev_mv_setup(6)
  set.seed(4)
  # a jump of rank one moves both traits by the same amount
  J <- matrix(1, 2, 2)
  sim <- simulate_evolving_metacommunity(s$ls, s$anc,
                                         speciation_control("fission", rate = 0.03,
                                                            jump_var = J),
                                         evolution = evolution_control(matrix(0, 2, 2)),
                                         env = s$env, dispersal = 0.005, n_steps = 120)
  tab <- sim$species
  both <- tab$species %in% colnames(sim$abundance) &
    tab$parent %in% colnames(sim$abundance)
  expect_true(any(both))
  child <- tab$species[both]
  parent <- tab$parent[both]
  d_svl <- sim$traits$svl[child, "optimum"] - sim$traits$svl[parent, "optimum"]
  d_head <- sim$traits$head[child, "optimum"] - sim$traits$head[parent, "optimum"]
  expect_equal(d_svl, d_head)
  expect_true(any(abs(d_svl) > 0))
  # without jumps the optimum stays that of the ancestor
  set.seed(4)
  still <- simulate_evolving_metacommunity(s$ls, s$anc, speciation_control("fission", rate = 0.03),
                                           evolution = evolution_control(matrix(0, 2, 2)),
                                           env = s$env, dispersal = 0.005, n_steps = 60)
  expect_true(all(still$traits$svl$optimum == 50))
  expect_true(all(still$traits$head$optimum == 0))
})

test_that("ecological speciation with several traits uses the Mahalanobis distance", {
  ev <- list(whiten = list(evolution_whiten(matrix(c(4, 1, 1, 2), 2))))
  z <- rbind(c(0, 0), c(2, 1), c(-1, 3))
  P <- matrix(c(4, 1, 1, 2), 2)
  manual <- outer(1:3, 1:3, Vectorize(function(a, b) {
    d <- z[a, ] - z[b, ]
    sqrt(drop(t(d) %*% solve(P) %*% d))
  }))
  expect_equal(unname(mv_distance(ev, 1, z)), manual)
  # a singular P leaves out the direction without variance
  ev0 <- list(whiten = list(evolution_whiten(diag(c(4, 0)))))
  expect_equal(unname(mv_distance(ev0, 1, rbind(c(0, 0), c(2, 9))))[1, 2], 1)
  s <- ev_mv_setup(6)
  set.seed(1)
  sim <- simulate_evolving_metacommunity(s$ls, s$anc,
                                         speciation_control("ecological", threshold = 1,
                                                            min_generations = 5,
                                                            migration_max = 0.01),
                                         evolution = evolution_control(diag(c(20, 2)),
                                                                       env_var = diag(2)),
                                         env = s$env, dispersal = 0.0005, n_steps = 200)
  expect_gt(nrow(sim$events), 0)
})

test_that("niche overlap on several axes is the product over axes", {
  st <- new.env()
  st$ev <- list(mv = TRUE, q = 2, sigma = rbind(c(2, 3), c(4, Inf)))
  st$alive <- 1:2
  st$at <- matrix(c(1, 0.5, 0.5, 1), 2)
  hook <- niche_crowding(st)
  n <- matrix(c(10, 5), 1)
  p <- array(c(1, 2, 0, 5), c(1, 2, 2))
  overlap <- exp(-(1 - 2)^2 / (2 * (4 + 16)))
  expect_equal(hook(n, p)[1, ], c(10 + 5 * 0.5 * overlap, 10 * 0.5 * overlap + 5))
})

test_that("evolving metacommunities check settings for several traits", {
  s <- ev_mv_setup()
  sc <- speciation_control("fission", rate = 0.01)
  expect_error(simulate_evolving_metacommunity(s$ls, s$anc, sc), "needs `env`")
  expect_error(simulate_evolving_metacommunity(s$ls, s$anc$svl, sc, env = s$env),
               "list of traits")
  expect_error(simulate_evolving_metacommunity(s$ls, s$anc, sc, env = s$env[1:5, ]),
               "one row per site")
  bad <- s$anc
  rownames(bad$head) <- "spx"
  expect_error(simulate_evolving_metacommunity(s$ls, bad, sc, env = s$env), "different species")
  expect_error(simulate_evolving_metacommunity(s$ls, s$anc$svl, sc,
                                               evolution = evolution_control(diag(2))),
               "several axes")
  expect_error(simulate_evolving_metacommunity(s$ls, s$anc$svl,
                                               speciation_control("point", rate = 1e-4,
                                                                  jump_var = diag(2))),
               "matrix `jump_var`")
  expect_error(simulate_evolving_metacommunity(s$ls, s$anc,
                                               speciation_control("point", rate = 1e-4,
                                                                  jump_var = diag(3)),
                                               env = s$env), "jump_var")
})
