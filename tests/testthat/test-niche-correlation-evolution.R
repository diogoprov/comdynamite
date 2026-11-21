cor2 <- function(rho) matrix(c(1, rho, rho, 1), 2)

nc_traits <- function(optima = c(50, 50), sigma = c(10, 10), r_max = 5, species = "sp1") {
  k <- length(species)
  list(a = data.frame(optimum = rep(optima[1], k), sigma = sigma[1], r_max = r_max,
                      row.names = species),
       b = data.frame(optimum = rep(optima[2], k), sigma = sigma[2], r_max = 1,
                      row.names = species))
}

nc_state <- function(traits, env, G, E = matrix(0, 2, 2), ...) {
  ev <- evolution_control(genetic_var = G, env_var = E, ne_ratio = Inf, ...)
  ev$traits <- traits
  ev$env <- env
  evolution_state_mv(ev, nrow(env), nrow(traits[[1]]), rownames(traits[[1]]))
}

test_that("growth and selection use Omega = D R D", {
  sigma <- c(10, 4)
  rho <- 0.6
  tr <- niche_correlation(nc_traits(sigma = sigma), cor2(rho))
  env <- cbind(a = c(50, 58, 44), b = c(50, 53, 51))
  G <- matrix(c(6, 1, 1, 2), 2)
  E <- matrix(c(3, -0.5, -0.5, 1), 2)
  st <- nc_state(tr, env, G, E)
  expect_true(st$correlated)
  om <- diag(sigma) %*% cor2(rho) %*% diag(sigma)
  s <- st$sel[[1]]
  expect_equal(s$A, solve(om + G + E))
  expect_equal(s$K, G %*% solve(om + G + E))
  expect_equal(s$c, 5 * sqrt(det(om) / det(om + G + E)))
  # the growth of a population is the mean growth of its individuals
  pa <- array(rep(c(50, 50), each = 3), c(3, 1, 2))
  r <- mv_growth(st, pa, unname(env), 1L)
  d <- sweep(unname(env), 2, c(50, 50))
  expect_equal(r[, 1], s$c * exp(-0.5 * rowSums((d %*% solve(om + G + E)) * d)))
  set.seed(1)
  z <- matrix(rnorm(4e5), ncol = 2) %*% chol(G + E)
  ind <- vapply(seq_len(3), function(i) {
    dz <- sweep(-z, 2, d[i, ], "+")
    mean(5 * exp(-0.5 * rowSums((dz %*% solve(om)) * dz)))
  }, numeric(1))
  expect_equal(r[, 1], ind, tolerance = 0.02)
  # the response to selection is G (Omega + P)^-1 d, and immigrants are filtered alike
  za <- mv_select(st, pa, pa, unname(env), matrix(10L, 3, 1), 1L)
  expect_equal(matrix(za, 3) - 50, d %*% t(G %*% solve(om + G + E)))
  w <- mv_fitness(st, matrix(50, 1, 2), unname(env), 1L)
  expect_equal(w[1, ], exp(-0.5 * rowSums((d %*% solve(om + G + E)) * d)))
  # without variance, the growth of growth_rates() with the same correlation
  st0 <- nc_state(tr, env, 0)
  expect_equal(mv_growth(st0, pa, unname(env), 1L)[, 1], unname(growth_rates(env, tr)[, 1]))
  # an axis with an infinite breadth exerts no selection, as in growth_rates()
  tri <- niche_correlation(nc_traits(sigma = c(10, Inf)), cor2(rho))
  sti <- nc_state(tri, env, 0)
  expect_equal(mv_growth(sti, pa, unname(env), 1L)[, 1], unname(growth_rates(env, tri)[, 1]))
  expect_equal(sti$sel[[1]]$s, 1L)
})

test_that("a correlated response moves the trait of the other axis", {
  tr <- nc_traits()
  env <- cbind(a = 60, b = 50)
  pa <- array(50, c(1, 1, 2))
  G <- diag(2) * 4
  step <- function(rho) {
    st <- nc_state(niche_correlation(tr, cor2(rho)), env, G)
    drop(mv_select(st, pa, pa, unname(env), matrix(10L, 1, 1), 1L)) - 50
  }
  expect_gt(step(0)[1], 0)
  expect_equal(step(0)[2], 0)
  # with a positive correlation, a departure on the first axis selects the second trait away
  # from its own environment
  expect_lt(step(0.7)[2], 0)
  expect_gt(step(-0.7)[2], 0)
})

test_that("the identity gives the results of independent axes in simulate_metacommunity", {
  set.seed(6)
  l <- landscape_grid(side = 5, n_axes = 2, axis_cor = 0.5)
  tr <- list(species_traits(3), species_traits(3, rmax_range = c(1, 1)))
  a <- competition_matrix(3)
  ev <- evolution_control(genetic_var = matrix(c(5, 2, 2, 5), 2), env_var = diag(2),
                          plasticity = 0.2, immigrant_selection = TRUE)
  run <- function(traits) {
    set.seed(7)
    simulate_metacommunity(l, traits, a, dispersal = 0.1, n_steps = 60, evolution = ev,
                           record_steps = c(30, 60))
  }
  base <- run(tr)
  shared <- run(niche_correlation(tr, diag(2)))
  each <- run(niche_correlation(tr, stats::setNames(rep(list(diag(2)), 3), rownames(tr[[1]]))))
  expect_gt(sum(base$abundance), 0)
  for (other in list(shared, each)) {
    expect_identical(other$abundance, base$abundance)
    expect_identical(other$history, base$history)
    expect_identical(other$evolution, base$evolution)
    expect_identical(other$migrant_load, base$migrant_load)
  }
  expect_identical(attr(shared$traits, "niche_cor"), diag(2))
  cor <- run(niche_correlation(tr, cor2(0.8)))
  expect_false(identical(cor$abundance, base$abundance))
  expect_false(isTRUE(all.equal(cor$evolution$genetic, base$evolution$genetic)))
  # without genetic variance, the abundances of the same niche without evolution
  ev0 <- evolution_control(genetic_var = 0, ne_ratio = Inf)
  trc <- niche_correlation(tr, cor2(0.8))
  set.seed(8)
  s1 <- simulate_metacommunity(l, trc, a, dispersal = 0.1, n_steps = 30, evolution = ev0)
  set.seed(8)
  s2 <- simulate_metacommunity(l, trc, a, dispersal = 0.1, n_steps = 30)
  expect_equal(s1$abundance, s2$abundance)
})

test_that("departures of the same sign do better with a positive correlation", {
  l <- landscape_grid(side = 4, disp_scale = 0.5)
  env <- cbind(a = rep(65, 16), b = rep(65, 16))
  tr <- nc_traits()
  a <- matrix(0.05, 1, 1, dimnames = list("sp1", "sp1"))
  ev <- evolution_control(genetic_var = diag(2) * 0.01, ne_ratio = Inf)
  start <- matrix(30L, 16, 1, dimnames = list(NULL, "sp1"))
  run <- function(rho) {
    set.seed(3)
    sim <- simulate_metacommunity(l, niche_correlation(tr, cor2(rho)), a, dispersal = 0.01,
                                  n_steps = 60, env = env, evolution = ev, initial = start,
                                  seed_steps = numeric(0))
    list(n = sum(sim$abundance), sim = sim)
  }
  # standardised departures (1.5, 1.5): growth 5 exp(-2.25 / (1 + rho))
  pos <- run(0.8)
  expect_gt(pos$n, 30)
  expect_equal(run(0)$n, 0)
  expect_equal(run(-0.5)$n, 0)
  # the load counts the correlation: departures of the same sign weigh less
  m <- maladaptation(pos$sim)
  om <- diag(10, 2) %*% cor2(0.8) %*% diag(10, 2) + diag(2) * 0.01
  d <- as.matrix(m[1, c("deviation_a", "deviation_b")])
  expect_equal(m$distance[1], sqrt(sum((d %*% solve(om)) * d)))
  expect_lt(m$distance[1], sqrt(sum(d^2 / (100 + 0.01))))
})

nc_evolving_setup <- function(side = 6) {
  ls <- landscape_grid(side, disp_scale = 0.5)
  list(ls = ls, env = cbind(svl = ls$env, head = seq(-3, 3, length.out = side^2)),
       anc = list(svl = data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1"),
                  head = data.frame(optimum = 0, sigma = 3, r_max = 1, row.names = "sp1")))
}

test_that("the identity gives the results of independent axes with speciation", {
  s <- nc_evolving_setup()
  ev <- evolution_control(matrix(c(10, 2, 2, 1), 2), env_var = diag(2),
                          immigrant_selection = TRUE)
  run <- function(traits, seed = 2) {
    set.seed(seed)
    simulate_evolving_metacommunity(s$ls, traits, speciation_control("fission", rate = 0.02),
                                    evolution = ev, env = s$env, dispersal = 0.005,
                                    n_steps = 120, record_steps = c(60, 120))
  }
  base <- run(s$anc)
  same <- run(niche_correlation(s$anc, diag(2)))
  expect_gt(ncol(base$abundance), 1)
  expect_identical(attr(same$traits, "niche_cor"), diag(2))
  attr(same$traits, "niche_cor") <- NULL
  expect_identical(same, base)
  # a shared correlation applies to the new species and changes the outcome
  cor <- run(niche_correlation(s$anc, cor2(0.8)))
  expect_gt(nrow(cor$species), 1)
  expect_identical(attr(cor$traits, "niche_cor"), cor2(0.8))
  expect_false(identical(cor$lineages, base$lineages))
  m <- maladaptation(cor)
  expect_true(all(m$load >= 0 & m$load < 1))
  expect_s3_class(as_metacommunity(cor), "mc_sim")
  expect_identical(attr(as_metacommunity(cor)$traits, "niche_cor"), cor2(0.8))
})

test_that("daughters inherit the correlation of their mother", {
  s <- nc_evolving_setup()
  tree <- ape::read.tree(text = "(a:50,b:50);")
  anc <- lapply(s$anc, function(t) {
    t <- t[c(1, 1), ]
    rownames(t) <- c("a", "b")
    t
  })
  attr(anc[[1]], "tree") <- tree
  cors <- list(a = cor2(0.7), b = cor2(-0.4))
  set.seed(5)
  sim <- simulate_evolving_metacommunity(
    s$ls, niche_correlation(anc, cors), speciation_control("fission", rate = 0.03),
    evolution = evolution_control(matrix(c(10, 2, 2, 1), 2), env_var = diag(2)),
    env = s$env, dispersal = 0.005, n_steps = 120)
  out <- attr(sim$traits, "niche_cor")
  extant <- colnames(sim$abundance)
  expect_named(out, extant)
  expect_gt(length(setdiff(extant, c("a", "b"))), 0)
  # the ancestor at the start of each extant species
  tab <- sim$species
  root <- vapply(extant, function(sp) {
    while (!is.na(tab$parent[tab$species == sp])) sp <- tab$parent[tab$species == sp]
    sp
  }, character(1))
  expect_true(all(c("a", "b") %in% root))
  for (sp in extant) expect_identical(out[[sp]], cors[[root[[sp]]]])
  m <- maladaptation(sim)
  expect_true(all(is.finite(m$distance)))
})

test_that("niche overlap counts the correlation", {
  # crowding with identity matrices is that of independent axes
  st <- new.env()
  st$alive <- c(1L, 3L, 4L)
  st$at <- matrix(c(1, 0.5, 0.2, 0.4, 1, 0.3, 0.6, 0.7, 1), 3) * 0.05
  sigma <- cbind(c(10, 8, 5, 12), c(3, Inf, 2, 4))
  st$ev <- list(mv = TRUE, q = 2, sigma = sigma, sel = rep(list(list(cor = diag(2))), 4))
  set.seed(1)
  n <- matrix(rpois(15, 3), 5, 3)
  p <- array(rnorm(30, 0, 5), c(5, 3, 2))
  old <- niche_crowding(st)(n, p)
  new <- niche_crowding_correlated(st)(n, p)
  expect_equal(new[n > 0], old[n > 0])
  # competition_matrix(type = "trait")
  tr <- list(data.frame(optimum = c(10, 20, 35), sigma = c(5, 8, 6), r_max = 5,
                        row.names = c("x", "y", "z")),
             data.frame(optimum = c(1, 3, 2), sigma = c(1, 2, 1.5), r_max = 1,
                        row.names = c("x", "y", "z")))
  base <- competition_matrix(type = "trait", traits = tr, inter_max = 1)
  expect_equal(competition_matrix(type = "trait", traits = niche_correlation(tr, diag(2)),
                                  inter_max = 1), base)
  rho <- 0.7
  a <- competition_matrix(type = "trait", traits = niche_correlation(tr, cor2(rho)),
                          inter_max = 1)
  om <- function(j) diag(c(tr[[1]]$sigma[j], tr[[2]]$sigma[j])) %*% cor2(rho) %*%
    diag(c(tr[[1]]$sigma[j], tr[[2]]$sigma[j]))
  d <- c(10, 2)
  expect_equal(a["x", "y"], 0.05 * exp(-0.5 * drop(d %*% solve(om(1) + om(2)) %*% d)))
  expect_true(isSymmetric(a))
  # x and y differ in the same direction on both axes: more overlap with a positive correlation
  expect_gt(a["x", "y"], base["x", "y"])
})

test_that("invalid correlated niches are refused", {
  tr <- nc_traits()
  expect_error(niche_correlation(tr, cor2(1.2)), "positive definite")
  expect_error(niche_correlation(tr, list(other = diag(2))), "one matrix per species")
  l <- landscape_grid(side = 3, n_axes = 2)
  env <- l$env_axes
  colnames(env) <- c("a", "b")
  a <- matrix(0.05, 1, 1, dimnames = list("sp1", "sp1"))
  ev <- evolution_control(genetic_var = diag(2), env_var = diag(2))
  # a breadth of 0 has no correlated niche
  zero <- niche_correlation(nc_traits(sigma = c(0, 10)), cor2(0.5))
  expect_error(simulate_metacommunity(l, zero, a, n_steps = 2, evolution = ev, env = env),
               "correlated axes")
  expect_error(simulate_evolving_metacommunity(l, zero, speciation_control(threshold = 50),
                                               evolution = ev, n_steps = 2, env = env),
               "correlated axes")
  expect_error(competition_matrix(type = "trait", traits = zero), "correlated axes")
  # an attribute that did not come from niche_correlation()
  bad <- tr
  attr(bad, "niche_cor") <- diag(3)
  expect_error(simulate_metacommunity(l, bad, a, n_steps = 2, evolution = ev, env = env),
               "2 x 2")
  attr(bad, "niche_cor") <- list(other = diag(2))
  expect_error(simulate_evolving_metacommunity(l, bad, speciation_control(threshold = 50),
                                               evolution = ev, n_steps = 2, env = env),
               "lacks matrices")
  # mutation_var is still refused with several traits
  expect_error(simulate_metacommunity(l, niche_correlation(tr, cor2(0.5)), a, n_steps = 2,
                                      evolution = evolution_control(1, mutation_var = 0.1),
                                      env = env),
               "mutation_var")
})
