# a genetic architecture of a few loci (`loci` of evolution_control()): the allele
# frequencies against results of population genetics with known answers

one_sp <- function(n, m) matrix(n, m, 1, dimnames = list(NULL, "sp1"))
loci_traits <- function(optimum, sigma = 1e6, r_max = 5) {
  data.frame(optimum = optimum, sigma = sigma, r_max = r_max, row.names = "sp1")
}
# two species on a small grid with two kinds of sites, as in the replication of Leibold et al.
loci_grid <- function(side = 4) {
  xy <- expand.grid(x = seq_len(side) - 1, y = seq_len(side) - 1)
  env <- rep(c(40, 60), length.out = side^2)
  sp <- c("sp1", "sp2")
  list(ls = landscape_points(xy, env, disp_scale = 0.5), env = env,
       traits = data.frame(optimum = c(40, 60), sigma = 20, r_max = 5, row.names = sp),
       alpha = matrix(0.01, 2, 2, dimnames = list(sp, sp)),
       initial = cbind(sp1 = (env == 40) * 100L, sp2 = (env == 60) * 100L))
}
# frequencies between 0 and 1, NA where absent, and a mean and a variance that follow from them
expect_loci_state <- function(sim, loci) {
  e <- sim$evolution
  p <- e$allele_frequency
  absent <- sim$abundance == 0
  expect_identical(dim(p), c(dim(sim$abundance), as.integer(loci$n)))
  expect_identical(dimnames(p)[1:2], dimnames(sim$abundance))
  expect_true(all(is.na(p[rep(absent, loci$n)])))
  expect_false(anyNA(p[rep(!absent, loci$n)]))
  expect_true(all(p >= 0 & p <= 1, na.rm = TRUE))
  origin <- if (is.null(loci$origin)) 0 else loci$origin
  expect_equal(e$genetic, origin + 2 * loci$effect * apply(p, 1:2, sum))
  expect_equal(e$genetic_variance, 2 * loci$effect^2 * apply(p * (1 - p), 1:2, sum))
  expect_null(e$linkage_disequilibrium)
}

test_that("neutral drift loses 1 / (2 Ne) of the heterozygosity per generation", {
  set.seed(11)
  m <- 2000
  steps <- 20
  # one locus at p = 0.5: the optimum is half of the range
  ev <- evolution_control(ne_ratio = 0.5, env = rep(0, m), traits = loci_traits(1),
                          loci = list(n = 1, effect = 1, mutation = 0))
  out <- mc_dynamics(matrix(5, m, 1), matrix(0.05), dispersal = 0, n_steps = steps,
                     seed_steps = numeric(0), initial = one_sp(80L, m),
                     record_steps = seq_len(steps), evolution = ev)
  n <- out$history[, 1, ]
  ok <- apply(n > 0, 1, all)
  p <- out$evolution$allele_frequency[ok, 1, 1]
  # 2 Ne = 2 * ne_ratio * N alleles are sampled each generation
  expected <- 0.5 * mean(apply(1 - 1 / n[ok, ], 1, prod))
  expect_equal(mean(2 * p * (1 - p)), expected, tolerance = 0.05)
  expect_lt(expected, 0.42)
  # without drift the frequencies stay where they were
  ev$ne_ratio <- Inf
  ev$env <- rep(0, 50)
  still <- mc_dynamics(matrix(5, 50, 1), matrix(0.05), dispersal = 0, n_steps = 5,
                       seed_steps = numeric(0), initial = one_sp(80L, 50), evolution = ev)
  expect_equal(as.vector(still$evolution$allele_frequency), rep(0.5, 50), tolerance = 1e-9)
})

test_that("mutation and drift balance at the heterozygosity of the two-allele model", {
  set.seed(12)
  m <- 500
  mu <- 0.005
  # five loci without selection are five replicates; they start fixed or at p = 0.5
  loci <- list(n = 5, effect = 1, mutation = mu)
  ev <- evolution_control(ne_ratio = 0.5, env = rep(0, m), traits = loci_traits(5),
                          loci = loci)
  out <- mc_dynamics(matrix(5, m, 1), matrix(0.05), dispersal = 0, n_steps = 300,
                     seed_steps = numeric(0), initial = one_sp(80L, m), evolution = ev)
  present <- out$abundance[, 1] > 0
  p <- out$evolution$allele_frequency[present, 1, ]
  ne <- 0.5 * mean(out$abundance[present, 1])
  # symmetric mutation between two alleles and binomial sampling of 2 Ne alleles (standard
  # population genetics): H' = (1 - 1 / (2 Ne)) (1 - (1 - 2 mu)^2 (1 - 2 H)) / 2, with the
  # equilibrium below, close to 4 Ne mu / (1 + 8 Ne mu) for a small mu
  keep <- 1 - 1 / (2 * ne)
  a <- (1 - 2 * mu)^2
  exact <- keep * (1 - a) / 2 / (1 - keep * a)
  expect_equal(exact, 4 * ne * mu / (1 + 8 * ne * mu), tolerance = 0.03)
  expect_equal(mean(2 * p * (1 - p)), exact, tolerance = 0.08)
  # mutation is symmetric: the mean frequency goes to one half
  expect_equal(mean(p), 0.5, tolerance = 0.05)
})

test_that("selection changes the frequencies by Wright's formula for the mean fitness", {
  set.seed(13)
  env <- c(3, 5, 9, 14)
  loci <- list(n = 4, effect = 2, mutation = 0, origin = -3)
  # log of the mean fitness of a population with the frequency p of the second locus (the
  # first is fixed), in the environment e, with the plastic part `plastic` of the phenotype
  log_w <- function(p, e, plastic) {
    z <- -3 + 2 * 2 * (1 + p)
    v <- 9 + 2 * 4 * p * (1 - p) + 1
    -log(v) / 2 - (e - z - plastic)^2 / (2 * v)
  }
  # c = (2 + 3) / 4 = 1.25: one locus fixed and one at p = 0.25
  for (b in c(0, 0.5)) {
    ev <- evolution_control(ne_ratio = Inf, env_var = 1, plasticity = b, env = env,
                            traits = loci_traits(2, sigma = 3), loci = loci)
    out <- mc_dynamics(matrix(5, 4, 1), matrix(0.001), dispersal = 0, n_steps = 1,
                       seed_steps = numeric(0), initial = one_sp(500L, 4), evolution = ev)
    plastic <- b * (env - mean(env))
    h <- 1e-6
    slope <- (log_w(0.25 + h, env, plastic) - log_w(0.25 - h, env, plastic)) / (2 * h)
    expect_equal(unname(out$evolution$allele_frequency[, 1, 2]) - 0.25,
                 0.25 * 0.75 / 2 * slope, tolerance = 1e-7)
    expect_identical(unname(out$evolution$allele_frequency[, 1, c(1, 3, 4)]),
                     matrix(c(1, 0, 0), 4, 3, byrow = TRUE))
    expect_loci_state(out, loci)
  }
  # selection per step is divided by the steps per generation
  half <- evolution_control(ne_ratio = Inf, env_var = 1, steps_per_generation = 2, env = env,
                            traits = loci_traits(2, sigma = 3), loci = loci)
  out2 <- mc_dynamics(matrix(5, 4, 1), matrix(0.001), dispersal = 0, n_steps = 1,
                      seed_steps = numeric(0), initial = one_sp(500L, 4), evolution = half)
  slope <- (log_w(0.25 + 1e-6, env, 0) - log_w(0.25 - 1e-6, env, 0)) / 2e-6
  expect_equal(unname(out2$evolution$allele_frequency[, 1, 2]) - 0.25,
               0.25 * 0.75 / 2 * slope / 2, tolerance = 1e-7)
  # frequencies stay in [0, 1], and only the populations asked for change
  p <- array(c(0.01, 0.99), c(1, 1, 2))
  far <- function(d) as.vector(loci_select(p, matrix(TRUE), matrix(d), matrix(1), 1, loci))
  # this far from the optimum, the second term rules and favours the rarer allele
  expect_identical(far(-90), c(1, 0))
  expect_identical(far(90), c(1, 0))
  expect_identical(loci_select(p, matrix(FALSE), matrix(1), matrix(1), 1, loci), p)
})

test_that("with the polymorphic loci at one half, the mean changes as Lande's equation", {
  set.seed(20)
  env <- c(3, 5, 9)
  loci <- list(n = 4, effect = 2, mutation = 0, origin = -3)
  # c = (3 + 3) / 4 = 1.5: one locus fixed and one at p = 0.5, G = 2 * 4 * 0.25
  for (b in c(0, 0.5)) {
    ev <- evolution_control(ne_ratio = Inf, env_var = 1, plasticity = b, env = env,
                            traits = loci_traits(3, sigma = 3), loci = loci)
    out <- mc_dynamics(matrix(5, 3, 1), matrix(0.001), dispersal = 0, n_steps = 1,
                       seed_steps = numeric(0), initial = one_sp(500L, 3), evolution = ev)
    phen <- 3 + b * (env - mean(env))
    expect_equal(unname(out$evolution$genetic[, 1]), 3 + 2 * (env - phen) / (9 + 2 + 1),
                 tolerance = 1e-12)
  }
})

test_that("stabilising selection at the optimum removes the variance of the loci", {
  set.seed(21)
  loci <- list(n = 4, effect = 2, mutation = 0, origin = -3)
  # in the engine: one locus at p = 0.25, in the environment of the optimum
  ev <- evolution_control(ne_ratio = Inf, env_var = 1, env = c(2, 2),
                          traits = loci_traits(2, sigma = 3), loci = loci)
  run <- function(steps) {
    mc_dynamics(matrix(5, 2, 1), matrix(0.001), dispersal = 0, n_steps = steps,
                seed_steps = numeric(0), initial = one_sp(500L, 2), evolution = ev)$evolution
  }
  one <- run(1)
  many <- run(60)
  expect_lt(one$allele_frequency[1, 1, 2], 0.25)
  expect_lt(one$genetic_variance[1, 1], 1.5)
  expect_lt(many$allele_frequency[1, 1, 2], one$allele_frequency[1, 1, 2])
  expect_lt(many$genetic_variance[1, 1], 0.5 * one$genetic_variance[1, 1])
  # two loci on either side of one half: both move towards fixation and the mean stays
  p <- p0 <- array(c(0.3, 0.7), c(1, 1, 2))
  two <- list(n = 2, effect = 2, mutation = 0, origin = 0)
  for (i in 1:40) {
    width <- 9 + loci_variance(p, two) + 1
    p <- loci_select(p, matrix(TRUE), 4 - loci_mean(p, two), width, 1, two)
  }
  expect_lt(p[1], 0.2)
  expect_gt(p[2], 0.8)
  expect_equal(loci_mean(p, two)[1, 1], 4)
  expect_lt(loci_variance(p, two)[1, 1], 0.75 * loci_variance(p0, two)[1, 1])
  # far from the optimum (d^2 > V) the second term favours variance instead
  rare <- array(0.1, c(1, 1, 1))
  one_locus <- list(n = 1, effect = 0.5, mutation = 0, origin = 0)
  far <- loci_select(rare, matrix(TRUE), matrix(10), matrix(10), 1, one_locus)
  expect_gt(far[1] - 0.1, 0.1 * 0.9 * 0.5 * 10 / 10)
  near <- loci_select(rare, matrix(TRUE), matrix(1), matrix(10), 1, one_locus)
  expect_lt(near[1] - 0.1, 0.1 * 0.9 * 0.5 * 1 / 10)
})

test_that("many loci of small effect respond as the infinitesimal model", {
  n_loci <- 200
  loci <- list(n = n_loci, effect = 0.2, mutation = 0, origin = 0)
  p <- array(0.5, c(1, 1, n_loci))
  g0 <- loci_variance(p, loci)[1, 1]
  expect_equal(g0, 4)
  theta <- 45
  z_inf <- z0 <- loci_mean(p, loci)[1, 1]
  expect_equal(z0, 40)
  for (i in 1:5) {
    z <- loci_mean(p, loci)
    width <- 25 + loci_variance(p, loci) + 1
    p <- loci_select(p, matrix(TRUE), theta - z, width, 1, loci)
    z_inf <- z_inf + g0 * (theta - z_inf) / (25 + g0 + 1)
  }
  expect_equal(loci_mean(p, loci)[1, 1] - z0, z_inf - z0, tolerance = 0.005)
  expect_gt(z_inf - z0, 2)
  # the recursion above is the engine with a constant genetic variance and no drift
  ev <- evolution_control(genetic_var = g0, ne_ratio = Inf, env_var = 1, founder = FALSE,
                          env = theta, traits = loci_traits(z0, sigma = 5))
  out <- mc_dynamics(matrix(5, 1, 1), matrix(0.001), dispersal = 0, n_steps = 5,
                     seed_steps = numeric(0), initial = one_sp(500L, 1), evolution = ev)
  expect_equal(out$evolution$genetic[1, 1], z_inf)
  expect_equal(loci_mean(p, loci)[1, 1], out$evolution$genetic[1, 1], tolerance = 0.001)
  # the variance has hardly changed, as the infinitesimal model assumes
  expect_equal(loci_variance(p, loci)[1, 1], g0, tolerance = 0.01)
})

test_that("the frequencies of the pool give the optimum of each species", {
  loci <- list(n = 10, effect = 5, mutation = 1e-4, origin = 0)
  p <- loci_initial(c(40, 60, 0, 100, 47), loci)
  expect_identical(dim(p), c(5L, 10L))
  expect_identical(p[1, ], rep(c(1, 0), c(4, 6)))
  expect_identical(p[2, ], rep(c(1, 0), c(6, 4)))
  expect_identical(p[3, ], rep(0, 10))
  expect_identical(p[4, ], rep(1, 10))
  expect_equal(p[5, ], c(1, 1, 1, 1, 0.7, 0, 0, 0, 0, 0))
  expect_equal(loci_mean(array(p, c(5, 1, 10)), loci)[, 1], c(40, 60, 0, 100, 47))
  expect_equal(loci_variance(array(p, c(5, 1, 10)), loci)[, 1],
               c(0, 0, 0, 0, 2 * 25 * 0.7 * 0.3))
  # an origin other than 0, and optima that are whole numbers of alleles but for rounding
  shifted <- list(n = 3, effect = 0.1, mutation = 0, origin = -0.3)
  expect_identical(loci_initial(c(-0.3, 0.1, 0.3), shifted),
                   rbind(c(0, 0, 0), c(1, 1, 0), c(1, 1, 1)))
  expect_error(loci_initial(c(40, 101), loci), "between `origin` \\(0\\) and .* \\(100\\)")
  expect_error(loci_initial(-1, loci), "optimum of every species")
  s <- loci_grid()
  s$traits$optimum[2] <- 120
  expect_error(simulate_metacommunity(s$ls, s$traits, s$alpha, n_steps = 2,
                                      evolution = evolution_control(loci = loci)),
               "optimum of every species")
  # populations at the start carry the frequencies of the pool, without sampling
  s <- loci_grid()
  set.seed(14)
  sim <- simulate_metacommunity(s$ls, s$traits, s$alpha, dispersal = 0, n_steps = 1,
                                seed_steps = integer(0), initial = s$initial,
                                evolution = evolution_control(
                                  loci = list(n = 10, effect = 5, mutation = 0),
                                  ne_ratio = Inf
                                ))
  here <- sim$abundance > 0
  expect_equal(sim$evolution$genetic[here], c(40, 60)[col(here)[here]])
  expect_true(all(sim$evolution$genetic_variance[here] == 0))
})

test_that("loci are checked and refuse the settings they cannot go with", {
  loci <- list(n = 10, effect = 5, mutation = 1e-4)
  ev <- evolution_control(loci = loci, ne_ratio = 0.5, env_var = 4)
  expect_s3_class(ev, "mc_evolution")
  expect_identical(ev$loci, list(n = 10, effect = 5, mutation = 1e-4, origin = 0))
  expect_identical(ev$genetic_var, 0)
  expect_output(print(ev), "G = from the loci .* loci: 10 of effect 5")
  expect_null(evolution_control(1)$loci)
  expect_false("loci" %in% names(evolution_control(1)))
  expect_warning(evolution_control(genetic_var = 3, loci = loci), "`genetic_var` is ignored")
  expect_error(evolution_control(loci = loci, mutation_var = 0.1), "`mutation_var`")
  expect_error(evolution_control(loci = loci, variance_selection = TRUE),
               "`variance_selection`")
  expect_error(evolution_control(loci = loci, variance_selection = "sexual"),
               "`variance_selection`")
  expect_error(evolution_control(loci = loci, immigrant_selection = TRUE),
               "`immigrant_selection = TRUE`")
  expect_error(evolution_control(loci = loci, env_var = diag(2)), "one trait")
  for (bad in list(10, list(10, 5, 1e-4), list(n = 10, effect = 5), "ten",
                   c(loci, list(dominance = 1)), c(loci, list(n = 2)))) {
    expect_error(evolution_control(loci = bad), "`loci` must be `NULL` or a list")
  }
  expect_error(evolution_control(loci = list(n = 1.5, effect = 5, mutation = 0)), "loci\\$n")
  expect_error(evolution_control(loci = list(n = 0, effect = 5, mutation = 0)), "loci\\$n")
  expect_error(evolution_control(loci = list(n = 2, effect = 0, mutation = 0)),
               "loci\\$effect")
  expect_error(evolution_control(loci = list(n = 2, effect = 1, mutation = 0.6)),
               "loci\\$mutation")
  expect_error(evolution_control(loci = list(n = 2, effect = 1, mutation = -0.1)),
               "loci\\$mutation")
  expect_error(evolution_control(loci = list(n = 2, effect = 1, mutation = 0, origin = NA)),
               "loci\\$origin")
  # several traits and the engine with speciation
  ls <- landscape_grid(3, disp_scale = 0.5)
  tr <- list(a = data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1"),
             b = data.frame(optimum = 0, sigma = Inf, r_max = 1, row.names = "sp1"))
  expect_error(simulate_metacommunity(ls, tr, matrix(0.01, 1, 1, dimnames = list("sp1", "sp1")),
                                      env = cbind(a = ls$env, b = 0), n_steps = 2,
                                      evolution = ev), "several traits")
  expect_error(simulate_evolving_metacommunity(ls, tr$a, speciation_control("geographic", Inf),
                                               evolution = ev, n_steps = 2),
               "`loci` are not yet available")
  s <- loci_grid()
  expect_error(simulate_metacommunity(
    s$ls, s$traits, s$alpha, n_steps = 2, recruitment = "lottery", microsites = 50,
    evolution = evolution_control(loci = loci, steps_per_generation = 2)
  ), "steps_per_generation = 1")
})

test_that("gene flow, founders, seeds and disturbance keep the loci consistent", {
  s <- loci_grid()
  loci <- list(n = 10, effect = 5, mutation = 0.01)
  for (founder in c(TRUE, FALSE)) {
    set.seed(15)
    sim <- simulate_metacommunity(s$ls, s$traits, s$alpha, dispersal = 0.05, n_steps = 150,
                                  seed_steps = c(20, 40), initial = s$initial,
                                  disturbance = 0.05, record_steps = c(75, 150),
                                  evolution = evolution_control(loci = loci, ne_ratio = 0.5,
                                                                env_var = 4,
                                                                founder = founder))
    expect_loci_state(sim, loci)
    # mutation and gene flow built a variance from a pool without any
    expect_gt(mean(sim$evolution$genetic_variance, na.rm = TRUE), 1)
    expect_identical(unname(sim$evolution$history[, , "step150"]),
                     unname(sim$evolution$genetic))
  }
  # a pool with a polymorphic locus, on an implicit landscape, with pruned species
  set.seed(16)
  ls <- landscape_implicit(10)
  traits <- data.frame(optimum = c(0.2, 0.55, 0.9), sigma = 0.3, r_max = c(5, 5, 0.01),
                       row.names = paste0("sp", 1:3))
  alpha <- matrix(0.02, 3, 3, dimnames = list(rownames(traits), rownames(traits)))
  wide <- list(n = 4, effect = 0.125, mutation = 1e-3)
  sim <- simulate_metacommunity(ls, traits, alpha, dispersal = 0.05, n_steps = 120,
                                seed_steps = c(1, 5), extirpation = 0.02,
                                evolution = evolution_control(loci = wide, ne_ratio = 1))
  expect_loci_state(sim, wide)
  expect_true(all(sim$abundance[, "sp3"] == 0))
  # the same seed gives the same run
  set.seed(16)
  again <- simulate_metacommunity(landscape_implicit(10), traits, alpha, dispersal = 0.05,
                                  n_steps = 120, seed_steps = c(1, 5), extirpation = 0.02,
                                  evolution = evolution_control(loci = wide, ne_ratio = 1))
  expect_identical(again$evolution$allele_frequency, sim$evolution$allele_frequency)
})

test_that("the lottery keeps the loci consistent, with and without overlapping generations", {
  s <- loci_grid()
  loci <- list(n = 10, effect = 5, mutation = 0.01)
  settings <- expand.grid(mortality = c(1, 0.3), founder = c(TRUE, FALSE), nu = c(0.5, Inf))
  for (i in seq_len(nrow(settings))) {
    set.seed(17)
    sim <- simulate_metacommunity(s$ls, s$traits, s$alpha, dispersal = 0.05, n_steps = 120,
                                  seed_steps = c(20, 40), initial = s$initial,
                                  disturbance = 0.05, recruitment = "lottery",
                                  microsites = 100, mortality = settings$mortality[i],
                                  evolution = evolution_control(
                                    loci = loci, ne_ratio = settings$nu[i], env_var = 4,
                                    founder = settings$founder[i]
                                  ))
    expect_loci_state(sim, loci)
    expect_gt(mean(sim$evolution$genetic_variance, na.rm = TRUE), 1)
  }
  # without mutation, drift or selection (a wide niche), a lottery leaves the pool as it is
  set.seed(18)
  s$traits$sigma <- 1e8
  sim <- simulate_metacommunity(s$ls, s$traits, s$alpha, dispersal = 0, n_steps = 30,
                                seed_steps = integer(0), initial = s$initial,
                                recruitment = "lottery", microsites = 100, mortality = 0.3,
                                evolution = evolution_control(
                                  loci = list(n = 10, effect = 5, mutation = 0),
                                  ne_ratio = 0.5
                                ))
  here <- sim$abundance > 0
  expect_identical(sim$evolution$genetic[here], c(40, 60)[col(here)[here]])
})

test_that("drift in a lottery loses heterozygosity at the rate of the recruits and survivors", {
  set.seed(19)
  m <- 1500
  steps <- 15
  ev <- evolution_control(ne_ratio = 0.5, env = rep(0, m), traits = loci_traits(1),
                          loci = list(n = 1, effect = 1, mutation = 0))
  # discrete generations: every microsite is refilled, by a sample of 2 * 0.5 * 60 alleles
  out <- mc_dynamics(matrix(5, m, 1), matrix(0), dispersal = 0, n_steps = steps,
                     seed_steps = numeric(0), initial = one_sp(60L, m), evolution = ev,
                     recruitment = "lottery", microsites = 60, mortality = 1)
  p <- out$evolution$allele_frequency[, 1, 1]
  expect_equal(mean(2 * p * (1 - p)), 0.5 * (1 - 1 / 60)^steps, tolerance = 0.05)
  # overlapping generations lose it more slowly per step
  slow <- mc_dynamics(matrix(5, m, 1), matrix(0), dispersal = 0, n_steps = steps,
                      seed_steps = numeric(0), initial = one_sp(60L, m), evolution = ev,
                      recruitment = "lottery", microsites = 60, mortality = 0.2)
  q <- slow$evolution$allele_frequency[, 1, 1]
  expect_gt(mean(2 * q * (1 - q)), mean(2 * p * (1 - p)) + 0.02)
  expect_lt(mean(2 * q * (1 - q)), 0.49)
})

test_that("individuals are sampled with the genetic variance of the loci of each population", {
  s <- loci_grid()
  loci <- list(n = 10, effect = 5, mutation = 0.01)
  set.seed(22)
  sim <- simulate_metacommunity(s$ls, s$traits, s$alpha, dispersal = 0.05, n_steps = 150,
                                seed_steps = integer(0), initial = s$initial,
                                record_steps = c(75, 150),
                                evolution = evolution_control(loci = loci, ne_ratio = 0.5,
                                                              env_var = 4))
  g <- sim$evolution$genetic_variance
  sites <- rownames(sim$abundance)
  state <- evolved_populations(sim, sites, c("sp1", "sp2"))
  here <- sim$abundance > 0
  expect_equal(state$var[here], g[here] + 4)
  expect_true(all(state$var[!here] == 4))
  # with `mutation_var`, the variance of each population also takes the place of `genetic_var`
  set.seed(23)
  mut <- simulate_metacommunity(s$ls, s$traits, s$alpha, dispersal = 0.05, n_steps = 60,
                                seed_steps = integer(0), initial = s$initial,
                                record_steps = c(30, 60),
                                evolution = evolution_control(50, ne_ratio = 0.5, env_var = 4,
                                                              mutation_var = 0.5))
  there <- mut$abundance > 0
  now <- evolved_populations(mut, sites, c("sp1", "sp2"))$var
  expect_equal(now[there], mut$evolution$genetic_variance[there] + 4)
  expect_false(isTRUE(all.equal(now[there], rep(54, sum(there)))))
  # before the last step it is not kept: the variance of the species is used
  expect_true(all(evolved_populations(mut, sites, c("sp1", "sp2"), step = 30)$var == 54))
  # the last step, also when it was recorded, but no step before it
  expect_identical(evolved_populations(sim, sites, c("sp1", "sp2"), step = 150)$var, state$var)
  expect_error(evolved_populations(sim, sites, c("sp1", "sp2"), step = 75), "last step only")
  # the variance among the sampled individuals of a large population follows it
  big <- which(sim$abundance == max(sim$abundance), arr.ind = TRUE)[1, ]
  ind <- sample_individuals(sim, n = 5000, source = "evolution", intraspecific_sd = 0)
  mine <- ind[ind$site == sites[big[1]] & ind$species == colnames(g)[big[2]], ]
  expect_gt(nrow(mine), 50)
  expect_gt(g[big[1], big[2]], 20)
  expect_equal(var(mine$trait), g[big[1], big[2]] + 4, tolerance = 0.5)
})
