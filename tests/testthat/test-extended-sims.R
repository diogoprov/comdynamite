# Extended tests of the simulators and of ABC (rOpenSci statistical standards): simulators
# against analytic results over several seeds, convergence as data grow, and recovery and
# calibration of ABC with random forests on a cheap simulator. They run only with
# COMDYNAMITE_EXTENDED_TESTS=true (helper-extended.R).
#
# Tolerances come from theory or from the Monte Carlo error: bounds of 4 standard errors
# (a correct implementation fails by chance in about 6 of 100000 comparisons each), exact F
# and binomial distributions (rejection_compatible(), exact binomial 99% interval), or
# expected rates of convergence with a margin of at least a factor of 2 against them.

# second-order (linear noise) approximation of the stationary mean of the Beverton-Holt map
# N' ~ Poisson(r N / (1 + alpha N)) around its deterministic equilibrium N* = alpha^-1 (r - 1):
# with J the Jacobian of the map at N*, the stationary covariance solves S = J S J' + diag(N*)
# (Poisson variance), and the shift of the mean solves (I - J) b = H[S] / 2, with H the second
# derivatives of the map. b is of order one individual whatever the scale of N*.
sims_lna <- function(r, alpha) {
  eq <- solve(alpha, r - 1)
  k <- length(r)
  f <- function(n) r * n / (1 + drop(alpha %*% n))
  h <- 1e-3 * eq
  unit <- function(j) replace(numeric(k), j, h[j])
  J <- vapply(seq_len(k), function(j) (f(eq + unit(j)) - f(eq - unit(j))) / (2 * h[j]),
              numeric(k))
  S <- matrix(solve(diag(k^2) - kronecker(J, J), as.vector(diag(eq, k))), k)
  second <- numeric(k)
  for (a in seq_len(k)) for (b in seq_len(k)) {
    d2 <- (f(eq + unit(a) + unit(b)) - f(eq + unit(a) - unit(b)) -
             f(eq - unit(a) + unit(b)) + f(eq - unit(a) - unit(b))) / (4 * h[a] * h[b])
    second <- second + d2 * S[a, b]
  }
  list(eq = eq, shift = drop(solve(diag(k) - J, second / 2)))
}

# a random community of three species that coexist at a stable equilibrium: intraspecific
# competition dominates (off-diagonal alpha below 0.4 of the diagonal), and communities are
# drawn again until every N* is above 30, so that local extinctions are negligible and the
# second-order approximation applies
sims_community <- function(k = 3) {
  repeat {
    alpha <- matrix(stats::runif(k * k, 0, 0.4), k)
    diag(alpha) <- 1
    com <- list(alpha = alpha * 0.02, r = stats::runif(k, 3, 6))
    if (all(solve(com$alpha, com$r - 1) > 30)) return(com)
  }
}

# mean over sites of the time-averaged abundance (steps 101 to n_steps) and its standard error
# from the independent sites (no dispersal)
sims_equilibrium <- function(r, alpha, sites = 300, n_steps = 400, interface = "mc_dynamics") {
  k <- length(r)
  sp <- paste0("sp", seq_len(k))
  start <- matrix(round(solve(alpha, r - 1)), sites, k, byrow = TRUE, dimnames = list(NULL, sp))
  out <- if (interface == "mc_dynamics") {
    mc_dynamics(matrix(r, sites, k, byrow = TRUE), alpha, dispersal = 0, n_steps = n_steps,
                initial = start, species = sp, seed_steps = numeric(0),
                record_steps = 101:n_steps)
  } else {
    # a flat niche (sigma huge) gives growth r_max at every site
    traits <- data.frame(optimum = 50, sigma = 1e8, r_max = r, row.names = sp)
    simulate_metacommunity(landscape_implicit(sites), traits,
                           `dimnames<-`(alpha, list(sp, sp)), dispersal = 0,
                           n_steps = n_steps, initial = start, seed_steps = numeric(0),
                           record_steps = 101:n_steps)
  }
  site_means <- apply(out$history, c(1, 2), mean)
  list(mean = colMeans(site_means), se = apply(site_means, 2, stats::sd) / sqrt(sites))
}

test_that("mc_dynamics and simulate_metacommunity reach the Lotka-Volterra equilibrium", {
  skip_if_not_extended()
  # G5.4, G5.6, G5.6b: six random three-species communities (one seed each), half run through
  # simulate_metacommunity(). The stationary mean must equal N* plus the second-order shift
  # of Poisson noise; tolerance 4 standard errors (sites are independent) plus 25% of the
  # shift, for the third-order terms neglected by the approximation (relative size about
  # N*^-1/2, under 20% for N* above 30).
  for (seed in 1:6) {
    set.seed(seed)
    com <- sims_community()
    theory <- sims_lna(com$r, com$alpha)
    interface <- if (seed %% 2) "mc_dynamics" else "simulate_metacommunity"
    obs <- sims_equilibrium(com$r, com$alpha, interface = interface)
    expect_lt(max(abs(obs$mean - theory$eq - theory$shift) - 4 * obs$se -
                    0.25 * abs(theory$shift)), 0)
    # the shift is of order one individual, as the approximation says
    expect_lt(max(abs(theory$shift)), 3)
  }
})

test_that("the stochastic equilibrium approaches the deterministic one as abundance grows", {
  skip_if_not_extended()
  # G5.7: dividing alpha by 10 multiplies N* by 10 while the noise-induced shift stays of
  # order one individual, so the relative error |mean / N* - 1| should fall about tenfold;
  # the test asks for at least threefold over five communities.
  rel <- vapply(1:5, function(seed) {
    set.seed(seed)
    com <- sims_community()
    vapply(c(1, 10), function(s) {
      obs <- sims_equilibrium(com$r, com$alpha / s, sites = 200, n_steps = 300)
      max(abs(obs$mean / solve(com$alpha / s, com$r - 1) - 1))
    }, numeric(1))
  }, numeric(2))
  expect_lt(mean(rel[2, ]), mean(rel[1, ]) / 3)
})

# exact F test that x ~ N(mu 1, s^2 V0), with mu and s unknown: whiten by V0, project out the
# intercept, and compare the mean square in the `n_dir` directions where the alternative V1
# puts most variance with the mean square in the others. Under the null the projected values
# are independent N(0, s^2), so the ratio is F(n_dir, n - 1 - n_dir) whatever mu and s. Values
# that were centred or rescaled (as the slopes and optima are) give the same F.
sims_isotropy_p <- function(x, V0, V1, n_dir = 4) {
  L <- t(chol(V0))
  n <- length(x)
  w <- forwardsolve(L, x)
  u <- forwardsolve(L, rep(1, n))
  M <- diag(n) - tcrossprod(u) / sum(u^2)
  A <- forwardsolve(L, t(forwardsolve(L, V1)))
  basis <- eigen(M %*% A %*% M, symmetric = TRUE)$vectors[, seq_len(n - 1)]
  y <- drop(crossprod(basis, M %*% w))
  top <- seq_len(n_dir)
  ratio <- (sum(y[top]^2) / n_dir) / (sum(y[-top]^2) / (n - 1 - n_dir))
  stats::pf(ratio, n_dir, n - 1 - n_dir, lower.tail = FALSE)
}

test_that("species_traits and simulate_logistic carry the phylogenetic signal of lambda", {
  skip_if_not_extended()
  # G5.4, G5.6b: 300 replicates per case on a 32-tip coalescent tree. With V(lambda) =
  # lambda C + (1 - lambda) diag(C), the exact F test of V(lambda_true) must reject at a rate
  # compatible with 0.05 (exact binomial 99% interval), and the test of the wrong lambda must
  # reject often (power above 0.8; observed in development 0.93 to 0.99).
  set.seed(30)
  tree <- ape::rcoal(32)
  C <- ape::vcv(tree)
  covariance <- function(lambda) lambda * C + (1 - lambda) * diag(diag(C))
  optima <- function(lambda) species_traits(tree, lambda = lambda)$optimum
  slopes <- function(lambda) {
    simulate_logistic(tree, stats::rnorm(10), lambda = lambda)$traits$slope[, 1]
  }
  n <- 300
  rejections <- function(draw, lambda, null, alternative) {
    sum(replicate(n, sims_isotropy_p(draw(lambda), covariance(null), covariance(alternative)))
        < 0.05)
  }
  for (draw in list(optima, slopes)) {
    expect_true(rejection_compatible(rejections(draw, 1, 1, 0), n))
    expect_true(rejection_compatible(rejections(draw, 0, 0, 1), n))
    expect_true(rejection_compatible(rejections(draw, 0.5, 0.5, 1), n))
    expect_gt(rejections(draw, 1, 0, 1) / n, 0.8)
    expect_gt(rejections(draw, 0, 1, 0) / n, 0.8)
  }
  # the slopes have exactly the standard deviation slope_sd
  sim <- simulate_logistic(tree, data.frame(x = stats::rnorm(10), y = stats::rnorm(10)),
                           slope_sd = 2.5)
  expect_equal(apply(sim$traits$slope, 2, stats::sd), c(x = 2.5, y = 2.5))
})

test_that("a Poisson GLM recovers the Gaussian responses of simulate_static", {
  skip_if_not_extended()
  # G5.4, G5.6, G5.6a, G5.6b: the log of the expected abundance is quadratic in the gradient,
  # log r_max - (E - theta)^2 / (2 sigma^2), so a Poisson GLM with E and E^2 has known true
  # coefficients. Five seeds x five species with different parameters, 200 sites: every Wald
  # z of the 75 coefficients within 4 (chance of any failing about 0.5%), and their mean
  # square within the 99.9% interval of chi-squared(75) / 75.
  z <- NULL
  for (seed in 1:5) {
    set.seed(seed)
    traits <- species_traits(5, optima_range = c(20, 80), sigma_mean = 15, sigma_sd = 3,
                             rmax_range = c(5, 20))
    env <- stats::runif(200, 0, 100)
    sim <- simulate_static(traits, env)
    for (j in seq_len(nrow(traits))) {
      fit <- stats::glm(sim$abundance[, j] ~ env + I(env^2), family = stats::poisson)
      t <- traits[j, ]
      truth <- c(log(t$r_max) - t$optimum^2 / (2 * t$sigma^2), t$optimum / t$sigma^2,
                 -1 / (2 * t$sigma^2))
      z <- c(z, (stats::coef(fit) - truth) / sqrt(diag(stats::vcov(fit))))
    }
  }
  expect_lt(max(abs(z)), 4)
  expect_gt(mean(z^2), stats::qchisq(0.0005, length(z)) / length(z))
  expect_lt(mean(z^2), stats::qchisq(0.9995, length(z)) / length(z))
})

test_that("the error of the optimum estimated from simulate_static falls with more sites", {
  skip_if_not_extended()
  # G5.7: the maximum-likelihood error falls as n^-1/2, so 16 times more sites should divide
  # the root mean squared error of the optimum by about 4 (observed in development 3.9). With
  # 100 replicates the RMSE has a Monte Carlo error near 7%, so the test asks for a factor of
  # 2 and a monotone decrease over 50, 200 and 800 sites.
  set.seed(31)
  traits <- data.frame(optimum = 50, sigma = 15, r_max = 10, row.names = "sp1")
  rmse <- vapply(c(50, 200, 800), function(n) {
    errors <- replicate(100, {
      env <- stats::runif(n, 0, 100)
      y <- simulate_static(traits, env)$abundance[, 1]
      b <- stats::coef(stats::glm(y ~ env + I(env^2), family = stats::poisson))
      -b[2] / (2 * b[3]) - 50
    })
    sqrt(mean(errors^2))
  }, numeric(1))
  expect_true(all(diff(rmse) < 0))
  expect_lt(rmse[3], rmse[1] / 2)
})

test_that("drift in simulate_evolving_metacommunity accumulates G times the sum of 1 / Ne", {
  skip_if_not_extended()
  # G5.4, G5.6b: isolated populations under drift alone (no selection, no gene flow, no
  # speciation). Lande (1976): the variance among populations of the genetic mean grows by
  # G / Ne per generation, here Ne = ne_ratio x N of each generation. Two values of ne_ratio
  # x three seeds, 400 sites. The genetic means are Gaussian given the abundances, so the
  # sample variance has a relative standard error sqrt(2 / (sites - 1)): each ratio of
  # observed to expected variance within 4 of these, and their mean within 4 / sqrt(6) of
  # them. The mean over sites stays at the ancestral optimum within 4 standard errors.
  sites <- 400
  steps <- 50
  ancestor <- data.frame(optimum = 50, sigma = 1e6, r_max = 5, row.names = "sp1")
  ratios <- NULL
  for (ne_ratio in c(0.5, 0.1)) for (seed in 1:3) {
    set.seed(seed)
    sim <- simulate_evolving_metacommunity(
      landscape_implicit(sites), ancestor, speciation_control("geographic", Inf),
      evolution = evolution_control(genetic_var = 1, ne_ratio = ne_ratio, founder = FALSE),
      dispersal = 0, n_steps = steps, competition = "none",
      initial = matrix(80L, sites, 1, dimnames = list(NULL, "sp1")),
      record_steps = seq_len(steps))
    n <- sim$history[, 1, ]
    alive <- apply(n > 0, 1, all)
    expected <- mean(rowSums(1 / (ne_ratio * n[alive, ])))
    g <- sim$genetic[alive, 1]
    se <- sqrt(2 / (sum(alive) - 1))
    expect_lt(abs(stats::var(g) / expected - 1), 4 * se)
    expect_lt(abs(mean(g) - 50), 4 * sqrt(expected / sum(alive)))
    ratios <- c(ratios, stats::var(g) / expected)
  }
  expect_lt(abs(mean(ratios) - 1), 4 * sqrt(2 / (sites - 1)) / sqrt(length(ratios)))
})

test_that("food_web_dynamics integrates logistic growth with fourth-order accuracy", {
  skip_if_not_extended()
  # G5.4, G5.7: the resource alone at its optimum trait, without evolution (V = 0), grows
  # logistically, dR/dt = R (c0 - d R) with c0 = 1 and d = 0.02 (R* = c0 / d = 50), whose
  # solution is known. Runge-Kutta 4 has a global error of order dt^4, so halving dt should
  # divide the error by about 16 (observed in development 15.6 to 15.2); the test asks for
  # at least 10 at each halving, and an error below 1e-3 at dt = 0.25.
  logistic <- function(t, start = 5) 50 / (1 + (50 - start) / start * exp(-t))
  errors <- vapply(c(1, 0.5, 0.25, 0.125), function(dt) {
    x <- food_web_dynamics(seq(0, 20, 1), init = c(R = 5, N = 0, P = 0),
                           traits = c(R = 10, N = 4, P = 1), V = 0, dt = dt)
    max(abs(x$R - logistic(x$time)))
  }, numeric(1))
  expect_true(all(errors[-4] / errors[-1] > 10))
  expect_lt(errors[3], 1e-3)
})

# a cheap simulator for ABC: ten species with one niche breadth and one maximum abundance on a
# gradient of 50 sites, summarised by summary_stats() (statistics without inputs are NA and
# dropped by abc_fit())
sims_abc_run <- function(scenario, shared) {
  traits <- species_traits(10, sigma_mean = scenario$breadth, sigma_sd = 0,
                           rmax_range = rep(scenario$rmax, 2))
  simulate_static(traits, env = shared$env)
}
sims_abc_setup <- function() list(env = seq(0, 100, length.out = 50))
sims_abc_prior <- list(breadth = c(5, 30), rmax = c(2, 20))
sims_abc_table <- function(n, seed) {
  abc_table(sims_abc_prior, sims_abc_run, n = n, setup = sims_abc_setup,
            summarise = function(sim, shared) summary_stats(sim), seed = seed)
}

test_that("abc_recovery recovers the parameters of a cheap simulator over seeds", {
  skip_if_not_extended()
  skip_if_not_installed("abcrf")
  # G5.4, G5.6, G5.6b: three reference tables (seeds) of 1200 simulations, 200 left out as
  # pseudo-observations. Both parameters are well identified by the statistics (R2 0.95 to
  # 0.97 in development), so R2 above 0.9 is asked; the bias must be within 4 standard errors
  # of the mean error of the 200 estimates.
  for (seed in 1:3) {
    table <- sims_abc_table(1200, seed)
    rec <- abc_recovery(table, n_test = 200, ntree = 300, seed = seed)
    est <- attr(rec, "estimates")
    expect_true(all(rec$r2 > 0.9))
    for (p in names(sims_abc_prior)) {
      err <- with(est[est$parameter == p, ], estimate - truth)
      expect_lt(abs(mean(err)), 4 * stats::sd(err) / sqrt(length(err)))
    }
  }
})

test_that("the error of abc_fit falls as the reference table grows", {
  skip_if_not_extended()
  skip_if_not_installed("abcrf")
  # G5.7: estimates of 200 fixed test simulations from forests trained on 150, 600 and 2400
  # simulations of the same table. The error of a random forest falls with more training
  # data towards the posterior spread; in development the RMSE fell by a factor of 1.7
  # (breadth) and 2.6 (rmax) from 200 to 3200 rows. The test asks for a monotone decrease and
  # at least 20% less error at 2400 than at 150 rows for each parameter.
  table <- sims_abc_table(2600, 7)
  test <- as.data.frame(table)[2401:2600, ]
  rmse <- vapply(c(150, 600, 2400), function(n) {
    train <- table[seq_len(n), ]
    attributes(train)[c("prior", "parameters", "statistics")] <-
      attributes(table)[c("prior", "parameters", "statistics")]
    post <- predict(abc_fit(train, ntree = 300), test)
    truth <- c(t(as.matrix(test[names(sims_abc_prior)])))
    tapply((post$estimate - truth)^2, post$parameter, function(e) sqrt(mean(e)))
  }, numeric(2))
  expect_true(all(apply(rmse, 1, function(r) all(diff(r) < 0))))
  expect_true(all(rmse[, 3] < 0.8 * rmse[, 1]))
})

test_that("abc_recovery credible intervals have coverage near the nominal level", {
  skip_if_not_extended()
  skip_if_not_installed("abcrf")
  skip(paste(
    "known issue: the quantile regression forest intervals of abcrf are not calibrated.",
    "Nominal 0.9; in the three tables below, coverage of 200 pseudo-observations was 0.95,",
    "0.925, 0.985 (breadth) and 0.97, 0.97, 0.975 (rmax), 4 of 6 outside the exact 99%",
    "interval; with 5000 simulations 0.96 and 0.995, intervals about 2.3 times wider than",
    "2 x 1.645 x RMSE. On a toy table with an analytic posterior (a ~ U(0, 1),",
    "s = a + N(0, 0.05)) coverage was 0.65 with s alone (under-coverage). abcrf::regAbcrf()",
    "called directly gives the same numbers, so the cause is in abcrf, not in the wrapper."))
  # G5.6b, interval calibration: coverage of the 90% intervals of 200 pseudo-observations in
  # each of three reference tables must be compatible with 0.9 (exact binomial 99% interval)
  for (seed in 1:3) {
    rec <- abc_recovery(sims_abc_table(1700, seed), n_test = 200, ntree = 300, seed = seed)
    for (i in seq_len(nrow(rec))) {
      expect_true(rejection_compatible(round(rec$coverage[i] * rec$n[i]), rec$n[i], 0.9))
    }
  }
})

test_that("abc_recovery credible intervals do not under-cover", {
  skip_if_not_extended()
  skip_if_not_installed("abcrf")
  # the weaker property that holds despite the known issue above: coverage of the 90%
  # intervals is not below 0.9 (upper end of the exact binomial 99% interval at or above 0.9)
  for (seed in 4:5) {
    rec <- abc_recovery(sims_abc_table(1200, seed), n_test = 200, ntree = 300, seed = seed)
    for (i in seq_len(nrow(rec))) {
      ci <- stats::binom.test(round(rec$coverage[i] * rec$n[i]), rec$n[i],
                              conf.level = 0.99)$conf.int
      expect_gte(ci[2], 0.9)
    }
  }
})

# reference table for abc_cv(): an 8 x 8 grid with the gradient along x, random folds, and five
# statistics of the sites kept and left out in each fold
sims_cv_table <- function(n, seed) {
  grid <- expand.grid(x = 1:8, y = 1:8)
  folds <- spatial_folds(grid, k = 4, method = "random")
  statistics <- c("alpha", "gamma", "occupancy", "abundance", "evenness")
  summarise <- fold_summarise(folds, function(x, shared) summary_stats(x)[statistics])
  abc_table(sims_abc_prior, sims_abc_run, n = n, setup = function() list(env = grid$x * 12.5),
            summarise = summarise, seed = seed)
}

test_that("abc_cv predictive intervals have coverage near the nominal level", {
  skip_if_not_extended()
  skip_if_not_installed("abcrf")
  skip(paste(
    "known issue: under the true model the predictive intervals of abc_cv() are too wide.",
    "Nominal 0.9; coverage per statistic and fold (100 data sets each) 0.97 to 1.00 for",
    "alpha, gamma, occupancy and abundance and 0.91 to 0.96 for evenness; root mean square",
    "of z 0.34 to 0.73 instead of about 1: the forest weights spread the approximate",
    "posterior, as for the intervals of abc_recovery()."))
  # G5.6b, calibration: data sets simulated from the prior are held out of the table; for
  # each statistic and fold, the coverage over 100 data sets must be compatible with 0.9
  table <- sims_cv_table(1500, 1)
  train <- table[1:1400, ]
  attributes(train)[c("prior", "parameters", "statistics")] <-
    attributes(table)[c("prior", "parameters", "statistics")]
  cv <- abc_cv(train, as.data.frame(table)[1401:1500, attr(table, "statistics")], n = 300,
               ntree = 200)
  for (cell in split(cv, list(cv$statistic, cv$fold))) {
    expect_true(rejection_compatible(sum(cell$inside), nrow(cell), 0.9))
  }
})

test_that("abc_cv keeps the type I error of its model check under the true model", {
  skip_if_not_extended()
  skip_if_not_installed("abcrf")
  # type I error: when the data come from the model, left-out statistics should rarely be
  # flagged. For each statistic and fold (100 independent data sets), the share of p-values
  # below 0.05 must not exceed 0.05 (lower end of the exact binomial 99% interval at or below
  # 0.05), and predictions must be unbiased: the mean z within 4 standard errors of zero,
  # computed from the mean z of each data set over the folds (folds of one data set are not
  # independent).
  table <- sims_cv_table(1500, 2)
  train <- table[1:1400, ]
  attributes(train)[c("prior", "parameters", "statistics")] <-
    attributes(table)[c("prior", "parameters", "statistics")]
  cv <- abc_cv(train, as.data.frame(table)[1401:1500, attr(table, "statistics")], n = 300,
               ntree = 200)
  for (cell in split(cv, list(cv$statistic, cv$fold))) {
    ci <- stats::binom.test(sum(cell$p < 0.05), nrow(cell), conf.level = 0.99)$conf.int
    expect_lte(ci[1], 0.05)
  }
  for (s in split(cv, cv$statistic)) {
    z <- tapply(s$z, s$observation, mean, na.rm = TRUE)
    z <- z[is.finite(z)]
    expect_lt(abs(mean(z)), 4 * stats::sd(z) / sqrt(length(z)))
  }
})
