# rOpenSci statistical standards tested in this file (see the srr tags in R/)
#' @srrstats {G5.6, G5.7} Tested in this file, as described in R/srr-stats-general.R.
NULL

# Extended tests of the evolutionary and DAMOCLES functions (rOpenSci statistical standards):
# parameter recovery over several seeds, convergence as data grow, and the type I error and
# power of the tests. They run only with COMDYNAMITE_EXTENDED_TESTS=true (helper-extended.R).
#
# Recovery is checked against theory wherever possible rather than with tuned tolerances:
# - likelihood-ratio coverage: when the model is the one that generated the data, twice the
#   gap between the maximum log-likelihood and the log-likelihood at the true parameters is
#   asymptotically chi-squared with as many degrees of freedom as parameters, so the true
#   values fall outside the 95% likelihood region in about 5% of replicates
#   (rejection_compatible() with an exact binomial 99% interval);
# - exact sampling distributions (Poisson counts of changes, chi-squared sums of squared
#   contrasts, binomial counts of transitions), with bounds of 4 standard errors or 99.9%
#   quantiles, so that a correct implementation fails by chance less than once in 1000 runs.

# an ultrametric tree of height one, from treeshapesim when installed
evol_tree <- function(n) {
  if (requireNamespace("treeshapesim", quietly = TRUE)) {
    tree <- treeshapesim::sim_tree_shape(n)
  } else {
    tree <- ape::rcoal(n)
  }
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  tree
}

# Gaussian log-likelihood of x with mean zero and covariance V
mvn_loglik <- function(x, V) {
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) return(-Inf)
  z <- backsolve(R, x, transpose = TRUE)
  -0.5 * (length(x) * log(2 * pi) + 2 * sum(log(diag(R))) + sum(z^2))
}

# maximum of a log-likelihood of log-parameters from several starts
max_loglik <- function(f, starts) {
  best <- -Inf
  for (s in starts) {
    o <- tryCatch(stats::optim(s, function(p) -f(p), method = "Nelder-Mead",
                               control = list(reltol = 1e-10, maxit = 2000)),
                  error = function(e) NULL)
    if (!is.null(o) && -o$value > best) best <- -o$value
  }
  best
}

# --- DAMOCLES --------------------------------------------------------------------------------

test_that("damocles_fit() recovers gamma and mu: likelihood coverage and convergence", {
  skip_if_not_extended()
  # G5.4, G5.6, G5.6b: 100 replicates (each its own tree and seed) of 20 communities on 60
  # tips, gamma = mu = 1 shared; the true rates must lie inside the 95% likelihood region at
  # the nominal rate. mu = 0 is on the boundary, which only makes the region conservative.
  set.seed(101)
  n_rep <- 100
  outside <- vapply(seq_len(n_rep), function(i) {
    tree <- evol_tree(60)
    comm <- simulate_damocles(tree, 1, 1, n = 20)
    fit <- damocles_fit(comm, tree, shared = TRUE)
    lr <- 2 * (attr(fit, "loglik_total") - sum(damocles_loglik(comm, tree, 1, 1)))
    lr > stats::qchisq(0.95, 2)
  }, logical(1))
  expect_true(rejection_compatible(sum(outside), n_rep))
  # G5.7: the error of gamma shrinks with more communities; 40 replicates with 2 and with 40
  # communities. Asymptotically the error falls as 1 / sqrt(20), but with few communities mu
  # often sits on its boundary at 0 and the fall is slower (a pilot of 3 x 30 replicates gave
  # ratios of medians of 2.4 to 3.8), so we ask for a factor of 1.75 only.
  set.seed(102)
  log_error <- function(k) {
    vapply(seq_len(40), function(i) {
      tree <- evol_tree(60)
      fit <- damocles_fit(simulate_damocles(tree, 1, 1, n = k), tree, shared = TRUE)
      abs(log(fit$gamma[1]))
    }, numeric(1))
  }
  e2 <- stats::median(log_error(2))
  e40 <- stats::median(log_error(40))
  expect_lt(e40, e2 / 1.75)
  expect_lt(e40, 0.2)
})

test_that("damocles_lrt() holds its size and has power against rates that differ by clade", {
  skip_if_not_extended()
  # Type I error: 100 replicates of 5 communities on 60 tips with gamma = mu = 1 everywhere,
  # constant rates against gamma varying between the crown clade of one daughter of the root
  # and the rest (1 df). data-raw/validation_damocles.R found 0.06 with 100 replicates.
  lrt_p <- function(effect) {
    tree <- evol_tree(60)
    clade <- ape::extract.clade(tree, ape::Ntip(tree) + 2)$tip.label
    reg <- clade_regimes(tree, list(A = clade))
    gamma <- if (effect) ifelse(reg == "A", 3, 0.5) else 1
    comm <- simulate_damocles(tree, gamma, 1, n = 5)
    alternative <- damocles_fit(comm, tree, shared = TRUE, regimes = reg, vary = "gamma")
    damocles_lrt(damocles_fit(comm, tree, shared = TRUE), alternative)$p_value[2]
  }
  set.seed(103)
  p_null <- vapply(seq_len(100), function(i) lrt_p(FALSE), numeric(1))
  expect_true(rejection_compatible(sum(p_null < 0.05), 100))
  # Power: gamma 3 in the clade and 0.5 elsewhere, 30 replicates. The validation script found
  # 0.93 with 100 replicates; with that power, fewer than 21 of 30 rejections (0.7) has
  # binomial probability below 0.001.
  set.seed(104)
  p_alt <- vapply(seq_len(30), function(i) lrt_p(TRUE), numeric(1))
  expect_gte(mean(p_alt < 0.05), 0.7)
})

test_that("damocles_test() holds its size under the DAMOCLES model and detects clustering", {
  skip_if_not_extended()
  # Type I error with the true rates: 4 trees of 60 tips, 50 communities each simulated with
  # gamma = mu = 1, 99 null communities; the lower-tail p-values of MPD and MNTD against the
  # DAMOCLES null are then exact ranks, so P(p <= 0.05) = 0.05 (ties only make it smaller).
  set.seed(105)
  res <- do.call(rbind, lapply(1:4, function(i) {
    tree <- evol_tree(60)
    damocles_test(simulate_damocles(tree, 1, 1, n = 50), tree, n_sim = 99, gamma = 1, mu = 1)
  }))
  ok <- !is.na(res$p_mpd_damocles)
  expect_gt(sum(ok), 180)
  expect_true(rejection_compatible(sum(res$p_mpd_damocles[ok] <= 0.05), sum(ok)))
  expect_true(rejection_compatible(sum(res$p_mntd_damocles[ok] <= 0.05), sum(ok)))
  # the same with rates estimated, shared by the 50 communities of each tree (a parametric
  # bootstrap; with 50 communities the estimates are close enough for a nominal size)
  set.seed(106)
  est <- do.call(rbind, lapply(1:4, function(i) {
    tree <- evol_tree(60)
    damocles_test(simulate_damocles(tree, 1, 1, n = 50), tree, n_sim = 99, shared = TRUE)
  }))
  ok <- !is.na(est$p_mpd_damocles)
  expect_true(rejection_compatible(sum(est$p_mpd_damocles[ok] <= 0.05), sum(ok)))
  # Power: communities made of one clade of 8 to 15 tips (strong clustering), 20 trees.
  set.seed(107)
  p_clust <- vapply(seq_len(20), function(i) {
    repeat {
      tree <- evol_tree(60)
      sizes <- vapply(ape::Ntip(tree) + seq_len(tree$Nnode), function(nd) {
        length(ape::extract.clade(tree, nd)$tip.label)
      }, numeric(1))
      nodes <- ape::Ntip(tree) + which(sizes >= 8 & sizes <= 15)
      if (length(nodes)) break
    }
    tips <- ape::extract.clade(tree, nodes[1])$tip.label
    comm <- stats::setNames(as.numeric(tree$tip.label %in% tips), tree$tip.label)
    damocles_test(comm, tree, n_sim = 99, gamma = 1, mu = 1)$p_mpd_damocles
  }, numeric(1))
  expect_gte(mean(p_clust <= 0.05), 0.9)
})

# --- continuous traits -----------------------------------------------------------------------

test_that("trait_bm() has unit rate and recovers Pagel's lambda", {
  skip_if_not_extended()
  # G5.6b: with lambda = 1 the phylogenetic independent contrasts of 5 trees of 200 tips are
  # iid N(0, 1), so their sum of squares is chi-squared with 5 * 199 df: the rate estimate
  # (mean square) must lie between its 0.05% and 99.95% quantiles.
  set.seed(201)
  pics <- unlist(lapply(1:5, function(i) {
    tree <- evol_tree(200)
    ape::pic(trait_bm(tree), tree)
  }))
  df <- length(pics)
  expect_gte(sum(pics^2), stats::qchisq(0.0005, df))
  expect_lte(sum(pics^2), stats::qchisq(0.9995, df))
  # G5.4, G5.6: lambda = 0.6, 100 replicates of 100 tips; the rate is profiled out and the
  # true lambda must fall outside the 95% profile-likelihood interval at the nominal rate.
  set.seed(202)
  outside <- vapply(seq_len(100), function(i) {
    tree <- evol_tree(100)
    C <- ape::vcv.phylo(tree)
    x <- trait_bm(tree, lambda = 0.6)[rownames(C)]
    profile <- function(lambda) {
      V <- lambda * C + (1 - lambda) * diag(diag(C))
      s2 <- sum(backsolve(chol(V), x, transpose = TRUE)^2) / length(x)
      mvn_loglik(x, s2 * V)
    }
    best <- stats::optimize(profile, c(0, 1), maximum = TRUE)$objective
    best <- max(best, profile(0), profile(1))
    2 * (best - profile(0.6)) > stats::qchisq(0.95, 1)
  }, logical(1))
  expect_true(rejection_compatible(sum(outside), 100))
})

test_that("trait_ou() recovers alpha and sigma2, and the error shrinks with more tips", {
  skip_if_not_extended()
  # OU with root = theta = 0 has mean zero and covariance
  # sigma2 / (2 alpha) exp(-alpha (d_i + d_j - 2 s_ij)) (1 - exp(-2 alpha s_ij)), with d the
  # depth of the tips and s_ij their shared path.
  ou_loglik <- function(x, C) {
    d <- diag(C)
    gap <- outer(d, d, "+") - 2 * C
    function(p) {
      a <- exp(p[1])
      V <- exp(p[2]) / (2 * a) * exp(-a * gap) * (1 - exp(-2 * a * C))
      mvn_loglik(x, V)
    }
  }
  starts <- list(c(0, 0), log(c(5, 3)), log(c(0.5, 0.5)))
  # G5.4, G5.6, G5.6b: alpha = 2, sigma2 = 1, 80 replicates of 150 tips; the true pair must
  # fall outside the 95% likelihood region (2 df) at the nominal rate.
  set.seed(301)
  outside <- vapply(seq_len(80), function(i) {
    tree <- evol_tree(150)
    C <- ape::vcv.phylo(tree)
    x <- trait_ou(tree, alpha = 2, sigma2 = 1)[rownames(C), 1]
    f <- ou_loglik(x, C)
    2 * (max_loglik(f, starts) - f(c(log(2), 0))) > stats::qchisq(0.95, 2)
  }, logical(1))
  expect_true(rejection_compatible(sum(outside), 80))
  # G5.7: the median absolute log error of alpha falls with 50 against 300 tips (20 replicates
  # each); asymptotically by sqrt(6), we ask for a factor of 1.5.
  set.seed(302)
  alpha_error <- function(n) {
    stats::median(vapply(seq_len(20), function(i) {
      tree <- evol_tree(n)
      C <- ape::vcv.phylo(tree)
      x <- trait_ou(tree, alpha = 2, sigma2 = 1)[rownames(C), 1]
      f <- ou_loglik(x, C)
      fits <- lapply(starts, function(s) stats::optim(s, function(p) -f(p)))
      best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
      abs(best$par[1] - log(2))
    }, numeric(1)))
  }
  expect_lt(alpha_error(300), alpha_error(50) / 1.5)
  # the deprecated alias draws the same values
  tree <- evol_tree(20)
  set.seed(303)
  a <- trait_ou(tree, alpha = 2)
  set.seed(303)
  b <- expect_warning_value(traits_ou(tree, alpha = 2), "deprecated")
  expect_identical(a, b)
})

test_that("trait_rates() is unit-rate Brownian motion on the scaled tree and recovers sigma2", {
  skip_if_not_extended()
  # G5.6b: for each model and 5 seeds (200 tips), the contrasts of the tips on the tree scaled
  # by the rates are iid N(0, 1): the pooled sum of squares must lie within the 0.05% and
  # 99.95% chi-squared quantiles, and a Kolmogorov-Smirnov test must not reject at 0.001.
  for (model in c("shift", "gbm", "iid")) {
    set.seed(401)
    pics <- unlist(lapply(1:5, function(i) {
      x <- trait_rates(evol_tree(200), model, sigma2 = 2, shift_size = 5, rate_sd = 0.7)
      ape::pic(x$tips[x$scaled$tip.label], x$scaled)
    }))
    expect_gte(sum(pics^2), stats::qchisq(0.0005, length(pics)))
    expect_lte(sum(pics^2), stats::qchisq(0.9995, length(pics)))
    expect_gt(stats::ks.test(pics, "pnorm")$p.value, 0.001)
  }
  # G5.4: without shifts the rate is sigma2 everywhere, so the mean square of the contrasts on
  # the original tree recovers sigma2 = 2.5 within its chi-squared 99.9% interval, per seed.
  for (seed in 402:406) {
    set.seed(seed)
    tree <- evol_tree(200)
    x <- trait_rates(tree, "shift", sigma2 = 2.5, shifts = 0)
    ss <- sum(ape::pic(x$tips[tree$tip.label], tree)^2)
    expect_gte(ss / 2.5, stats::qchisq(0.0005, 199))
    expect_lte(ss / 2.5, stats::qchisq(0.9995, 199))
  }
})

test_that("history_rates() with trait_rates() gives Brownian motion with a rate per state", {
  skip_if_not_extended()
  # A habitat history from trait_mk(), rates 1 and 5 per state; the tips have covariance
  # sigma2_1 C_1 + sigma2_2 C_2, with C_k built from the time spent in state k. G5.4, G5.6,
  # G5.6b: 80 replicates of 120 tips; the true rates must fall outside the 95% likelihood
  # region (2 df) at the nominal rate.
  set.seed(501)
  outside <- vapply(seq_len(80), function(i) {
    tree <- evol_tree(120)
    h <- trait_mk(tree, q = 1, states = c("a", "b"))
    rates <- history_rates(h, "state", state_rates = c(a = 1, b = 5))
    x <- trait_rates(h, edge_rates = rates)$tips
    base <- tree
    base$edge <- h$edge
    base$tip.label <- h$tip.label
    Ck <- lapply(c("a", "b"), function(s) {
      base$edge.length <- h$mapped.edge[, s]
      ape::vcv.phylo(base)
    })
    x <- x[rownames(Ck[[1]])]
    f <- function(p) mvn_loglik(x, exp(p[1]) * Ck[[1]] + exp(p[2]) * Ck[[2]])
    best <- max_loglik(f, list(c(0, 0), c(0, log(5)), log(c(3, 3))))
    2 * (best - f(log(c(1, 5)))) > stats::qchisq(0.95, 2)
  }, logical(1))
  expect_true(rejection_compatible(sum(outside), 80))
})

# --- discrete traits ---------------------------------------------------------------------------

test_that("trait_mk() recovers q from its history and from the tips", {
  skip_if_not_extended()
  # Under the equal-rates model with two states every lineage leaves its state at rate q, so
  # the number of changes on a tree of total length L is Poisson(q L). G5.6b: 20 seeds of 100
  # tips with q = 0.8, each count within the 0.05% and 99.95% Poisson quantiles, and the
  # pooled count likewise.
  set.seed(601)
  counts <- t(vapply(seq_len(20), function(i) {
    tree <- evol_tree(100)
    h <- trait_mk(tree, q = 0.8)
    c(changes = sum(lengths(h$maps) - 1), length = sum(tree$edge.length))
  }, numeric(2)))
  mu <- 0.8 * counts[, "length"]
  expect_true(all(counts[, "changes"] >= stats::qpois(0.0005, mu)))
  expect_true(all(counts[, "changes"] <= stats::qpois(0.9995, mu)))
  expect_gte(sum(counts[, "changes"]), stats::qpois(0.0005, sum(mu)))
  expect_lte(sum(counts[, "changes"]), stats::qpois(0.9995, sum(mu)))
  # From the tips alone, by maximum likelihood with the root drawn with equal probabilities as
  # in trait_mk(). G5.4, G5.6: 100 replicates of 150 tips, q = 1; the true q must fall outside
  # the 95% likelihood interval at the nominal rate.
  mk_loglik <- function(tree, states) {
    tree <- ape::reorder.phylo(tree, "postorder")
    n <- ape::Ntip(tree)
    tip <- cbind(states[tree$tip.label] == "0", states[tree$tip.label] == "1") * 1
    function(q) {
      L <- matrix(1, n + tree$Nnode, 2)
      L[seq_len(n), ] <- tip
      log_scale <- 0
      for (e in seq_len(nrow(tree$edge))) {
        same <- 0.5 + 0.5 * exp(-2 * q * tree$edge.length[e])
        child <- L[tree$edge[e, 2], ]
        L[tree$edge[e, 1], ] <- L[tree$edge[e, 1], ] *
          c(same * child[1] + (1 - same) * child[2], (1 - same) * child[1] + same * child[2])
        if (e %% 20 == 0) {
          s <- max(L[tree$edge[e, 1], ])
          L[tree$edge[e, 1], ] <- L[tree$edge[e, 1], ] / s
          log_scale <- log_scale + log(s)
        }
      }
      log(sum(0.5 * L[n + 1, ])) + log_scale
    }
  }
  set.seed(602)
  outside <- vapply(seq_len(100), function(i) {
    tree <- evol_tree(150)
    f <- mk_loglik(tree, trait_mk(tree, q = 1)$states)
    best <- stats::optimize(function(lq) f(exp(lq)), log(c(0.01, 100)),
                            maximum = TRUE)$objective
    2 * (best - f(1)) > stats::qchisq(0.95, 1)
  }, logical(1))
  expect_true(rejection_compatible(sum(outside), 100))
})

# --- quantitative genetics ---------------------------------------------------------------------

test_that("drift_test() holds its size under drift with known G and has power under selection", {
  skip_if_not_extended()
  # 8 traits with a known G whose eigenvalues span two orders of magnitude; 30 populations of
  # 20 individuals with P = G. Under drift the means are N(0, 0.5 G); under selection
  # N(0, 0.5 G^2 / smallest eigenvalue), which makes the slope 2 and keeps the divergence on
  # every component well above the sampling error of the means.
  set.seed(701)
  q <- 8
  rotation <- qr.Q(qr(matrix(stats::rnorm(q * q), q)))
  values <- 2^seq(3, -3.6, length.out = q)
  G <- rotation %*% diag(values) %*% t(rotation)
  B_sel <- rotation %*% diag(values^2 / min(values)) %*% t(rotation)
  one <- function(B) {
    means <- matrix(stats::rnorm(30 * q), 30) %*% chol(0.5 * B)
    ind <- means[rep(1:30, each = 20), ] + matrix(stats::rnorm(600 * q), 600) %*% chol(G)
    drift_test(ind, group = rep(1:30, each = 20), W = G)
  }
  # Type I error: 300 replicates (validation_drift.R found 3.6% of 60 runs on clades).
  p_drift <- vapply(seq_len(300), function(i) one(G)$p_value, numeric(1))
  expect_true(rejection_compatible(sum(p_drift < 0.05), 300))
  # Power: 50 replicates; the slope is estimated with a standard error near 0.1 here, so a
  # slope of 2 should be rejected almost always (we ask for 90%).
  sel <- lapply(seq_len(50), function(i) one(B_sel))
  expect_gte(mean(vapply(sel, `[[`, numeric(1), "p_value") < 0.05), 0.9)
  expect_true(all(vapply(sel, `[[`, numeric(1), "slope") > 1))
})

test_that("maladaptation() follows the response to selection of Lande (1976)", {
  skip_if_not_extended()
  # 16 isolated sites (no dispersal, no seeds after the first step), one species with niche
  # breadth 10, G = 4 and no environmental variance; half the sites jump from 50 to 55 at step
  # 21. Without drift the deviation shrinks by a factor 1 - a each step, a = G / (sigma^2 + P),
  # exactly; with drift it does so on average, since the drift noise has mean zero.
  land <- landscape_grid(4, disp_scale = 0.5)
  tr <- data.frame(optimum = 50, sigma = 10, r_max = 4, row.names = "sp1")
  alpha <- matrix(0.01, dimnames = list("sp1", "sp1"))
  env <- env_pulse(rep(50, 16), n_steps = 40, sites = 1:8, start = 21, duration = Inf,
                   value = 55)
  a <- 4 / (100 + 4)
  run <- function(ne_ratio) {
    simulate_metacommunity(land, tr, alpha, dispersal = 0, n_steps = 40, env_time = env,
                           seed_steps = 1, record_steps = 20:40,
                           initial = matrix(100, 16, 1, dimnames = list(NULL, "sp1")),
                           evolution = evolution_control(4, ne_ratio = ne_ratio,
                                                         founder = FALSE))
  }
  # G5.4: deterministic recursion
  set.seed(801)
  m <- maladaptation(run(Inf), step = 21:40)
  m <- m[m$site %in% paste0("site", 1:8), ]
  d <- tapply(m$deviation, m$step, mean)
  expect_equal(as.vector(d[-1] / d[-length(d)]), rep(1 - a, length(d) - 1), tolerance = 1e-8)
  expect_equal(m$deviation, 5 * (1 - a)^(m$since_change + 1), tolerance = 1e-8)
  expect_equal(m$load, 1 - exp(-m$deviation^2 / (2 * (100 + 4))), tolerance = 1e-10)
  # G5.6b: with drift (Ne = 0.05 N), over 8 runs, the mean deviation of the 64 changed
  # populations 10 steps after the change is within 4 standard errors of 5 (1 - a)^11
  set.seed(802)
  dev <- unlist(lapply(1:8, function(i) {
    mm <- maladaptation(run(0.05), step = 31)
    mm$deviation[mm$site %in% paste0("site", 1:8)]
  }))
  expect_gt(length(dev), 60)
  expect_lt(abs(mean(dev) - 5 * (1 - a)^11), 4 * stats::sd(dev) / sqrt(length(dev)))
})

# --- metapopulations and tip rates -------------------------------------------------------------

test_that("persistence() recovers colonisation and extirpation and converges with more sites", {
  skip_if_not_extended()
  markov <- function(sites, times, c_true, e_true) {
    x <- array(0L, c(sites, 3, times))
    x[, , 1] <- stats::rbinom(3 * sites, 1, c_true / (c_true + e_true))
    for (t in 2:times) {
      x[, , t] <- stats::rbinom(3 * sites, 1, ifelse(x[, , t - 1] > 0, 1 - e_true, c_true))
    }
    x
  }
  # G5.4, G5.6, G5.6b: 5 seeds x 3 species, 100 sites, 30 intervals, c = 0.15 and e = 0.25 per
  # interval; each estimate is a binomial proportion, so it must lie within 4 binomial
  # standard errors of the truth (given the number of exposures).
  for (seed in 901:905) {
    set.seed(seed)
    x <- markov(100, 31, 0.15, 0.25)
    p <- persistence(x)
    absent <- rowSums(colSums(x[, , -31] <= 0))
    present <- rowSums(colSums(x[, , -31] > 0))
    expect_true(all(abs(p$colonisation - 0.15) < 4 * sqrt(0.15 * 0.85 / absent)))
    expect_true(all(abs(p$extirpation - 0.25) < 4 * sqrt(0.25 * 0.75 / present)))
  }
  # G5.7: the root mean squared error of epsilon falls with 20 against 320 sites (30 replicates
  # each); asymptotically by a factor of 4, we ask for 2.5.
  set.seed(906)
  rmse <- function(sites) {
    truth <- stats::qlogis(0.15 / 0.4)
    sqrt(mean(unlist(lapply(seq_len(30), function(i) {
      (persistence(markov(sites, 21, 0.15, 0.25))$epsilon - truth)^2
    }))))
  }
  expect_lt(rmse(320), rmse(20) / 2.5)
})

test_that("dr_rates() recovers the speciation rate of pure-birth trees", {
  skip_if_not_extended()
  # The equal-splits measures of the tips of a binary tree add up to its total length L (the
  # weights 2^-(j-1) of the tips below each edge add up to one), so the harmonic mean of DR is
  # n / L. Under pure birth, L is close to a sum of n - 1 exponential waiting times of rate
  # lambda (k lineages for an Exp(k lambda) time), so lambda * L is about Gamma(n - 1, 1).
  # G5.4, G5.6b: 20 trees of 200 tips from ape::rphylo() with lambda = 1.5; the pooled
  # lambda * L must lie within 4 standard deviations of its mean (one waiting time more or
  # less shifts it by less than 0.3 standard deviations).
  set.seed(1001)
  n <- 200
  total <- vapply(seq_len(20), function(i) {
    tree <- ape::rphylo(n, 1.5, 0)
    dr <- dr_rates(tree)
    expect_equal(1 / mean(1 / dr), n / sum(tree$edge.length), tolerance = 1e-10)
    n * mean(1 / dr)  # the total length L
  }, numeric(1))
  k <- 20 * (n - 1)
  expect_lt(abs(1.5 * sum(total) - k), 4 * sqrt(k))
  # G5.7: the relative error of the harmonic mean of DR falls with tree size, as 1 / sqrt(n):
  # 30 trees of 50 against 800 tips, a factor of 4 asymptotically; we ask for 2.5.
  set.seed(1002)
  rel_error <- function(n) {
    sqrt(mean(vapply(seq_len(30), function(i) {
      (1 / mean(1 / dr_rates(ape::rphylo(n, 1.5, 0))) / 1.5 - 1)^2
    }, numeric(1))))
  }
  expect_lt(rel_error(800), rel_error(50) / 2.5)
})
