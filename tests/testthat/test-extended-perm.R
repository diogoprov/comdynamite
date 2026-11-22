# rOpenSci statistical standards tested in this file (see the srr tags in R/)
#' @srrstats {RE1.4} Tested in this file, as described in R/srr-stats-regression.R.
NULL

# Extended tests of the permutation and null-model tests: pcc_test(), phylo_beta_test(),
# pcps_test(), node_test() with node_contrasts(), phylo_vector(), trait_partition(), and the
# restricted nulls torus_permutations(), rigid_permutations(), patch_permutations() and
# surrogates of the predictors (null = "variogram"). Run with COMDYNAMITE_EXTENDED_TESTS=true.
#
# Tolerances. A permutation test whose null hypothesis holds rejects at the 5% level with
# probability at most 0.05, and exactly 0.05 up to the discreteness of 99 permutations
# (rejection when p <= 0.05, that is at most 4 permutations of 99 reach the observed value,
# probability 5/100). Replicates are independent data sets, so the number of rejections is
# binomial: `rejection_compatible()` asks that 0.05 lies in the exact 99% interval of the
# observed rate. The larger of two p-values (`p_max`) is conservative by construction when
# one null hypothesis is false and the other true, so for it, and for restricted nulls whose
# claim is validity rather than exactness, we only ask that the rate is not above 0.05
# (`not_inflated()`: the lower end of the 99% interval is at most 0.05). Power is checked
# with the lower end of the 99% interval against a stated floor. Only the p-value of the root
# node enters the tests of node_test(), so that each replicate gives one independent value.
#
# Set COMDYNAMITE_EXTENDED_VERBOSE=true to print the observed rates.

not_inflated <- function(rejections, n, level = 0.05) {
  stats::binom.test(rejections, n, conf.level = 0.99)$conf.int[1] <= level
}
at_least <- function(rejections, n, floor) {
  stats::binom.test(rejections, n, conf.level = 0.99)$conf.int[1] >= floor
}
note <- function(label, x) {
  if (identical(Sys.getenv("COMDYNAMITE_EXTENDED_VERBOSE"), "true")) {
    message(label, ": ", paste(names(x), signif(x, 3), sep = " = ", collapse = ", "))
  }
}
make_tree <- function(n) {
  tree <- if (requireNamespace("treeshapesim", quietly = TRUE)) {
    treeshapesim::sim_tree_shape(n)
  } else {
    ape::rcoal(n)
  }
  tree$tip.label <- paste0("sp", seq_len(n))
  tree
}
rejections <- function(p) sum(p <= 0.05, na.rm = TRUE)
n_valid <- function(p) sum(!is.na(p))

# Data of Table 1 of Peres-Neto & Kembel (2015), as in data-raw/validation_pcc.R: Gaussian
# responses to a uniform gradient, incidence data. With `phylogeny`, optima evolve on the tree
# (lambda = 1, the null hypothesis of the permutation of species is false); without it they do
# not (lambda = 0, true). With `environment`, the gradient that made the data is analysed (the
# null hypothesis of the permutation of sites is false); without it, an unrelated uniform
# predictor is (true).
whole_tree_data <- function(n_sites, n_species, environment, phylogeny, sigma_mean = 10) {
  repeat {
    tree <- make_tree(n_species)
    traits <- species_traits(tree, lambda = as.numeric(phylogeny), optima_range = c(-1, 101),
                             sigma_mean = sigma_mean, sigma_sd = 10, rmax_range = c(0, 30))
    env <- stats::runif(n_sites, 0, 100)
    comm <- simulate_static(traits, env, incidence = TRUE)$abundance
    used <- if (environment) env else stats::runif(n_sites, 0, 100)
    keep <- rowSums(comm) > 0
    comm <- comm[keep, colSums(comm[keep, , drop = FALSE]) > 0, drop = FALSE]
    if (nrow(comm) >= 10 && ncol(comm) >= 5) break
  }
  list(comm = comm, tree = ape::keep.tip(tree, colnames(comm)), env = used[keep])
}

whole_tree_p <- function(d, n_perm = 99) {
  pcc <- pcc_test(d$comm, d$tree, d$env, n_perm = n_perm)
  beta <- phylo_beta_test(d$comm, d$tree, d$env, n_perm = n_perm)
  pcps <- pcps_test(d$comm, d$tree, d$env, axes = 1, n_perm = n_perm)$tests
  c(pcc_site = pcc$p_site, pcc_species = pcc$p_species, pcc_max = pcc$p_max,
    beta_site = beta$p_site, beta_species = beta$p_species, beta_max = beta$p_max,
    pcps_site = pcps$p_site[1], pcps_species = pcps$p_species[1], pcps_max = pcps$p_max[1],
    pcps1_site = pcps$p_site[2], pcps1_species = pcps$p_species[2], pcps1_max = pcps$p_max[2])
}

test_that("whole-tree tests hold the level when neither null hypothesis is false", {
  skip_if_not_extended()
  # Type I error of the permutation of sites and of species, and of p_max, for pcc_test(),
  # phylo_beta_test() and pcps_test() (whole dissimilarity and first axis), with free
  # permutations. 200 data sets of 60 sites and 24 species, 99 permutations.
  set.seed(4101)
  n_rep <- 200
  p <- replicate(n_rep, whole_tree_p(whole_tree_data(60, 24, FALSE, FALSE)))
  rates <- apply(p, 1, function(x) rejections(x) / n_valid(x))
  note("neither", rates)
  for (name in grep("_(site|species)$", rownames(p), value = TRUE)) {
    expect_true(rejection_compatible(rejections(p[name, ]), n_valid(p[name, ])), label = name)
  }
  for (name in grep("_max$", rownames(p), value = TRUE)) {
    expect_true(not_inflated(rejections(p[name, ]), n_valid(p[name, ])), label = name)
  }
})

test_that("whole-tree tests hold the level of the species test and have power for sites", {
  skip_if_not_extended()
  # Type I error of the permutation of species when composition follows the analysed
  # gradient without phylogenetic signal in the optima, with power of the permutation of
  # sites for the same data; p_max not above 0.05. 200 data sets, 60 sites, 24 species, 99
  # permutations. Floor of power: 0.9, since the gradient made the data (data-raw/
  # validation_pcc.R and validation_phylo_beta.R gave 100% with 100 sites).
  set.seed(4102)
  n_rep <- 200
  p <- replicate(n_rep, whole_tree_p(whole_tree_data(60, 24, TRUE, FALSE)))
  rates <- apply(p, 1, function(x) rejections(x) / n_valid(x))
  note("environment", rates)
  for (name in grep("_species$", rownames(p), value = TRUE)) {
    expect_true(rejection_compatible(rejections(p[name, ]), n_valid(p[name, ])), label = name)
  }
  for (name in grep("_max$", rownames(p), value = TRUE)) {
    expect_true(not_inflated(rejections(p[name, ]), n_valid(p[name, ])), label = name)
  }
  for (name in c("pcc_site", "beta_site", "pcps_site")) {
    expect_true(at_least(rejections(p[name, ]), n_valid(p[name, ]), 0.9), label = name)
  }
})

test_that("whole-tree tests reject when both the gradient and the phylogeny act", {
  skip_if_not_extended()
  # Power of p_max when optima evolve on the tree and the gradient made the data. 100 data
  # sets of 100 sites and 40 species, 99 permutations. data-raw/validation_pcc.R and
  # validation_phylo_beta.R gave 85% to 91% with 100 sites and 50 species, so for those two
  # the floor of the 99% interval is 0.5, well below those rates. pcps_test() has no stored
  # power for this scenario, so its floor is 0.2, four times the level.
  set.seed(4103)
  n_rep <- 100
  p <- replicate(n_rep, whole_tree_p(whole_tree_data(100, 40, TRUE, TRUE)))
  rates <- apply(p, 1, function(x) rejections(x) / n_valid(x))
  note("both", rates)
  floors <- c(pcc_max = 0.5, beta_max = 0.5, pcps_max = 0.2)
  for (name in names(floors)) {
    expect_true(at_least(rejections(p[name, ]), n_valid(p[name, ]), floors[[name]]), label = name)
  }
})

test_that("power of the whole-tree tests grows with the number of sites and of species", {
  skip_if_not_extended()
  # G5.7, in two parts, for pcc_test(), phylo_beta_test() and pcps_test() (whole
  # dissimilarity). (1) Power of the permutation of sites grows with sites (15, 30, 60) when
  # the analysed predictor is the gradient plus normal noise of standard deviation 100, a
  # weak relation, 20 species. (2) Power of p_max grows with species (8, 16, 32) when optima
  # evolve on the tree and the gradient made the data, 60 sites: with the gradient analysed,
  # the permutation of sites rejects almost always, so p_max is limited by the permutation
  # of species, whose information is the number of species and not the number of sites.
  # 80 data sets per size, 99 permutations. The largest size must have more power than the
  # smallest (one-sided test of two proportions at the 1% level), and no step may fall by
  # more than 2.58 standard errors of a difference of two proportions near 0.5.
  set.seed(4107)
  n_rep <- 80
  three <- function(d, column) {
    c(pcc = pcc_test(d$comm, d$tree, d$env, n_perm = 99)[[column]],
      beta = phylo_beta_test(d$comm, d$tree, d$env, n_perm = 99)[[column]],
      pcps = pcps_test(d$comm, d$tree, d$env, axes = 1, n_perm = 99)$tests[[column]][1])
  }
  sites <- vapply(c(15, 30, 60), function(n) {
    rowMeans(replicate(n_rep, {
      d <- whole_tree_data(n, 20, TRUE, FALSE)
      d$env <- d$env + stats::rnorm(length(d$env), 0, 100)
      three(d, "p_site")
    }) <= 0.05)
  }, numeric(3))
  species <- vapply(c(8, 16, 32), function(k) {
    rowMeans(replicate(n_rep, three(whole_tree_data(60, k, TRUE, TRUE), "p_max")) <= 0.05)
  }, numeric(3))
  colnames(sites) <- paste0("sites_", c(15, 30, 60))
  colnames(species) <- paste0("species_", c(8, 16, 32))
  se <- sqrt(2 * 0.25 / n_rep)
  for (test in rownames(sites)) {
    note(paste("power of p_site", test), sites[test, ])
    note(paste("power of p_max", test), species[test, ])
    for (power in list(sites[test, ], species[test, ])) {
      x <- round(power * n_rep)
      expect_lt(stats::prop.test(x[c(3, 1)], c(n_rep, n_rep), alternative = "greater")$p.value,
                0.01, label = test)
      expect_true(all(diff(power) > -2.58 * se), label = test)
    }
  }
})

# Contrasts of the root node, with optima on each side of the root drawn around 30 and 70
# (effect) or from one uniform distribution (no effect)
root_node_data <- function(n_sites, n_species, environment, clade_effect) {
  repeat {
    tree <- make_tree(n_species)
    side <- node_contrasts(tree)[, 1]
    optimum <- if (clade_effect) {
      50 + 20 * sign(side) + stats::rnorm(n_species, 0, 10)
    } else {
      stats::runif(n_species, -1, 101)
    }
    traits <- data.frame(optimum = optimum, sigma = 20, r_max = stats::runif(n_species, 0, 30),
                         row.names = tree$tip.label)
    env <- stats::runif(n_sites, 0, 100)
    comm <- simulate_static(traits, env, incidence = TRUE)$abundance
    used <- if (environment) env else stats::runif(n_sites, 0, 100)
    keep <- rowSums(comm) > 0
    comm <- comm[keep, , drop = FALSE]
    nodes <- node_contrasts(tree)[, 1, drop = FALSE]
    present <- colSums(comm) > 0
    if (nrow(comm) >= 10 && any(nodes[present, 1] > 0) && any(nodes[present, 1] < 0)) break
  }
  list(comm = comm, nodes = nodes, env = used[keep])
}

test_that("node_test() holds the level at the root node and detects a clade effect", {
  skip_if_not_extended()
  # Type I error of the permutation of sites (unrelated predictor) and of species (gradient
  # analysed, optima independent of the tree) at the root node, and power of p_max when the
  # two clades of the root have optima around 30 and 70. 150 data sets per scenario, 60 sites,
  # 24 species, 99 permutations. Floor of power 0.8: the clades differ by two standard
  # deviations of their optima.
  set.seed(4104)
  n_rep <- 150
  run <- function(environment, clade_effect) {
    replicate(n_rep, {
      d <- root_node_data(60, 24, environment, clade_effect)
      unlist(node_test(d$comm, d$nodes, d$env, n_perm = 99)$nodes[1, c("p_site", "p_species",
                                                                        "p_max")])
    })
  }
  sites <- run(FALSE, FALSE)
  species <- run(TRUE, FALSE)
  clade <- run(TRUE, TRUE)
  note("node site null", c(p_site = rejections(sites["p_site", ]) / n_valid(sites["p_site", ]),
                           p_max = rejections(sites["p_max", ]) / n_valid(sites["p_max", ])))
  note("node species null", c(p_species = rejections(species["p_species", ]) /
                                n_valid(species["p_species", ]),
                              p_max = rejections(species["p_max", ]) /
                                n_valid(species["p_max", ])))
  note("node clade effect", c(p_max = rejections(clade["p_max", ]) / n_valid(clade["p_max", ])))
  expect_true(rejection_compatible(rejections(sites["p_site", ]), n_valid(sites["p_site", ])))
  expect_true(rejection_compatible(rejections(species["p_species", ]),
                                   n_valid(species["p_species", ])))
  expect_true(not_inflated(rejections(sites["p_max", ]), n_valid(sites["p_max", ])))
  expect_true(not_inflated(rejections(species["p_max", ]), n_valid(species["p_max", ])))
  expect_true(at_least(rejections(clade["p_max", ]), n_valid(clade["p_max", ]), 0.8))
})

test_that("phylo_vector() holds the level of its global test and has power", {
  skip_if_not_extended()
  # Type I error of the global test: permutation of sites with an unrelated predictor and
  # of species with optima independent of the tree, and power of p_max with optima that
  # evolve on the tree and the gradient analysed. 200 data sets per scenario, 60 sites, 24
  # species, 99 permutations; power with 100 sites and 40 species and a floor of 0.5, as for
  # the whole-tree tests.
  set.seed(4105)
  global <- function(d) {
    out <- phylo_vector(d$comm, d$tree, d$env, n_perm = 99, n_sim = 10)
    c(p_site = out$p_site, p_species = out$p_species, p_max = out$p_max)
  }
  neither <- replicate(200, global(whole_tree_data(60, 24, FALSE, FALSE)))
  environment <- replicate(200, global(whole_tree_data(60, 24, TRUE, FALSE)))
  both <- replicate(100, global(whole_tree_data(100, 40, TRUE, TRUE)))
  note("phylo_vector neither", rowMeans(neither <= 0.05))
  note("phylo_vector environment", rowMeans(environment <= 0.05))
  note("phylo_vector both", rowMeans(both <= 0.05))
  expect_true(rejection_compatible(rejections(neither["p_site", ]), 200))
  expect_true(rejection_compatible(rejections(neither["p_species", ]), 200))
  expect_true(rejection_compatible(rejections(environment["p_species", ]), 200))
  expect_true(not_inflated(rejections(neither["p_max", ]), 200))
  expect_true(not_inflated(rejections(environment["p_max", ]), 200))
  expect_true(at_least(rejections(both["p_max", ]), 100, 0.5))
})

# Generator of data-raw/validation_fractions.R: a grid of 8 x 8 sites, `space` with the six
# broadest Moran's eigenvector maps, and Poisson counts of Gaussian responses to a gradient
# that combines env (weight alpha) and a smooth field in the span of `space` (weight beta).
# With `hidden`, composition follows another trait than the one analysed.
fraction_data <- function(alpha, beta, mems, hidden = FALSE, k = 20) {
  n_sites <- nrow(mems)
  field <- function(cols) {
    as.vector(scale(mems[, cols, drop = FALSE] %*% stats::rnorm(length(cols))))
  }
  repeat {
    env <- stats::rnorm(n_sites)
    g <- if (alpha == 0 && beta == 0) {
      rep(0, n_sites)
    } else {
      (alpha * env + beta * field(1:3)) / sqrt(alpha^2 + beta^2)
    }
    trait <- stats::setNames(stats::runif(k, 0, 100), paste0("sp", seq_len(k)))
    driver <- if (hidden) stats::setNames(stats::runif(k, 0, 100), names(trait)) else trait
    lambda <- 5 * exp(-outer(50 + 20 * g, driver, "-")^2 / (2 * 15^2)) + 0.3
    comm <- matrix(stats::rpois(length(lambda), lambda), nrow(lambda),
                   dimnames = list(NULL, names(trait)))
    if (all(rowSums(comm) > 0) && sum(colSums(comm) > 0) >= 5) break
  }
  list(comm = comm, trait = trait, env = cbind(e = env))
}

test_that("trait_partition() holds the level of the fractions a and c and has power for a", {
  skip_if_not_extended()
  skip_if_not_installed("adespatial")
  # Type I error of the permutation of sites for `a` (composition follows space only), for
  # `c` (composition follows an unstructured env only) and of the permutation of species for
  # `a` (composition follows env through another trait), and power of p_max of `a` when env
  # acts through the analysed trait. 150 data sets per scenario, 99 permutations. Floor of
  # power 0.8: data-raw/validation_fractions.R gave complete power in this scenario.
  set.seed(4106)
  n_rep <- 150
  coords <- as.matrix(expand.grid(x = 1:8, y = 1:8))
  mems <- as.matrix(adespatial::dbmem(coords, silent = TRUE))
  space <- mems[, 1:6]
  run <- function(alpha, beta, hidden = FALSE) {
    replicate(n_rep, {
      d <- fraction_data(alpha, beta, mems, hidden)
      f <- trait_partition(d$comm, d$trait, d$env, space, n_perm = 99)$fractions
      rownames(f) <- f$fraction
      c(a_site = f["a", "p_site"], a_species = f["a", "p_species"], a_max = f["a", "p_max"],
        c_site = f["c", "p_site"])
    })
  }
  space_only <- run(0, 1)
  env_only <- run(1, 0)
  hidden <- run(1, 0, hidden = TRUE)
  note("partition space only", rowMeans(space_only <= 0.05))
  note("partition env only", rowMeans(env_only <= 0.05))
  note("partition hidden trait", rowMeans(hidden <= 0.05))
  expect_true(rejection_compatible(rejections(space_only["a_site", ]), n_rep))
  expect_true(rejection_compatible(rejections(env_only["c_site", ]), n_rep))
  expect_true(rejection_compatible(rejections(hidden["a_species", ]), n_rep))
  expect_true(not_inflated(rejections(hidden["a_max", ]), n_rep))
  expect_true(at_least(rejections(env_only["a_max", ]), n_rep, 0.8))
})

# Gaussian random field on a grid of side x side sites, standardised, with opposite edges
# adjacent (periodic) or not
grid_field <- function(xy, side, bandwidth, periodic) {
  dx <- abs(outer(xy[, 1], xy[, 1], "-"))
  dy <- abs(outer(xy[, 2], xy[, 2], "-"))
  if (periodic) {
    dx <- pmin(dx, side - dx)
    dy <- pmin(dy, side - dy)
  }
  z <- drop(exp(-(dx^2 + dy^2) / (2 * bandwidth^2)) %*% stats::rnorm(nrow(xy)))
  (z - mean(z)) / stats::sd(z)
}

# Communities that follow one autocorrelated field and a predictor that is another,
# independent field with the same autocorrelation (data-raw/validation_pcps.R)
field_data <- function(xy, side, bandwidth, periodic, n_species = 20) {
  repeat {
    tree <- make_tree(n_species)
    traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100), rmax_range = c(0, 30))
    hidden <- grid_field(xy, side, bandwidth, periodic)
    env <- grid_field(xy, side, bandwidth, periodic)
    comm <- simulate_static(traits, (hidden - min(hidden)) / diff(range(hidden)) * 100,
                            incidence = TRUE)$abundance
    present <- colSums(comm) > 0
    if (sum(present) >= 5 && sum(rowSums(comm) > 0) >= 20) break
  }
  comm <- comm[, present, drop = FALSE]
  list(comm = comm, tree = ape::keep.tip(tree, colnames(comm)), env = env,
       occupied = rowSums(comm) > 0)
}

test_that("torus_permutations() hold the level with autocorrelated fields on a torus", {
  skip_if_not_extended()
  # Type I error of the permutation of sites with an autocorrelated predictor unrelated to
  # communities that follow another autocorrelated field, both periodic on a grid of 10 x 10
  # sites (bandwidth 1.5), where the moves of the torus keep the autocorrelation exactly. The
  # torus must hold the level in pcps_test(), pcc_test() and phylo_beta_test(); free
  # permutations of the same data must reject more than 5% (data-raw/validation_pcps.R:
  # 52-65%), which shows the data carry the problem the torus solves. 150 data sets, 99
  # permutations.
  set.seed(4108)
  n_rep <- 150
  side <- 10
  xy <- as.matrix(expand.grid(x = seq_len(side), y = seq_len(side)))
  p <- replicate(n_rep, {
    d <- field_data(xy, side, 1.5, TRUE)
    moves <- torus_permutations(side, n_perm = 99)
    o <- d$occupied
    c(pcps_free = pcps_test(d$comm[o, ], d$tree, d$env[o], axes = 1, n_perm = 99)$tests$p_site[1],
      pcps_torus = pcps_test(d$comm, d$tree, d$env, axes = 1, n_perm = 99,
                             site_permutations = moves)$tests$p_site[1],
      pcc_torus = pcc_test(d$comm, d$tree, d$env, n_perm = 99, site_permutations = moves)$p_site,
      beta_torus = phylo_beta_test(d$comm, d$tree, d$env, n_perm = 99,
                                   site_permutations = moves)$p_site)
  })
  rates <- rowMeans(p <= 0.05)
  note("torus, periodic fields", rates)
  expect_true(at_least(rejections(p["pcps_free", ]), n_rep, 0.05 + 1e-9))
  for (name in c("pcps_torus", "pcc_torus", "beta_torus")) {
    expect_true(rejection_compatible(rejections(p[name, ]), n_rep), label = name)
  }
})

test_that("variogram surrogates hold the level with autocorrelated fields with edges", {
  skip_if_not_extended()
  # Type I error of the site test with null = "variogram" in pcps_test() and
  # phylo_beta_test(), with fields that stop at the edges of a grid of 10 x 10 sites
  # (bandwidth 1.5), where the torus rejects too often (data-raw/validation_pcps.R: 10%).
  # Surrogates approximate the autocorrelation, so the claim is validity: the rate must not
  # be above 0.05. The torus rate on the same data is printed for comparison. 120 data sets,
  # 99 surrogates.
  set.seed(4109)
  n_rep <- 120
  side <- 10
  xy <- as.matrix(expand.grid(x = seq_len(side), y = seq_len(side)))
  p <- replicate(n_rep, {
    d <- field_data(xy, side, 1.5, FALSE)
    o <- d$occupied
    c(pcps_variogram = pcps_test(d$comm[o, ], d$tree, d$env[o], axes = 1, n_perm = 99,
                                 null = "variogram", coords = xy[o, ])$tests$p_site[1],
      beta_variogram = phylo_beta_test(d$comm[o, ], d$tree, d$env[o], n_perm = 99,
                                       null = "variogram", coords = xy[o, ])$p_site,
      pcps_torus = pcps_test(d$comm, d$tree, d$env, axes = 1, n_perm = 99,
                             site_permutations = torus_permutations(side, n_perm = 99)
                             )$tests$p_site[1])
  })
  rates <- rowMeans(p <= 0.05)
  note("fields with edges", rates)
  for (name in c("pcps_variogram", "beta_variogram")) {
    expect_true(not_inflated(rejections(p[name, ]), n_rep), label = name)
  }
})

# Two blocks of 5 x 5 sites placed at random without overlap (data-raw/validation_edges.R)
random_blocks <- function(xy) {
  repeat {
    centres <- xy[sample(nrow(xy), 2), ]
    one <- pmax(abs(xy[, 1] - centres[1, 1]), abs(xy[, 2] - centres[1, 2])) <= 2
    two <- pmax(abs(xy[, 1] - centres[2, 1]), abs(xy[, 2] - centres[2, 2])) <= 2
    if (!any(one & two) && sum(one) >= 15 && sum(two) >= 15) return(one * 1 - two * 1)
  }
}

test_that("patch and rigid permutations hold the level of biogeographic events", {
  skip_if_not_extended()
  # Type I error of the test of an event at the root node of node_test() when communities
  # follow a hidden autocorrelated environment and the event has no effect, in a square of 14
  # x 14 sites with edges, as in data-raw/validation_edges.R and validation_rigid.R (3-6% and
  # 4-5% there, against 44-58% with free permutations). Restricted nulls claim validity, so
  # their rates must not be above 0.05, and free permutations must be above it. Power with an
  # effect of 40 (the event moves the environment of its sites by 40% of its range) must be
  # above the level for both nulls. 120 data sets per scenario, 99 permutations.
  set.seed(4110)
  n_rep <- 120
  xy <- as.matrix(expand.grid(1:14, 1:14))
  d2 <- as.matrix(stats::dist(xy))^2
  run <- function(effect) {
    replicate(n_rep, {
      repeat {
        tree <- make_tree(12)
        hidden <- drop(exp(-d2 / (2 * 3^2)) %*% stats::rnorm(nrow(xy)))
        hidden <- (hidden - min(hidden)) / diff(range(hidden)) * 100
        event <- cbind(event = random_blocks(xy))
        traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100),
                                 rmax_range = c(0, 30))
        comm <- simulate_static(traits, hidden + effect * event[, 1], incidence = TRUE)$abundance
        comm <- comm[, colSums(comm) > 0, drop = FALSE]
        if (ncol(comm) < 4) next
        nodes <- node_contrasts(ape::keep.tip(tree, colnames(comm)))[, 1, drop = FALSE]
        break
      }
      noise <- cbind(noise = stats::rnorm(nrow(xy)))
      p <- function(...) {
        node_test(comm, nodes, noise, event, n_perm = 99, ...)$biogeography$p_partial_site[1, 1]
      }
      c(free = p(), patch = p(area_permutations = patch_permutations(xy, event, 99)),
        rigid = p(area_permutations = suppressWarnings(rigid_permutations(xy, event, 99))))
    })
  }
  null <- run(0)
  effect <- run(40)
  note("events, no effect", rowMeans(null <= 0.05, na.rm = TRUE))
  note("events, effect 40", rowMeans(effect <= 0.05, na.rm = TRUE))
  expect_true(at_least(rejections(null["free", ]), n_valid(null["free", ]), 0.05 + 1e-9))
  for (name in c("patch", "rigid")) {
    expect_true(not_inflated(rejections(null[name, ]), n_valid(null[name, ])), label = name)
    expect_true(at_least(rejections(effect[name, ]), n_valid(effect[name, ]), 0.05 + 1e-9),
                label = name)
  }
})
