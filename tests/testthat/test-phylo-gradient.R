# worked example of Appendix A of Peres-Neto et al. (2012)
example_a <- function() {
  comm <- rbind(c(1, 1, 0, 0, 0, 0), c(1, 1, 0, 0, 0, 0), c(1, 1, 0, 0, 0, 0), c(1, 0, 1, 0, 0, 0),
                c(1, 0, 1, 0, 0, 0), c(0, 0, 0, 1, 1, 0), c(0, 0, 0, 1, 0, 1), c(0, 0, 0, 1, 1, 0),
                c(0, 0, 0, 1, 1, 0), c(0, 0, 0, 0, 1, 1))
  colnames(comm) <- paste0("sp", 1:6)
  list(comm = comm, trait = setNames(c(0.9, 1.3, 1.1, 3.2, 3.4, 3.7), colnames(comm)),
       env = c(1.9, 1.9, 1.1, 1.5, 2, 3.1, 3.8, 3.1, 3.6, 3.7),
       space = c(0.1, 0.3, 0.5, 0.5, 0.8, 0.7, 0.2, 0.4, 0.7, 0.5))
}

simulated <- function(seed = 1, n = 50, k = 20) {
  set.seed(seed)
  tree <- ape::rcoal(k)
  env <- runif(n, 0, 100)
  optima <- ape::rTraitCont(tree)
  optima <- 100 * (optima - min(optima)) / diff(range(optima))
  comm <- t(exp(-outer(optima, env, "-")^2 / (2 * 15^2)))
  colnames(comm) <- names(optima)
  list(tree = tree, env = env, comm = comm)
}

test_that("standardised values match the appendix", {
  a <- example_a()
  expect_equal(round(weighted_standardise(as.matrix(a$trait), colSums(a$comm))[, 1], 3),
               c(-1.125, -0.784, -0.954, 0.835, 1.005, 1.261), ignore_attr = TRUE)
  expect_equal(round(weighted_standardise(as.matrix(a$env), rowSums(a$comm))[1:3, 1], 3),
               c(-0.709, -0.709, -1.555))
})

test_that("the correlation and fractions match the original code in Octave", {
  a <- example_a()
  expect_equal(abs(corner_correlation(a$comm, as.matrix(a$trait), as.matrix(a$env))[1]), 0.938977,
               tolerance = 1e-5)
  out <- trait_partition(a$comm, a$trait, a$env, a$space, n_perm = 3)
  f <- out$fractions
  expect_equal(f$r2[1:3], c(0.8839515378, 0.8816771038, 0.0107264733), tolerance = 1e-9)
  expect_equal(f$r2_adjusted[c(4, 6)], c(0.9637275518, -0.0160919076), tolerance = 1e-9)
  expect_equal(out$component_share, 19.630 / 20, tolerance = 1e-3)
  expect_equal(sum(f$r2[f$fraction %in% c("a", "b", "c", "d")]), 1)
  expect_true(all(is.na(f$p_max[f$fraction %in% c("b", "d")])))
  expect_output(print(out), "mean component")
})

test_that("mean and variance shares add to one", {
  a <- example_a()
  mean_part <- trait_partition(a$comm, a$trait, a$env, a$space, n_perm = 1)
  dispersion <- trait_partition(a$comm, a$trait, a$env, a$space, "variance", n_perm = 1)
  expect_equal(mean_part$component_share + dispersion$component_share, 1)
  expect_identical(dispersion$component, "variance")
})

test_that("several traits add to the total variation", {
  a <- example_a()
  two <- cbind(a$trait, rev(a$trait))
  out <- trait_partition(a$comm, two, a$env, a$space, n_perm = 1)
  expect_lte(out$fractions$r2[1], 1)
  expect_equal(out$component_share, trait_partition(a$comm, two[, c(1, 1)], a$env, a$space,
                                                    n_perm = 1)$component_share / 2 +
                 trait_partition(a$comm, two[, c(2, 2)], a$env, a$space,
                                 n_perm = 1)$component_share / 2)
})

test_that("nulls have one row per permutation", {
  a <- example_a()
  set.seed(1)
  out <- corner_nulls(a$comm, as.matrix(a$trait), as.matrix(a$env), as.matrix(a$env),
                      t(replicate(2, sample.int(nrow(a$comm)))), 2)
  expect_equal(dim(out$site), c(2, 1))
  expect_equal(dim(out$taxa), c(2, 1))
})

test_that("weighted coordinates are centred and orthogonal under the weights", {
  distances <- matrix(c(0, 2, 4, 4, 2, 0, 4, 4, 4, 4, 0, 2, 4, 4, 2, 0), 4)
  weights <- c(5, 3, 2, 4)
  out <- weighted_pcoa(distances^2 / 2, weights)
  expect_equal(out$values, c(3.4772239, 0.5329868, 0.3775444), tolerance = 1e-6)
  expect_equal(colSums(weights * out$vectors), rep(0, 3))
  expect_equal(crossprod(out$vectors, weights / 14 * out$vectors), diag(out$values))
})

test_that("the phylogenetic vector follows a phylogenetic response", {
  d <- simulated()
  set.seed(2)
  out <- phylo_vector(d$comm, d$tree, d$env, n_perm = 99, n_sim = 50)
  expect_s3_class(out, "mc_phylo_vector")
  expect_lte(out$p_max, 0.05)
  expect_gt(length(out$selected), 0)
  expect_identical(names(out$pev), colnames(d$comm))
  expect_gte(abs(out$correlation[["pev"]]), 0.3)
  expect_output(print(out), "PEV-B")
  # signed weights never fit worse than absolute ones, for orthogonal eigenvectors
  sel <- out$selected
  absolute <- out$eigenvectors[, sel, drop = FALSE] %*% abs(out$eigenvector_correlation[sel])
  expect_gte(abs(out$correlation[["pev"]]),
             abs(corner_correlation(d$comm, absolute, as.matrix(d$env))[1]) - 1e-12)
  x <- trait_partition(d$comm, out$pev, d$env, d$env + rnorm(length(d$env)), n_perm = 3)
  expect_gt(x$fractions$r2[1], 0)
})

test_that("no vector comes from a predictor without relation", {
  d <- simulated()
  set.seed(3)
  out <- phylo_vector(d$comm, d$tree, rnorm(length(d$env)), n_perm = 99, level = 0.001, n_sim = 5)
  expect_null(out$pev)
  expect_length(out$selected, 0)
  expect_output(print(out), "0 retained")
})

test_that("absent species are left out", {
  d <- simulated()
  d$comm[, 1] <- 0
  out <- phylo_vector(d$comm, d$tree, d$env, n_perm = 19, n_sim = 5,
                      site_permutations = t(replicate(19, sample(nrow(d$comm)))))
  expect_equal(nrow(out$eigenvectors), ncol(d$comm) - 1)
  expect_error(phylo_vector(d$comm, d$tree, d$env, n_perm = 9, site_permutations = cbind(1:3)))
})

partition_data <- function() {
  set.seed(2)
  trait <- cbind(t1 = seq(0, 100, 10), t2 = rev(seq(0, 100, 10))^1.5)
  rownames(trait) <- letters[1:11]
  env <- cbind(e1 = runif(40, 0, 100), e2 = runif(40, 0, 100))
  space <- cbind(runif(40), runif(40))
  comm <- outer(env[, 1], trait[, 1], function(e, t) exp(-(e - t)^2 / 800)) + 0.01
  colnames(comm) <- rownames(trait)
  list(comm = comm, trait = trait, env = env, space = space)
}

test_that("semipartials are the loss of R2 of separate weighted regressions", {
  d <- partition_data()
  out <- trait_partition(d$comm, d$trait, d$env, d$space, n_perm = 3)
  w <- rowSums(d$comm)
  y <- trait_response(d$comm, d$trait, "mean")[, 1]
  x <- as.data.frame(cbind(d$env, d$space))
  names(x) <- c("e1", "e2", "s1", "s2")
  explained <- function(fit) sum(w * (fitted(fit) - weighted.mean(y, w))^2)
  full <- lm(y ~ ., x, weights = w)
  loss <- explained(full) - explained(lm(y ~ . - e1, x, weights = w))
  slope <- coef(lm(y ~ e1, x, weights = w))[["e1"]]
  expect_equal(out$semipartial["e1", "t1"], sign(slope) * sqrt(loss / sum(d$comm)))
  expect_equal(dimnames(out$semipartial), list(c("e1", "e2", "space1", "space2"), c("t1", "t2")))
  expect_equal(out$p_partial_max, pmax(out$p_partial_site, out$p_partial_species))
})

test_that("the sign of a semipartial follows the predictor", {
  d <- partition_data()
  set.seed(1)
  one <- trait_partition(d$comm, d$trait, d$env, d$space, n_perm = 1)$semipartial
  set.seed(1)
  other <- trait_partition(d$comm, d$trait, -d$env, d$space, n_perm = 1)$semipartial
  expect_equal(one[1:2, ], -other[1:2, ])
})

test_that("collinear predictors give no semipartials", {
  d <- partition_data()
  out <- trait_partition(d$comm, d$trait, d$env, cbind(d$env[, 1], d$space), n_perm = 2)
  expect_true(all(is.na(out$semipartial)))
  expect_true(all(is.na(out$p_partial_max)))
})

test_that("the statistic of a and c is the share of the component left by the other set", {
  d <- partition_data()
  out <- trait_partition(d$comm, d$trait, d$env, d$space, n_perm = 2)
  f <- out$fractions
  rownames(f) <- f$fraction
  expect_equal(f["a", "statistic"],
               (f["abc", "r2_component"] - f["bc", "r2_component"]) / (1 - f["bc", "r2_component"]))
  expect_equal(f["c", "statistic"],
               (f["abc", "r2_component"] - f["ab", "r2_component"]) / (1 - f["ab", "r2_component"]))
  expect_equal(f[c("abc", "ab", "bc", "b"), "statistic"], f[c("abc", "ab", "bc", "b"), "r2"])
  expect_true(is.na(f["d", "statistic"]))
  # one environmental variable: its semipartial test is the test of a
  one <- trait_partition(d$comm, d$trait[, 1], d$env[, 1], d$space, n_perm = 2)
  expect_equal(one$p_partial_site[1, 1], one$fractions$p_site[one$fractions$fraction == "a"])
  expect_equal(dimnames(one$semipartial), list(c("env1", "space1", "space2"), "trait1"))
})

test_that("a fraction that cannot be estimated has no p-value", {
  a <- example_a()
  set.seed(1)
  saturated <- matrix(rnorm(90), 10)
  out <- trait_partition(a$comm, a$trait, a$env, saturated, n_perm = 3)
  f <- out$fractions
  rownames(f) <- f$fraction
  expect_equal(f["bc", "r2_component"], 1)
  expect_true(all(is.na(f["a", c("statistic", "p_site", "p_species", "p_max")])))
  expect_false(anyNA(f["c", c("p_site", "p_species", "p_max")]))
  expect_true(is.na(permutation_p(c(0.1, 0.2), NA)))
  expect_equal(permutation_p(c(0.1, NA, 0.3), 0.2), 2 / 3)
})

shared_data <- function(structured = TRUE, side = 6, k = 12) {
  coords <- as.matrix(expand.grid(x = seq_len(side), y = seq_len(side)))
  mems <- NULL
  capture.output(mems <- as.matrix(adespatial::dbmem(coords, silent = TRUE)))
  space <- mems[, 1:3]
  env <- if (structured) 50 + 25 * mems[, 1] / max(abs(mems[, 1])) + rnorm(side^2, 0, 3) else {
    runif(side^2, 0, 100)
  }
  trait <- setNames(seq(0, 100, length.out = k), paste0("sp", seq_len(k)))
  comm <- outer(env, trait, function(e, t) exp(-(e - t)^2 / 800)) + 0.01
  colnames(comm) <- names(trait)
  list(comm = comm, trait = trait, env = cbind(e = env), space = space, coords = coords)
}

test_that("the shared fraction is tested with Moran spectral randomisation", {
  skip_if_not_installed("adespatial")
  d <- shared_data()
  set.seed(1)
  out <- trait_partition(d$comm, d$trait, d$env, d$space, n_perm = 19, shared_null = "msr",
                         coords = d$coords)
  f <- out$fractions
  rownames(f) <- f$fraction
  expect_equal(out$shared[c("null", "draws", "tested")], list(null = "msr", draws = 19L, tested = TRUE))
  expect_lte(out$shared$p_env_space, 0.05)
  expect_false(anyNA(f["b", c("p_site", "p_species", "p_max")]))
  expect_equal(f["b", "p_max"], max(f["b", c("p_site", "p_species")]))
  expect_true(is.na(f["d", "p_max"]))
  expect_output(print(out), "19 Moran spectral randomisations of env")
  # the other fractions do not depend on the request
  set.seed(1)
  plain <- trait_partition(d$comm, d$trait, d$env, d$space, n_perm = 19)
  expect_equal(plain$fractions[c("r2", "statistic")], f[c("r2", "statistic")], ignore_attr = TRUE)
  expect_true(all(is.na(plain$fractions$p_max[plain$fractions$fraction == "b"])))
  expect_equal(plain$shared, list(null = "none", draws = 0L, p_env_space = NA_real_, tested = FALSE))
  expect_output(print(plain), "permutations:")
  expect_false(any(grepl("shared", capture.output(print(plain)))))
})

test_that("the shared fraction stays untested when env is unrelated to space", {
  skip_if_not_installed("adespatial")
  d <- shared_data(structured = FALSE)
  set.seed(2)
  out <- trait_partition(d$comm, d$trait, d$env, d$space, n_perm = 19, shared_null = "msr",
                         coords = d$coords)
  expect_false(out$shared$tested)
  expect_gt(out$shared$p_env_space, 0.05)
  expect_true(all(is.na(out$fractions[out$fractions$fraction == "b", c("p_site", "p_species")])))
  expect_output(print(out), "not tested, env is unrelated to space")
  # a lenient level tests it
  set.seed(2)
  lenient <- trait_partition(d$comm, d$trait, d$env, d$space, n_perm = 19, shared_null = "msr",
                             coords = d$coords, level = 0.999)
  expect_true(lenient$shared$tested)
})

test_that("the shared fraction is tested with moves of env on a grid", {
  skip_if_not_installed("adespatial")
  d <- shared_data()
  moves <- torus_permutations(6)[1:19, ]
  set.seed(3)
  out <- trait_partition(d$comm, d$trait, d$env, d$space, n_perm = 19, site_permutations = moves,
                         shared_null = "sites")
  expect_equal(out$shared$null, "sites")
  expect_equal(out$shared$draws, 19L)
  expect_true(out$shared$tested)
  expect_false(anyNA(out$fractions$p_max[out$fractions$fraction == "b"]))
  expect_output(print(out), "19 moves of env alone")
  # empty sites stay in the landscape that env moves over
  empty <- d$comm
  empty[c(1, 8), ] <- 0
  with_empty <- trait_partition(empty, d$trait, d$env, d$space, n_perm = 19,
                                site_permutations = moves, shared_null = "sites")
  expect_true(with_empty$shared$tested)
})

test_that("the shared test explains what it needs", {
  skip_if_not_installed("adespatial")
  d <- shared_data()
  expect_error(trait_partition(d$comm, d$trait, d$env, d$space, n_perm = 2, shared_null = "sites"),
               "needs `site_permutations`")
  expect_error(trait_partition(d$comm, d$trait, d$env, d$space, n_perm = 2, shared_null = "msr"),
               "needs `coords`")
  expect_error(trait_partition(d$comm, d$trait, d$env, d$space, n_perm = 2, shared_null = "msr",
                               coords = d$coords[-1, ]), "one numeric row per site \\(36\\)")
  expect_error(trait_partition(d$comm, d$trait, d$env, d$space, n_perm = 2, level = 1), "level")
  local_mocked_bindings(has_adespatial = function() FALSE)
  expect_error(trait_partition(d$comm, d$trait, d$env, d$space, n_perm = 2, shared_null = "msr",
                               coords = d$coords), "needs the package adespatial")
})
