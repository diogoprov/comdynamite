small_data <- function() {
  set.seed(3)
  tree <- ape::rcoal(8)
  env <- cbind(a = runif(20), b = runif(20))
  comm <- matrix(rpois(160, 3) + 1, 20, dimnames = list(paste0("s", 1:20), tree$tip.label))
  rownames(env) <- rownames(comm)
  list(tree = tree, env = env, comm = comm)
}

test_that("sites are matched by name in any order", {
  d <- small_data()
  set.seed(1)
  one <- pcc_test(d$comm, d$tree, d$env, n_perm = 9)
  set.seed(1)
  other <- pcc_test(d$comm, d$tree, d$env[20:1, ], n_perm = 9)
  expect_equal(one$r2_mean, other$r2_mean)
  different <- d$env
  rownames(different)[1] <- "elsewhere"
  expect_error(pcc_test(d$comm, d$tree, different, n_perm = 9), "different names: elsewhere, s1")
  duplicated_names <- d$env
  rownames(duplicated_names)[2] <- "s1"
  expect_equal(check_sites(duplicated_names, d$comm, "x"), duplicated_names)
  expect_equal(check_sites(unname(d$env), d$comm, "x"), unname(d$env))
})

test_that("predictors with problems are named", {
  d <- small_data()
  env <- as.data.frame(d$env)
  with_id <- env
  with_id$place <- "x"
  expect_error(pcc_test(d$comm, d$tree, with_id, n_perm = 2), "not numeric: place")
  expect_error(pcc_test(d$comm, d$tree, letters[1:20], n_perm = 2), "`env` is not numeric")
  missing_value <- d$env
  missing_value[2, 1] <- NA
  expect_error(phylo_vector(d$comm, d$tree, missing_value[, 1], n_perm = 2),
               "`env` has missing values in 1 sites")
  constant <- d$env
  constant[, "b"] <- 1
  expect_error(node_test(d$comm, node_contrasts(d$tree), constant, n_perm = 2),
               "without variance in `env`: b")
  expect_error(pcc_test(d$comm, d$tree, cbind(d$env, 2 * d$env[, 1]), n_perm = 2), "collinear")
  expect_error(pcc_test(d$comm, d$tree, d$env[-1, ], n_perm = 2), "19 rows and `comm` has 20 sites")
  expect_error(phylo_vector(d$comm, d$tree, d$env, n_perm = 2), "single variable")
  expect_error(trait_partition(d$comm, setNames(1:8, d$tree$tip.label), NULL, d$env, n_perm = 2),
               "`env` is NULL")
})

test_that("communities with problems are explained", {
  d <- small_data()
  with_id <- data.frame(place = rownames(d$comm), d$comm)
  expect_error(phylo_beta(with_id, d$tree), "not numeric: place")
  missing_value <- d$comm
  missing_value[1, 1] <- NA
  expect_error(phylo_beta(missing_value, d$tree), "missing values")
  expect_error(phylo_beta(-d$comm, d$tree), "negative")
  empty <- d$comm
  empty[3, ] <- 0
  expect_error(phylo_beta(empty, d$tree), "have none: s3")
  only_outside <- d$comm
  only_outside[3, ] <- 0
  only_outside <- cbind(only_outside, stranger = 0)
  only_outside[3, "stranger"] <- 1
  expect_error(suppressWarnings(phylo_beta(only_outside, d$tree)),
               "1 have none of the species shared with the phylogeny: s3")
  expect_error(phylo_beta(d$comm, "((a,b),c);"), "`phylo` object")
})

test_that("permutations of sites are checked", {
  d <- small_data()
  expect_error(pcc_test(d$comm, d$tree, d$env, n_perm = 2, site_permutations = 1:20), "one column")
  repeated <- rbind(c(1, 1, 3:20))
  expect_error(pcc_test(d$comm, d$tree, d$env, n_perm = 2, site_permutations = repeated),
               "each site once")
})

test_that("traits with problems are explained", {
  d <- small_data()
  traits <- setNames(as.numeric(1:8), d$tree$tip.label)
  expect_error(trait_partition(d$comm, NULL, d$env, d$env, n_perm = 2), "`traits` is NULL")
  expect_error(trait_partition(d$comm, setNames(letters[1:8], d$tree$tip.label), d$env, d$env,
                               n_perm = 2), "not numeric")
  missing_value <- traits
  missing_value[2] <- NA
  expect_error(trait_partition(d$comm, missing_value, d$env, d$env, n_perm = 2), "missing values")
  expect_error(trait_partition(d$comm, traits * 0, d$env, d$env, n_perm = 2), "without variance")
})

test_that("species without individuals do not change the analyses", {
  d <- small_data()
  absent <- d$comm
  absent[, 1:2] <- 0
  pruned <- ape::drop.tip(d$tree, colnames(absent)[1:2])
  set.seed(2)
  one <- pcc_test(absent, d$tree, d$env, n_perm = 9)
  set.seed(2)
  two <- pcc_test(absent[, pruned$tip.label], pruned, d$env, n_perm = 9)
  expect_equal(one$r2_adjusted, two$r2_adjusted)
  expect_equal(one$p_species, two$p_species)
  # traits of absent species cannot reach present ones
  traits <- setNames(c(5, 5, 1, 1, 1, 1, 1, 2), colnames(absent))
  expect_s3_class(trait_partition(absent, traits, d$env[, 1], d$env[, 2], n_perm = 5),
                  "mc_trait_partition")
  absent[, 3:7] <- 0
  expect_error(pcc_test(absent, d$tree, d$env, n_perm = 2), "at least 3 species")
})

test_that("results keep and print the permutations", {
  d <- small_data()
  nodes <- node_contrasts(d$tree)
  out <- node_test(d$comm, nodes, d$env, n_perm = 9)
  expect_equal(out$permutations, list(species = 9, sites = 9, null = "free"))
  expect_output(print(out), "9 free of sites, 9 of species; smallest possible p: 0.1")
  perms <- t(replicate(4, sample(20)))
  restricted <- pcc_test(d$comm, d$tree, d$env, n_perm = 9, site_permutations = perms)
  expect_equal(restricted$permutations$null, "restricted")
  expect_output(print(restricted), "4 restricted of sites, 9 of species; smallest possible p: 0.2")
  expect_output(print(phylo_beta_test(d$comm, d$tree, d$env, n_perm = 2)), "permutations")
  expect_output(print(trait_partition(d$comm, setNames(as.numeric(1:8), d$tree$tip.label),
                                      d$env[, 1], d$env[, 2], n_perm = 2)), "permutations")
  expect_output(print(phylo_vector(d$comm, d$tree, d$env[, 1], n_perm = 2, level = 0.01)),
                "no vector retained")
})

test_that("nodes that were not tested are reported", {
  d <- small_data()
  nodes <- node_contrasts(d$tree)
  d$comm[, 1:2] <- 0
  expect_output(print(node_test(d$comm, nodes, d$env, n_perm = 2)), "nodes not tested")
})

test_that("`alpha` for a significance level is deprecated", {
  expect_warning(expect_equal(deprecated_level(0.05, 0.1), 0.1), "use `level`")
  expect_equal(deprecated_level(0.05, NULL), 0.05)
})

test_that("restricted permutations accept empty sites", {
  set.seed(4)
  tree <- ape::rcoal(6)
  comm <- matrix(rpois(96, 2), 16, dimnames = list(NULL, tree$tip.label))
  comm[c(1, 7), ] <- 0
  env <- cbind(a = rnorm(16))
  moves <- torus_permutations(4)[1:9, ]
  expect_error(pcc_test(comm, tree, env, n_perm = 3), "every site")
  expect_equal(pcc_test(comm, tree, env, n_perm = 3, site_permutations = moves)$permutations$sites, 9)
  expect_equal(phylo_beta_test(comm, tree, env, n_perm = 3,
                               site_permutations = moves)$permutations$sites, 9)
  traits <- cbind(t1 = rnorm(6))
  rownames(traits) <- tree$tip.label
  out <- trait_partition(comm, traits, env, cbind(s = rnorm(16)), n_perm = 3,
                         site_permutations = moves)
  expect_equal(out$permutations$null, "restricted")
  expect_error(trait_partition(comm, traits, env, cbind(s = rnorm(16)), n_perm = 3), "every site")
})

test_that("the simulator explains wrong input", {
  land <- landscape_implicit(3)
  traits <- species_traits(3)
  expect_error(simulate_metacommunity(land, traits, diag(2)), "one row and one column")
  expect_error(simulate_metacommunity(land, traits[, 1:2], diag(3)), "r_max")
  expect_error(mc_dynamics(matrix(c(2, 3), 1), diag(2) * 0.1, 0.1, n_steps = 2), "two sites")
  expect_silent(mc_dynamics(matrix(c(2, 3), 1), diag(2) * 0.1, 0, n_steps = 2))
  expect_error(logistic_gradients(5, n_sites = 6), "at least 8")
})

test_that("patch_permutations() explains repeated sites", {
  expect_error(patch_permutations(cbind(rep(1:3, each = 3), 0), cbind(rep(c(1, -1, 0), 3))),
               "share their coordinates")
})

test_that("phylo_vector() accepts empty sites and the old name of `env`", {
  set.seed(8)
  tree <- ape::rcoal(6)
  comm <- matrix(rpois(96, 2), 16, dimnames = list(NULL, tree$tip.label))
  comm[c(1, 7), ] <- 0
  env <- rnorm(16)
  moves <- torus_permutations(4)[1:9, ]
  expect_error(phylo_vector(comm, tree, env, n_perm = 3), "every site")
  out <- phylo_vector(comm, tree, env, n_perm = 3, n_sim = 3, site_permutations = moves)
  expect_equal(out$permutations$sites, 9)
  full <- comm + 1
  set.seed(9)
  one <- phylo_vector(full, tree, env, n_perm = 3, n_sim = 3)
  set.seed(9)
  expect_warning(other <- phylo_vector(full, tree, predictor = env, n_perm = 3, n_sim = 3),
                 "use `env`")
  expect_equal(one$p_max, other$p_max)
})
