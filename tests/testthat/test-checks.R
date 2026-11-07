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
  expect_equal(one$r2, other$r2)
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
               "`predictor` has missing values in 1 sites")
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
  expect_equal(one$p_taxa, two$p_taxa)
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
  expect_output(print(phylo_vector(d$comm, d$tree, d$env[, 1], n_perm = 2, alpha = 0.01)),
                "no vector retained")
})

test_that("nodes that were not tested are reported", {
  d <- small_data()
  nodes <- node_contrasts(d$tree)
  d$comm[, 1:2] <- 0
  expect_output(print(node_test(d$comm, nodes, d$env, n_perm = 2)), "nodes not tested")
})
