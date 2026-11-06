tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
comm <- matrix(c(1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1), 4,
               dimnames = list(NULL, c("a", "b", "c", "d")))

test_that("the order of the species does not matter", {
  shuffled <- comm[, c("d", "b", "a", "c")]
  expect_equal(phylo_beta(shuffled, tree), phylo_beta(comm, tree))
  trait <- c(a = 1, b = 2, c = 3, d = 5)
  env <- c(1, 4, 2, 7)
  expect_equal(trait_partition(shuffled, rev(trait), env, c(2, 1, 1, 3), n_perm = 1)$fractions$r2,
               trait_partition(comm, trait, env, c(2, 1, 1, 3), n_perm = 1)$fractions$r2)
})

test_that("species in only one source are named and dropped", {
  extra <- cbind(comm, e = c(1, 0, 0, 1))
  expect_warning(out <- phylo_beta(extra, tree), "1 species of `comm` not in the phylogeny dropped: e")
  expect_equal(out, phylo_beta(comm, tree))
  expect_warning(pruned <- phylo_beta(comm[, -4], tree), "not in `comm` dropped: d")
  expect_equal(pruned, phylo_beta(comm[, -4], ape::drop.tip(tree, "d")))
  expect_warning(trait_partition(comm, c(a = 1, b = 2, c = 3), 1:4, c(1, 3, 2, 2), n_perm = 1),
                 "`comm` not in `traits`")
  expect_warning(pcc_test(comm[, -4], tree, 1:4, n_perm = 2), "not in `comm` dropped: d")
  expect_warning(phylo_vector(comm[, -4], tree, 1:4, n_perm = 2), "not in `comm` dropped: d")
})

test_that("species without individuals raise no warning", {
  extra <- cbind(comm, e = 0)
  expect_no_warning(phylo_beta(extra, tree))
})

test_that("permutations without variance do not count", {
  expect_equal(permutation_p(c(1, NaN, 3), 2), 2 / 3)
})

test_that("long lists are cut", {
  expect_identical(species_list(letters[1:3]), "a, b, c")
  expect_match(species_list(letters), "and 16 more$")
})

test_that("names are required and must not repeat", {
  expect_error(phylo_beta(unname(comm), tree), "species names in its columns")
  expect_error(trait_partition(comm, c(1, 2, 3, 4), 1:4, c(1, 3, 2, 2), n_perm = 1),
               "needs species names")
  repeated <- comm
  colnames(repeated)[2] <- "a"
  expect_error(phylo_beta(repeated, tree), "repeated species: a")
  expect_error(suppressWarnings(phylo_beta(comm[, c("a", "b")] * 0 + 1, ape::rcoal(3))), "fewer than two")
  expect_error(suppressWarnings(phylo_beta(cbind(comm[, 1:2], z = 1)[, c(1, 3)] * c(0, 1, 1, 1), tree)),
               "fewer than two")
})

test_that("dropping species can empty a site", {
  lonely <- cbind(comm, e = 0)
  lonely[1, ] <- c(0, 0, 0, 0, 1)
  expect_error(suppressWarnings(phylo_beta(lonely, tree)), "at least one individual")
})
