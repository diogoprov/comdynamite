test_that("check_tree() and check_count() give messages that name the argument", {
  expect_error(check_tree("x"), "`tree` must be a `phylo` object")
  tree <- ape::rtree(5)
  expect_identical(check_tree(tree), tree)
  tree$edge.length <- NULL
  expect_error(check_tree(tree), "needs branch lengths")
  expect_silent(check_tree(tree, branch_lengths = FALSE))
  expect_error(check_count("a", "n"), "`n` must be one whole number, 1 or more")
  expect_error(check_count(-1, "nperm", min = 0), "0 or more")
  expect_silent(check_count(0, "nperm", min = 0))
})

test_that("functions with a tree stop with a readable message on a bad tree", {
  expect_error(trait_rates("x"), "`tree` must be a `phylo` object")
  expect_error(trait_mk("x", q = 1), "`tree` must be a `phylo` object")
  expect_error(traits_ou("x"), "`tree` must be a `phylo` object")
  expect_error(paint_regimes("x", 14), "`tree` must be a `phylo` object")
  expect_error(traits_ou(ape::rtree(5), n_traits = 0), "`n_traits`")
  expect_error(species_traits("x"), "number of species")
  expect_error(competition_matrix(5, "phylo"), "needs a 'phylo' tree")
})

test_that("other arguments are checked with their names", {
  set.seed(1)
  traits <- matrix(rnorm(10), 10, dimnames = list(paste0("sp", 1:10), NULL))
  expect_error(assemble_community(traits, size = "a"), "`size` must be one number")
  expect_error(assemble_community(traits, size = 11), "between 1 and the 10 species")
  ind <- data.frame(site = rep(c("a", "b"), each = 6), species = rep(c("x", "y"), 6),
                    trait = rnorm(12))
  expect_error(tstats(ind, nperm = -1), "`nperm`")
  expect_error(community_index(ind, stats::var, nperm = 1.5), "`nperm`")
  tree <- ape::rtree(4)
  comm <- matrix(1, 1, 2, dimnames = list(NULL, c("p", "q")))
  expect_error(damocles_loglik(comm, tree, 1, 1), "no species of `comm`")
})
