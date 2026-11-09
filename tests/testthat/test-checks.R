test_that("check_tree() and check_count() give messages that name the argument", {
  expect_error(check_tree("x"), "`tree` must be a `phylo` object")
  tree <- ape::rtree(5)
  expect_identical(check_tree(tree), tree)
  tree$edge.length <- NULL
  expect_error(check_tree(tree), "needs branch lengths")
  expect_silent(check_tree(tree, branch_lengths = FALSE))
  expect_error(check_count("a", "n"), "`n` must be one whole number, 1 or more")
  expect_error(check_count(-1, "n_perm", min = 0), "0 or more")
  expect_silent(check_count(0, "n_perm", min = 0))
})

test_that("functions with a tree stop with a readable message on a bad tree", {
  expect_error(trait_rates("x"), "`tree` must be a `phylo` object")
  expect_error(trait_mk("x", q = 1), "`tree` must be a `phylo` object")
  expect_error(trait_ou("x"), "`tree` must be a `phylo` object")
  expect_error(paint_regimes("x", 14), "`tree` must be a `phylo` object")
  expect_error(trait_ou(ape::rtree(5), n_traits = 0), "`n_traits`")
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
  expect_error(tstats(ind, n_perm = -1), "`n_perm`")
  expect_error(community_index(ind, stats::var, n_perm = 1.5), "`n_perm`")
  tree <- ape::rtree(4)
  comm <- matrix(1, 1, 2, dimnames = list(NULL, c("p", "q")))
  expect_error(damocles_loglik(comm, tree, 1, 1), "no species of `comm`")
})

test_that("renamed functions and arguments still work, with a warning", {
  set.seed(5)
  tree <- ape::rtree(6)
  set.seed(1); new <- trait_ou(tree)
  set.seed(1); expect_warning(old <- traits_ou(tree), "deprecated")
  expect_equal(old, new)
  expect_warning(x <- trait_rates(tree, "iid", root_value = 3), "`root_value` is deprecated")
  expect_equal(x$edge_rates, x$edge_rate)
  expect_warning(h <- trait_mk(tree, q = 0.1, root = "1"), "`root` is deprecated")
  ind <- data.frame(site = rep(c("a", "b"), each = 6), species = rep(c("x", "y"), 6),
                    trait = rnorm(12))
  expect_warning(res <- tstats(ind, nperm = 9), "`nperm` is deprecated")
  expect_equal(dim(attr(res, "null"))[1], 9)
  expect_warning(community_index(ind, stats::var, nperm = 0), "`nperm` is deprecated")
  expect_warning(tstats_summary(tstats(ind, n_perm = 49), level = 0.05), "coverage")
  comm <- matrix(c(1, 0, 1, 0, 1, 0), 1, dimnames = list(NULL, tree$tip.label))
  expect_warning(fit_damocles(comm, tree), "deprecated")
})

test_that("plot() of trait_rarefaction() takes `colours` and warns on `colour`", {
  skip_if_not_installed("ggplot2")
  set.seed(2)
  ind <- data.frame(site = rep(c("a", "b"), each = 20), species = rep(c("x", "y"), 20),
                    trait = rnorm(40))
  r <- trait_rarefaction(ind, n = 2:5, reps = 3)
  expect_s3_class(plot(r, colours = "red"), "ggplot")
  expect_warning(plot(r, colour = "red"), "`colour` is deprecated")
test_that("plots and rarefaction name the bad argument", {
  tree <- ape::rtree(5)
  expect_error(plot_node_test(list(), tree), "`x` must come from node_test()")
  expect_error(plot_phylo_beta(list(), matrix(1), tree), "`x` must come from phylo_beta_test()")
  ind <- data.frame(site = rep(c("a", "b"), each = 6), species = rep(c("x", "y"), 6),
                    trait = rnorm(12))
  expect_error(trait_rarefaction(ind, n = "a"), "`n` must be whole numbers")
})
