test_that("sim_shift_tree() returns a binary ultrametric tree with one regime per edge", {
  set.seed(1)
  sim <- sim_shift_tree(60, n_shifts = 2, ratio = 3)
  tree <- sim$tree
  expect_s3_class(tree, "phylo")
  expect_equal(ape::Ntip(tree), 60)
  expect_true(ape::is.binary(tree))
  expect_true(ape::is.ultrametric(tree, tol = 1e-8))
  expect_equal(max(ape::node.depth.edgelength(tree)), 1)
  expect_length(sim$regimes, nrow(tree$edge))
  expect_setequal(levels(sim$regimes), c("0", "1", "2"))
  expect_equal(sim$rates[1], 1)
  # The terminal edge of each tip carries the regime of the tip.
  terminal <- sim$regimes[match(seq_len(60), tree$edge[, 2])]
  expect_equal(as.integer(as.character(terminal)), sim$tip_regime)
  # Every regime is at least 5% of the tips.
  expect_true(all(table(sim$tip_regime) >= 3))
})

test_that("each shifted regime enters the tree through one edge", {
  set.seed(2)
  sim <- sim_shift_tree(80, n_shifts = 3, ratio = 2)
  tree <- sim$tree
  parent_edge <- match(tree$edge[, 1], tree$edge[, 2])
  # An edge enters a regime when its parent edge is in another regime, or when it leaves the
  # root in a shifted regime (the shift happened on a lineage that became a root edge).
  above <- ifelse(is.na(parent_edge), "0", as.character(sim$regimes[parent_edge]))
  entries <- as.character(sim$regimes) != above
  expect_equal(sort(as.character(sim$regimes[entries])), c("1", "2", "3"))
})

test_that("sim_shift_tree() without shifts has a single regime and is reproducible", {
  set.seed(3)
  a <- sim_shift_tree(30, n_shifts = 0)
  set.seed(3)
  b <- sim_shift_tree(30, n_shifts = 0)
  expect_identical(a, b)
  expect_equal(levels(a$regimes), "0")
})

test_that("faster regimes have more speciation per unit of branch length", {
  set.seed(4)
  ratios <- vapply(1:20, function(i) {
    sim <- sim_shift_tree(80, n_shifts = 1, ratio = 4)
    tree <- sim$tree
    internal <- tree$edge[, 2] > ape::Ntip(tree)
    births <- tapply(internal, sim$regimes, sum)
    length_in <- tapply(tree$edge.length, sim$regimes, sum)
    (births[["1"]] / length_in[["1"]]) / (births[["0"]] / length_in[["0"]])
  }, numeric(1))
  expect_gt(stats::median(ratios), 2)
})
