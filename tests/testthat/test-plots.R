skip_if_not_installed("ggplot2")

daphniid_result <- function(biogeography = daphniids$biogeography) {
  set.seed(1)
  node_test(daphniids$comm, node_contrasts(daphniids$tree), daphniids$env, biogeography, n_perm = 9)
}

test_that("each significant test gets a symbol, with the correction per set", {
  out <- daphniid_result()
  marks <- node_marks(out, "site", 0.2)
  p <- out$p_partial_site
  expect_equal(sum(marks$set == "environment"), sum(p <= 0.2, na.rm = TRUE))
  expect_equal(sum(marks$set == "biogeography"),
               sum(out$biogeography$p_partial_site <= 0.2, na.rm = TRUE))
  expected <- sum(p.adjust(p, "BH") <= 0.2, na.rm = TRUE)
  expect_equal(sum(marks$corrected[marks$set == "environment"]), expected)
  expect_true(all(c("node", "variable", "p", "p_adjusted") %in% names(marks)))
})

test_that("the three tests read their own p-values", {
  out <- daphniid_result(NULL)
  expect_equal(nrow(node_marks(out, "max", 0.2)), sum(out$p_partial_max <= 0.2, na.rm = TRUE))
  expect_equal(nrow(node_marks(out, "joint", 0.2)), sum(out$p_partial_joint <= 0.2, na.rm = TRUE))
})

test_that("nodes are drawn in both layouts, with or without symbols", {
  out <- daphniid_result()
  p <- expect_ggplot(plot_node_test(out, daphniids$tree, test = "site", level = 0.2))
  marks <- attr(p, "marks")
  expect_gt(nrow(marks), 0)
  # one symbol per node, set and level, predictors above the node and events below
  points <- ggplot2::layer_data(p, 4)
  expect_equal(nrow(points), nrow(unique(marks[c("node", "set", "corrected")])))
  expect_ggplot(plot_node_test(out, daphniids$tree, test = "site", type = "fan", level = 0.2,
                               legend = FALSE))
  expect_equal(nrow(attr(plot_node_test(out, daphniids$tree, level = 0), "marks")), 0)
  expect_warning(plot_node_test(out, daphniids$tree, col = c("red", "blue"), level = 0.2),
                 "`col` is deprecated")
  expect_warning(plot_node_test(out, daphniids$tree, cex = 1, level = 0.2), "`cex` is deprecated")
})

test_that("a tree that cannot hold the nodes is an error", {
  out <- daphniid_result(NULL)
  expect_error(plot_node_test(out, ape::rcoal(5)), "not the tree")
  caterpillar <- ape::stree(13, "left")
  caterpillar$tip.label <- daphniids$tree$tip.label
  out$nodes$n_species[1] <- 14
  expect_error(plot_node_test(out, caterpillar), "not the tree")
})

test_that("the partition plot returns the three fractions in order", {
  set.seed(1)
  trait <- stats::setNames(seq(0, 100, 10), letters[1:11])
  env <- runif(40, 0, 100)
  space <- cbind(runif(40), runif(40))
  comm <- outer(env, trait, function(e, t) exp(-(e - t)^2 / 800))
  colnames(comm) <- names(trait)
  out <- trait_partition(comm, trait, env, space, n_perm = 9)
  shown <- attr(expect_ggplot(plot_trait_partition(out)), "fractions")
  expect_equal(shown$fraction, c("a", "b", "c", "d"))
  expect_equal(shown$r2_adjusted, out$fractions$r2_adjusted[match(c("a", "b", "c", "d"),
                                                                  out$fractions$fraction)])
  expect_true(any(shown$r2_adjusted < 0))
  expect_equal(shown$shown, pmax(shown$r2_adjusted, 0))
  expect_ggplot(plot_trait_partition(out, value = "r2", test = "site"))
  expect_error(plot_trait_partition(list()), "trait_partition")
  expect_error(plot_trait_partition(out, labels = "a"), "two names")
  expect_warning(plot_trait_partition(out, cex = 1), "`cex` is deprecated")
  expect_warning(plot_trait_partition(out, test = "taxa"), "deprecated")
})

test_that("the heatmap holds predictors and events, with blanks where there is no value", {
  out <- daphniid_result()
  values <- attr(expect_ggplot(plot_node_heatmap(out, test = "site", level = 0.2)), "values")
  expect_equal(colnames(values), c(colnames(daphniids$env), colnames(daphniids$biogeography)))
  expect_equal(rownames(values), as.character(out$nodes$node))
  expect_true(anyNA(values))
  only_env <- daphniid_result(NULL)
  expect_equal(ncol(attr(plot_node_heatmap(only_env, level = 0), "values")),
               ncol(daphniids$env))
  expect_error(plot_node_heatmap(list()), "node_test")
  expect_warning(plot_node_heatmap(out, level = 0.2, cex = 1), "`cex` is deprecated")
})

test_that("labels of small fractions stay apart and inside the line", {
  edges <- c(0, 0.01, 0.02, 0.03, 1)
  at <- fraction_positions(edges, 0.1)
  expect_true(all(diff(at) >= 0.1 - 1e-9))
  crowded <- fraction_positions(c(0, 0.97, 0.98, 0.99, 1), 0.1)
  expect_true(all(diff(crowded) >= 0.1 - 1e-9))
  expect_lte(crowded[4] + 0.05, 1)
})

test_that("the second label goes to the side with more room", {
  x <- structure(list(component = "mean", fractions = data.frame(
    fraction = c("a", "b", "c", "d"), r2_adjusted = c(0.5, 0.1, 0.3, 0.1),
    p_max = c(0.01, NA, 0.2, NA))), class = "mc_trait_partition")
  p <- plot_trait_partition(x)
  expect_equal(attr(p, "fractions")$shown, c(0.5, 0.1, 0.3, 0.1))
  # c (0.3) on the right of the bar of space has more room than a (0.5) on its left? no
  expect_equal(ggplot2::layer_data(p, 2)$hjust[2], 1)
  # a tested shared fraction is coloured like the others
  x$fractions$p_max[2] <- 0.01
  expect_ggplot(plot_trait_partition(x))
})

test_that("plot_phylo_beta() orders sites and species", {
  set.seed(13)
  tree <- ape::rcoal(8)
  env <- seq(0, 100, length.out = 30)
  traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100), rmax_range = c(5, 30))
  comm <- simulate_static(traits, env, incidence = TRUE)$abundance
  keep <- rowSums(comm) > 0
  out <- phylo_beta_test(comm[keep, ], tree, env[keep], n_perm = 3)
  shown <- attr(expect_ggplot(plot_phylo_beta(out, comm[keep, ], tree)), "ordered")
  expect_equal(dim(shown), c(sum(keep), sum(colSums(comm) > 0)))
  expect_equal(sort(rowSums(shown)), sort(rowSums(comm[keep, ])))
  expect_ggplot(plot(out, comm[keep, ], tree, axis = 1))
  # the columns follow the tips of the ladderized tree
  ladder <- ape::ladderize(tree)
  tips <- ladder$edge[ladder$edge[, 2] <= 8, 2]
  expect_equal(colnames(shown), intersect(ladder$tip.label[tips], colnames(shown)))
  expect_error(plot_phylo_beta(out, comm[keep, ][-1, ], tree), "occupied sites")
  expect_error(plot_phylo_beta(out, comm[keep, ], tree, axis = 99), "`axis`")
  expect_warning(plot_phylo_beta(out, comm[keep, ], tree, cex = 1), "`cex` is deprecated")
  no_lengths <- tree
  no_lengths$edge.length <- NULL
  expect_ggplot(plot_phylo_beta(out, comm[keep, ], no_lengths))
})
