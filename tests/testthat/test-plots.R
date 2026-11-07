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
  pdf(NULL)
  on.exit(dev.off())
  out <- daphniid_result()
  marks <- plot_node_test(out, daphniids$tree, test = "site", level = 0.2)
  expect_gt(nrow(marks), 0)
  expect_silent(plot_node_test(out, daphniids$tree, test = "site", type = "fan", level = 0.2,
                               legend = FALSE))
  expect_equal(nrow(plot_node_test(out, daphniids$tree, level = 0)), 0)
})

test_that("symbols of a node stay apart, with predictors above events", {
  expect_equal(symbol_steps(c(1, -1, 1), 1), c(0.6, -0.6, 1.6) * 1.3)
})

test_that("labels go towards the root and stay readable", {
  phylogram <- label_place(5, 1, width = 2, gap = 1, angle = 0, room = 5)
  expect_equal(unlist(phylogram), c(x = 3, y = 1, srt = 0))
  expect_equal(label_place(5, 1, 2, 1, 0, room = 2)$x, 7)
  left <- label_place(-5, 0, width = 2, gap = 1, angle = pi, room = 5)
  expect_equal(c(left$x, left$srt), c(-3, 360))
  top <- label_place(0, 5, width = 2, gap = 1, angle = pi / 2, room = 5)
  expect_equal(c(top$x, top$y, top$srt), c(0, 3, 90))
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
  pdf(NULL)
  on.exit(dev.off())
  set.seed(1)
  trait <- stats::setNames(seq(0, 100, 10), letters[1:11])
  env <- runif(40, 0, 100)
  space <- cbind(runif(40), runif(40))
  comm <- outer(env, trait, function(e, t) exp(-(e - t)^2 / 800))
  colnames(comm) <- names(trait)
  out <- trait_partition(comm, trait, env, space, n_perm = 9)
  shown <- plot_trait_partition(out)
  expect_equal(shown$fraction, c("a", "b", "c", "d"))
  expect_equal(shown$r2_adjusted, out$fractions$r2_adjusted[match(c("a", "b", "c", "d"),
                                                                  out$fractions$fraction)])
  expect_true(any(shown$r2_adjusted < 0))
  expect_equal(shown$shown, pmax(shown$r2_adjusted, 0))
  expect_silent(plot_trait_partition(out, value = "r2", test = "site"))
})

test_that("the heatmap holds predictors and events, with blanks where there is no value", {
  pdf(NULL)
  on.exit(dev.off())
  out <- daphniid_result()
  values <- plot_node_heatmap(out, test = "site", level = 0.2)
  expect_equal(colnames(values), c(colnames(daphniids$env), colnames(daphniids$biogeography)))
  expect_equal(rownames(values), as.character(out$nodes$node))
  expect_true(anyNA(values))
  only_env <- daphniid_result(NULL)
  expect_equal(ncol(plot_node_heatmap(only_env, level = 0)), ncol(daphniids$env))
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
  pdf(NULL)
  on.exit(dev.off())
  x <- structure(list(component = "mean", fractions = data.frame(
    fraction = c("a", "b", "c", "d"), r2_adjusted = c(0.5, 0.1, 0.3, 0.1),
    p_max = c(0.01, NA, 0.2, NA))), class = "mc_trait_partition")
  expect_equal(plot_trait_partition(x)$shown, c(0.5, 0.1, 0.3, 0.1))
})
