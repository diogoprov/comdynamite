test_that("tree_layout() matches the coordinates of ape::plot.phylo()", {
  set.seed(3)
  tree <- ape::rtree(15)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  ape::plot.phylo(tree, plot = FALSE)
  lp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  l <- tree_layout(tree)
  expect_equal(l$x, lp$xx)
  expect_equal(l$y, lp$yy)
})

test_that("plot_rates() draws one L-shaped segment per branch and checks its input", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  tree <- ape::rtree(10)
  rates <- rexp(nrow(tree$edge))
  p <- plot_rates(tree, rates, tip_labels = TRUE, title = "x")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 2 * nrow(tree$edge))
  expect_equal(p$data$rate, rep(rates, 2))
  x <- trait_rates(tree, "iid", rate_sd = 0.5)
  expect_equal(plot_rates(tree, x)$data$rate, rep(x$edge_rates, 2))
  expect_s3_class(plot_rates(tree, rates - 1, log = FALSE), "ggplot")
  expect_error(plot_rates(tree, rates - 1), "positive")
  expect_error(plot_rates(tree, rates[-1]), "one value per branch")
  no_length <- tree
  no_length$edge.length <- NULL
  expect_error(plot_rates(no_length, rates), "branch lengths")
})

test_that("pcps_test() keeps what plot() needs and plot() checks its arguments", {
  skip_if_not_installed("ggplot2")
  set.seed(2)
  tree <- ape::rtree(10)
  comm <- matrix(rpois(20 * 10, 2), 20, dimnames = list(NULL, tree$tip.label))
  env <- cbind(temp = seq_len(20), rain = rnorm(20))
  res <- pcps_test(comm, tree, env, n_perm = 9)
  expect_equal(dim(res$species), c(10, ncol(res$scores)))
  expect_equal(length(res$share), ncol(res$scores))
  expect_equal(colnames(res$env), c("temp", "rain"))
  p <- plot(res, env = "rain")
  expect_s3_class(p, "ggplot")
  expect_equal(p$data$env, env[, "rain"])
  expect_s3_class(plot(res, env = NULL, species = FALSE), "ggplot")
  expect_s3_class(plot(res, axes = c(2, 3), species_labels = FALSE), "ggplot")
  expect_error(plot(res, axes = c(1, 1)), "two different axes")
  expect_error(plot(res, env = "pH"), "one column")
  old <- res
  old$species <- NULL
  expect_error(plot(old), "older version")
})

test_that("plot_rates() takes a result of trait_rates() alone", {
  skip_if_not_installed("ggplot2")
  set.seed(4)
  x <- trait_rates(ape::rtree(8), "gbm", rate_sd = 0.5)
  expect_equal(plot_rates(x)$data$rate, rep(x$edge_rates, 2))
  expect_error(plot_rates(x$tree), "`rates` is missing")
})

test_that("plot() of pcps_test() names an unnamed predictor", {
  skip_if_not_installed("ggplot2")
  set.seed(6)
  tree <- ape::rtree(10)
  comm <- matrix(rpois(200, 2), 20, dimnames = list(NULL, tree$tip.label))
  p <- plot(pcps_test(comm, tree, seq_len(20), n_perm = 9))
  expect_equal(p$scales$get_scales("colour")$name, "env1")
})

test_that("plot() of pcps_test() groups species into clades", {
  skip_if_not_installed("ggplot2")
  set.seed(8)
  tree <- ape::rcoal(12)
  comm <- matrix(rpois(30 * 12, 2), 30, dimnames = list(NULL, tree$tip.label))
  env <- cbind(temp = seq_len(30))
  res <- pcps_test(comm, tree, env, n_perm = 9)
  expect_s3_class(res$tree, "phylo")
  # k deepest clades: the first split is the root, and each further split is the crown nearest
  # the root
  two <- pcps_clades(2, res$tree, rownames(res$species))
  root_children <- res$tree$edge[res$tree$edge[, 1] == ape::Ntip(res$tree) + 1, 2]
  first <- if (root_children[1] <= 12) res$tree$tip.label[root_children[1]] else
    ape::extract.clade(res$tree, root_children[1])$tip.label
  expect_equal(nlevels(two), 2)
  expect_length(unique(two[rownames(res$species) %in% first]), 1)
  four <- pcps_clades(4, res$tree, rownames(res$species))
  expect_equal(nlevels(four), 4)
  expect_false(anyNA(four))
  # each clade of four is nested in one clade of two
  expect_true(all(tapply(as.character(two), four, function(z) length(unique(z))) == 1))
  p <- plot(res, clades = 4)
  expect_s3_class(p, "ggplot")
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomPolygon"), logical(1))))
  expect_equal(p$scales$get_scales("fill")$name, "temp")
  families <- stats::setNames(rep(c("Hylidae", "Leptodactylidae"), each = 5),
                              res$tree$tip.label[1:10])
  q <- plot(res, clades = families)
  expect_s3_class(q, "ggplot")
  expect_true(any(vapply(q$layers, function(l) {
    inherits(l$geom, "GeomPoint") && identical(l$aes_params$colour, "grey60")
  }, logical(1))))
  expect_s3_class(plot(res, clades = 3, env = NULL, species_labels = FALSE), "ggplot")
  expect_error(plot(res, clades = 1), "between 2 and")
  expect_error(plot(res, clades = c("a", "b")), "named by species")
  expect_error(plot(res, clades = c(x = "a")), "no species")
  old <- res
  old$tree <- NULL
  expect_error(plot(old, clades = 2), "older version")
})
