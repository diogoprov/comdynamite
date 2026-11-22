# Internal helpers of trait simulation, landscapes, DAMOCLES and permutations, called directly
# so that each of their branches runs whatever the random draws of the other tests.

test_that("helpers of species_traits() check and place the traits", {
  expect_identical(species_labels(3, FALSE), c("sp1", "sp2", "sp3"))
  expect_error(species_labels(0, FALSE), "`species` must be")
  expect_error(species_labels(3, TRUE), "needs a tree")
  expect_error(check_breadth_moments(0, 0), "`sigma_mean` must be greater than 0")
  expect_true(check_breadth_moments(1, 0))
  given <- given_traits(NULL, c(1, 2), NULL, c("a", "b"))
  expect_named(given, c("optima", "sigma", "r_max"))
  expect_null(given$optima)
  expect_identical(given$sigma, c(1, 2))
  expect_error(given_traits(NULL, c(1, 0), NULL, c("a", "b")), "`sigma` must be positive")
  expect_error(given_traits(NULL, NULL, c(1, -1), c("a", "b")), "`r_max` cannot be negative")
  expect_identical(scale_optima(c(0, 1, 2), NULL, TRUE, FALSE, c(0, 100)), c(0, 50, 100))
  expect_identical(scale_optima(c(0, 1, 2), c(5, 6, 9), FALSE, FALSE, c(0, 100)), c(5, 6, 9))
  expect_identical(scale_optima(c(0, 1, 2), c(5, 6, 9), TRUE, FALSE, c(0, 100)), c(0, 25, 100))
  expect_identical(scale_optima(c(2, 0, 1), NULL, TRUE, TRUE, c(0, 10)), c(10, 0, 5))
  expect_identical(with_signal(c(3, 1, 2), NULL, FALSE), c(3, 1, 2))
  set.seed(1)
  breadth <- draw_breadth(5, 0, 1)
  expect_length(breadth, 5)
  expect_gt(min(breadth), 0)
  expect_identical(trait_column(matrix(1:2, 2, 1, dimnames = list(c("a", "b"), NULL)), "x"),
                   c(a = 1L, b = 2L))
  expect_null(names(trait_column(data.frame(v = c(1, 2)), "x")))
  expect_error(trait_column(matrix(1:4, 2), "x"), "`x` needs one column")
  expect_identical(check_niche_cor(diag(2), "`cor`", 2L), diag(2))
  expect_error(check_niche_cor(diag(3), "`cor`", 2L), "must be a 2 x 2 numeric matrix")
})

test_that("growth_rates_axes() multiplies the responses to each axis", {
  traits <- data.frame(optimum = c(1, 3), sigma = 1, r_max = 2, row.names = c("a", "b"))
  env <- cbind(x = c(1, 3), y = c(3, 1))
  r <- growth_rates_axes(env, list(traits, traits))
  expect_identical(r, growth_rates(env[, 1], traits) * growth_rates(env[, 2], traits))
  expect_error(growth_rates_axes(env, list(traits)), "`env` has 2 axes")
})

test_that("helpers of the DAMOCLES functions check the design and the rates", {
  skip_if_not_installed("ape")
  set.seed(1)
  tree <- ape::rcoal(6)
  n_edge <- nrow(tree$edge)
  at_root <- tree$edge[, 1] == 7
  msg <- "the regime differs"
  expect_identical(levels(design_regimes(tree, rep(c("a", "b"), 5), at_root, msg)), c("a", "b"))
  expect_error(design_regimes(tree, "a", at_root, msg), "`regimes` needs one regime")
  rooted <- tree
  rooted$root.edge <- 1
  expect_error(design_regimes(rooted, ifelse(at_root & !duplicated(at_root), "a", "b"),
                              at_root, msg), "the regime differs")
  expect_identical(design_covariate(tree, rep(2, n_edge), at_root, msg), rep(0, n_edge))
  expect_error(design_covariate(tree, rep(0, n_edge), at_root, msg), "`edge_covariate` needs")
  expect_error(design_covariate(rooted, seq_len(n_edge), at_root, msg), "the covariate differs")
  expect_identical(design_names(factor(c("a", "b")), 2L, 1L, FALSE, TRUE, TRUE),
                   c("gamma_a", "gamma_b", "mu", "beta_mu", "sympatry"))
  expect_identical(design_names(NULL, 1L, 1L, TRUE, FALSE, FALSE),
                   c("gamma", "mu", "beta_gamma"))
  carried <- damocles_carry(c(1, 0), c(0, 1), 0, 1, 1)
  expect_identical(carried, list(absent = c(1, 0), present = c(0, 1)))
  given <- damocles_given_rates(tree, list(1, seq_len(n_edge)), 2, 0.5, NULL, NULL, 2)
  expect_named(given, c("parameters", "edge_rates", "sympatry", "loglik"))
  expect_identical(given$parameters, data.frame(gamma = c(1, NA), mu = c(2, 2)))
  expect_identical(given$sympatry, c(0.5, 0.5))
  expect_identical(given$loglik, NA_real_)
  expect_length(given$edge_rates, 2)
  expect_error(damocles_given_rates(tree, 1, 2, 0, "a", NULL, 2), "are for estimated rates")
  expect_error(damocles_given_rates(tree, 1, 1:3, 0, NULL, NULL, 2), "one per community")
  expect_error(damocles_given_rates(tree, 1, 2, NULL, NULL, NULL, 2), "`sympatry = NULL`")
  set.seed(2)
  comm <- simulate_damocles(tree, 1, 1, n = 2)
  fitted <- damocles_fitted_rates(damocles_presence(comm, tree), tree, "both", TRUE, 0, NULL,
                                  NULL, "both")
  expect_named(fitted, c("parameters", "edge_rates", "sympatry", "loglik"))
  expect_named(fitted$parameters, c("gamma", "mu"))
  set.seed(3)
  site <- damocles_presence(comm, tree)[1, , drop = FALSE]
  row <- damocles_test_row(site, tree, tip_distances(tree)[colnames(site), colnames(site)],
                           damocles_edge_rates(tree, 1, 1), 9, "both", 0)
  expect_named(row, c("richness", "mpd", "ses_mpd_damocles", "p_mpd_damocles", "ses_mpd_random",
                      "p_mpd_random", "mntd", "ses_mntd_damocles", "p_mntd_damocles",
                      "ses_mntd_random", "p_mntd_random", "n_valid"))
  expect_identical(row[["richness"]], sum(site))
})

test_that("helpers of clade_regimes() find the node and the edges of a clade", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  expect_identical(clade_node(tree, 6, "A"), 6)
  expect_identical(clade_node(tree, "c", "A"), 3L)
  expect_identical(clade_node(tree, c("c", "d")), ape::getMRCA(tree, c("c", "d")))
  expect_error(clade_node(tree, 99, "A"), "clade A must be one node")
  expect_error(clade_node(tree, "z", "A"), "tips of clade A not in `tree`")
  node <- ape::getMRCA(tree, c("a", "b"))
  expect_identical(clade_edges(tree, c("a", "b"), "A", FALSE), tree$edge[, 1] == node)
  expect_identical(clade_edges(tree, node, "A", TRUE),
                   tree$edge[, 1] == node | tree$edge[, 2] == node)
  expect_identical(clade_edges(tree, "c", "A", FALSE), tree$edge[, 2] == 3)
  expect_warning(edges <- clade_edges(tree, c("a", "c"), "A", FALSE), "not monophyletic")
  expect_identical(edges, rep(TRUE, 6))
})

test_that("the parts of damocles_partition_test() give no test without draws", {
  expect_identical(partition_bootstrap(NULL, NULL, NULL, 0, "both", "both", 0, NULL, 1, 1),
                   list(null = numeric(0), p = NA_real_))
  expect_identical(partition_random(NULL, NULL, 0, 2, 0.25, NULL, identity, 1, 1),
                   list(null = numeric(0), p = NA_real_))
})

test_that("helpers of the landscape functions check their arguments", {
  expect_identical(check_axis_cor_matrix(diag(2), 2L), diag(2))
  expect_error(check_axis_cor_matrix(matrix(1, 2, 2), 2L), "positive definite")
  expect_null(matrix_coords(NULL, 3))
  expect_identical(matrix_coords(cbind(x = 1:2, y = 3:4), 2), data.frame(x = 1:2, y = 3:4))
  expect_error(matrix_coords(cbind(1:2, 3:4), 3), "`coords` must have two columns")
  expect_null(matrix_dist(NULL, 3))
  expect_identical(matrix_dist(data.frame(a = c(0, 1), b = c(1, 0)), 2),
                   matrix(c(0, 1, 1, 0), 2))
  expect_error(matrix_dist(diag(2), 3), "`dist` must be a numeric matrix of 3 x 3")
  expect_length(torus_maps(FALSE, TRUE), 1)
  expect_length(torus_maps(TRUE, FALSE), 4)
  expect_length(torus_maps(TRUE, TRUE), 8)
  images <- torus_images(torus_maps(FALSE, FALSE), c(2, 3))
  expect_identical(dim(images), c(6L, 6L))
  expect_identical(images[1, ], 1:6)
  expect_identical(pulse_sites(c(TRUE, FALSE, TRUE), 3), c(1L, 3L))
  expect_identical(pulse_sites(2, 3), 2)
  expect_error(pulse_sites(TRUE, 3), "a logical `sites`")
  expect_error(pulse_sites(c(1, 1), 3), "distinct indices")
  expect_true(check_pulse_timing(2, Inf, NULL, 5))
  expect_true(check_pulse_timing(2, 1, 3, 5))
  expect_error(check_pulse_timing(6, Inf, NULL, 5), "after the last step \\(5\\)")
  expect_error(check_pulse_timing(2, 3, 3, 5), "`every` must be longer")
})

test_that("helpers of the river functions build the network and its covariance", {
  children <- list(c(2L, 3L), integer(0), integer(0))
  expect_identical(river_upwards(children, 1L), c(1L, 2L, 3L))
  expect_error(river_upwards(list(integer(0), 3L, 2L), 1L), "must form a tree")
  size <- river_size(children, c(1L, 2L, 3L), c(0, 2, 3))
  expect_identical(size, list(magnitude = c(2, 1, 1), upstream = c(5, 2, 3),
                              order = c(2L, 1L, 1L)))
  expect_equal(river_angle(children, c(1L, 2L, 3L), size$magnitude, 1L),
               c(pi / 2, pi / 3, 2 * pi / 3))
  net <- river_table(c(NA, 1, 1), c(0, 2, 3))
  expect_identical(river_given_coords(net, net[, c("downstream", "length")]), net)
  moved <- river_given_coords(net, data.frame(x = c(1, 2, 3), y = c(0, 0, 1)))
  expect_identical(moved$x, c(1, 2, 3))
  expect_identical(moved$y, c(0, 0, 1))
  expect_error(river_given_coords(net, data.frame(x = c("a", "b", "c"), y = 1)),
               "must be finite numbers")
  expect_identical(river_site_env(NULL, net, c(0, 30)), c(0, 20, 30))
  expect_identical(river_site_env(c(3, 2, 1), net, c(0, 30)), c(3, 2, 1))
  expect_error(river_site_env(1, net, c(0, 30)), "one finite number per site of `network` \\(3\\)")
  channel <- river_channel(river_downstream(net), 1, Inf)
  expect_identical(diag(channel), c(0, 0, 0))
  expect_equal(rowSums(channel), c(1, 1, 1))
  expect_null(covariance_component(NULL, "tail_up"))
  expect_identical(covariance_component(c(range = 2, sill = 1), "tail_up"), c(1, 2))
  expect_identical(covariance_component(c(1, 2), "tail_up"), c(1, 2))
  expect_error(covariance_component(c(1, 0), "tail_up"), "`tail_up` must be NULL or two")
})

test_that("helpers of rigid_permutations() enumerate the moves and collect the layouts", {
  expect_error(default_gap(matrix(0, 6, 6)), "share their coordinates")
  expect_identical(default_gap(as.matrix(stats::dist(1:6))), 3)
  coords <- as.matrix(expand.grid(x = 1:4, y = 1:4))
  lattice <- lattice_of(coords, 1.5 * (1 - 1e-6))
  cells <- lattice$cells[1:2, , drop = FALSE]
  one <- symmetry_moves(1L, cells, lattice, 0)
  expect_type(one, "list")
  expect_true(all(lengths(one) == 3))
  by_symmetry <- lapply(seq_along(lattice$maps), symmetry_moves, cells, lattice, 0)
  expect_identical(do.call(rbind, unlist(by_symmetry, recursive = FALSE)),
                   lattice_moves(1:2, lattice, 0))
  # distinct layouts other than the observed one, up to `n_perm`
  event <- c(1, 0, 0)
  acc <- new.env(parent = emptyenv())
  acc$perms <- vector("list", 2)
  acc$keys <- character(2)
  acc$found <- 0
  tries <- list(NULL, 1:3, c(2L, 1L, 3L), c(2L, 1L, 3L), c(3L, 2L, 1L), c(2L, 3L, 1L))
  collect_layouts(seq_along(tries), function(pick) tries[[pick]], acc, event,
                  layout_key(event, 1:3), 2)
  expect_identical(acc$found, 2)
  expect_identical(acc$perms, list(c(2L, 1L, 3L), c(3L, 2L, 1L)))
  set.seed(1)
  event <- as.numeric(coords[, "x"] <= 2 & coords[, "y"] <= 2)
  distances <- as.matrix(stats::dist(coords))
  drawn <- event_placements(event, "e", coords, distances, lattice, 3, 1.5, 0, Inf, 60)
  expect_identical(dim(drawn), c(3L, 16L))
  expect_null(attr(drawn, "moves"))
  groups <- list(which(event == 1))
  moves <- lapply(groups, landscape_moves, coords, lattice, 0)
  placed <- place_groups(1L, groups, moves, coords, distances, lattice, 1.5, Inf)
  expect_identical(sort(placed), 1:16)
})

test_that("helpers of trait_ou() arrange the parameters by regime", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  expect_identical(ou_maps(tree, NULL, 0), list(c("0" = 1), c("0" = 1), c("0" = 1), c("0" = 2)))
  expect_identical(vapply(ou_maps(tree, c("x", "y", "y", "x"), 0), names, ""),
                   c("x", "y", "y", "x"))
  expect_error(ou_maps(tree, NULL, c(0, 1)), "several optima need `regimes`")
  expect_error(ou_maps(tree, c("x", "y"), 0), "one regime per branch \\(4\\)")
  expect_error(ou_maps(tree, tree, 0), "needs a mapped history")
  expect_identical(ou_by_regime(2, "theta", c("x", "y")), list(x = 2, y = 2))
  expect_identical(ou_by_regime(c(y = 1, x = 2), "theta", c("x", "y")), list(x = 2, y = 1))
  expect_error(ou_by_regime(c(1, 2), "theta", c("x", "y")), "`theta` must be named")
  expect_identical(ou_optima(matrix(c(1, 2), 1), c("x", "y")), list(x = c(1, 2), y = c(1, 2)))
  expect_identical(ou_optima(matrix(1:4, 2, dimnames = list(c("y", "x"), NULL)), c("x", "y")),
                   list(x = c(2L, 4L), y = c(1L, 3L)))
  expect_identical(ou_optima(3, "x"), list(x = 3))
  scalar <- ou_rates(2, c("x", "y"), 3)
  expect_identical(scalar, list(rates = list(x = 2, y = 2), multivariate = FALSE, n_traits = 3,
                                trait_names = NULL))
  several <- ou_rates(diag(2), "x", 1)
  expect_true(several$multivariate)
  expect_identical(several$n_traits, 2L)
  expect_identical(several$rates, list(x = diag(2)))
  expect_error(ou_rates(list(x = diag(2), y = 1), c("x", "y"), 1), "not both")
  expect_error(ou_rates(list(x = diag(2), y = diag(3)), c("x", "y"), 1), "the same size")
  expect_error(ou_rates(diag(2), "x", 3), "`n_traits` is 3")
  expect_error(check_scalar_rates(list(x = 1, y = -1)), "`sigma2` must be zero or more")
  expect_identical(check_scalar_rates(list(x = 1)), list(x = 1))
  expect_identical(ou_n_traits(0, 2, FALSE), 2)
  expect_identical(ou_n_traits(matrix(0, 1, 3), 1, FALSE), 3L)
  expect_error(ou_n_traits(matrix(0, 1, 3), 2, TRUE), "`theta` has 3 columns")
  # without noise the traits stay at the root value when it is the optimum
  x <- ou_simulate(tree, ou_maps(tree, NULL, 0), list("0" = 1), list("0" = 0), 1, 1, 2, FALSE)
  expect_identical(x, matrix(1, 5, 2))
})

test_that("helpers of rate_matrix() check the variances and correlations", {
  expect_identical(rate_matrix_size(c(1, 2), 0, NULL), 2L)
  expect_identical(rate_matrix_size(1, diag(3), NULL), 3L)
  expect_identical(rate_matrix_size(list(a = c(1, 2, 3)), 0, NULL), 3L)
  expect_identical(rate_matrix_size(1, 0, 4), 4L)
  expect_error(rate_matrix_size(numeric(0), 0, NULL), "at least one trait")
  expect_identical(check_regime_list(list(a = 1), "variance", "correlation"), list(a = 1))
  expect_error(check_regime_list(list(1), "variance", "correlation"),
               "for \"correlation\", `variance` must be a list named by regime")
  expect_identical(common_correlation(0.5, 2), matrix(c(1, 0.5, 0.5, 1), 2))
  expect_error(common_correlation(2, 2), "one number in \\[-1, 1\\]")
  expect_identical(check_rate_correlation(diag(2), 2), diag(2))
  expect_error(check_rate_correlation(diag(3), 2), "symmetric 2 x 2 matrix")
  expect_identical(constrained_rates("diagonal", c(1, 4), 0, 2L, NULL, NULL), diag(c(1, 4)))
  expect_identical(proportional_rates(1, 0, 2L, c(a = 1, b = 3)),
                   list(a = diag(2), b = diag(2) * 3))
  expect_error(proportional_rates(1, 0, 2L, c(1, 3)), "named by regime")
  expect_equal(shared_rates(1, 0, 2L, list(a = c(2, 2))), list(a = diag(2) * 2))
  expect_error(shared_rates(1, 0, 2L, list(a = 1)), "`eigenvalues` needs 2 values")
})

test_that("helpers of trait_mk() build the rates and the hidden classes", {
  skip_if_not_installed("ape")
  expect_identical(mk_rate_matrix(2, c("a", "b")),
                   list(q = matrix(c(0, 2, 2, 0), 2, dimnames = list(c("a", "b"), c("a", "b"))),
                        states = c("a", "b")))
  expect_error(mk_rate_matrix(matrix(1, 2, 2), NULL), "the same state names")
  expect_error(mk_rate_matrix(1, "a"), "two or more different names")
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  expect_null(mk_edge_rates(tree, NULL, NULL))
  expect_identical(mk_edge_rates(tree, NULL, c(1, 2, 3, 4)), c(1, 2, 3, 4))
  set.seed(1)
  drawn <- mk_edge_rates(tree, 2, NULL)
  expect_length(drawn, 4)
  expect_gt(min(drawn), 0)
  expect_error(mk_edge_rates(tree, 1, c(1, 2, 3, 4)), "not both")
  expect_error(mk_edge_rates(tree, NULL, 1), "one number, zero or more, per row")
  expect_identical(mk_classes(c(1, 2)), c(R1 = 1, R2 = 2))
  expect_identical(mk_classes(c(slow = 1, fast = 2)), c(slow = 1, fast = 2))
  expect_error(mk_classes(1), "two or more multipliers")
  expect_error(mk_classes(c(a = 1, a = 2)), "different names")
  set.seed(1)
  q <- mk_rate_matrix(1, c("a", "b"))$q
  hidden <- mk_hidden(q, c("a", "b"), "a", c(slow = 0, fast = 2), 0.5)
  expect_named(hidden, c("q", "full", "labels", "start"))
  expect_identical(hidden$labels, c("a|slow", "b|slow", "a|fast", "b|fast"))
  expect_identical(unname(hidden$q), matrix(c(0, 0, 0.5, 0, 0, 0, 0, 0.5, 0.5, 0, 0, 2, 0, 0.5,
                                              2, 0), 4))
  expect_true(hidden$start %in% c("a|slow", "a|fast"))
  tree <- ape::reorder.phylo(tree, "cladewise")
  run <- mk_history(tree, hidden$q, hidden$start, rep(1, 4))
  out <- mk_observed(tree, run, hidden, c("a", "b"))
  expect_s3_class(out, "simmap")
  expect_named(out$classes, tree$tip.label)
  expect_identical(diag(out$hidden$Q), -rowSums(hidden$q))
})

test_that("helpers of trait_rates() place the shifts and give the rates", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  expect_identical(given_edge_rates(c(a = 1, b = 2, c = 3, d = 4), 4), c(1, 2, 3, 4))
  expect_error(given_edge_rates(c(1, 2), 4), "one positive number per row")
  expect_identical(shift_places(tree, 0, NULL, NULL),
                   list(edges = integer(0), position = numeric(0)))
  expect_identical(shift_places(tree, 1, c(2, 4), NULL),
                   list(edges = c(2L, 4L), position = c(0, 0)))
  expect_error(shift_places(tree, 1, c(2, 2), NULL), "different row numbers")
  set.seed(1)
  drawn <- shift_places(tree, 2, NULL, NULL)
  expect_length(drawn$edges, 2)
  set.seed(1)
  along <- shift_places(tree, 1, NULL, 2)
  expect_named(along, c("edges", "position"))
  expect_length(along$position, length(along$edges))
  expect_true(all(along$position <= tree$edge.length[along$edges]))
  expect_identical(shift_sizes(4, 2), c(4, 4))
  expect_identical(shift_sizes(function(n) seq_len(n), 3), 1:3)
  expect_identical(shift_sizes("a", 0), character(0))
  expect_error(shift_sizes(0, 2), "one positive multiplier per shift")
  cladewise <- ape::reorder.phylo(tree, "cladewise", index.only = TRUE)
  # a shift at the start of the first edge passes to the edges below it
  first <- shift_rates(tree, cladewise, 1, 1L, 0, 3, TRUE, FALSE)
  expect_identical(first$edge_rate, c(3, 3, 3, 1))
  expect_identical(first$shifts, data.frame(edge = 1L, node = tree$edge[1, 2], position = 0,
                                            multiplier = 3, rate = 3))
  expect_identical(shift_rates(tree, cladewise, 1, 1L, 0, 3, TRUE, TRUE)$edge_rate,
                   c(3, 1, 1, 1))
  expect_identical(shift_rates(tree, cladewise, 2, 1L, 0.5, 3, FALSE, FALSE)$edge_rate,
                   c(4, 6, 6, 2))
  constant <- gbm_rates(tree, cladewise, 2, 0)
  expect_identical(constant$edge_rate, rep(2, 4))
  expect_identical(constant$node_rate, stats::setNames(rep(2, 5), 1:5))
})

test_that("the steps of shift_tree_once() keep the lineages and their segments", {
  s <- shift_tree_state(2)
  expect_identical(s$alive, 1L)
  expect_identical(s$rates, 2)
  set.seed(1)
  shift_tree_shift(s, 3)
  expect_identical(s$shifts, 1L)
  expect_identical(s$rate, 6)
  expect_identical(s$rates, c(2, 6))
  expect_identical(s$regime, 1L)
  expect_identical(s$segments[1, ], c(0, 1, 0, 0, 0, 1))
  # with no extinction every event is a speciation
  s$t <- 1
  shift_tree_event(s, 0)
  expect_identical(s$alive, 2:3)
  expect_identical(s$from, c(1L, 2L, 2L))
  expect_identical(s$rate, c(6, 6, 6))
  expect_identical(s$segments[2, ], c(1, 2, 0, 1, 1, 1))
  # the table of segments grows when it is full
  for (i in 1:70) shift_tree_close(s, 2L, 2, 2)
  expect_identical(s$n, 72L)
  expect_identical(nrow(s$segments), 128L)
  expect_null(shift_tree_sample(s, numeric(0), numeric(0), 1))
  # a run that spent no time with the number of tips wanted is never sampled
  expect_null(shift_tree_sample(s, 1, 1, 1))
  # with extinction far more likely than speciation the lineage goes extinct
  e <- shift_tree_state(1)
  set.seed(1)
  shift_tree_event(e, 1e9)
  expect_length(e$alive, 0)
  expect_identical(e$segments[1, 6], 2)
  set.seed(2)
  s <- shift_tree_state(1)
  for (i in 1:3) shift_tree_event(s, 0)
  out <- shift_tree_sample(s, 0.5, 1.5, 1e-9)
  expect_s3_class(out$tree, "phylo")
  expect_true(attr(out, "capped"))
})
