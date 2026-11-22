# Branches of the functions split into helpers (phylogenetic tests, argument checks and plots)
# that the other test files do not reach.

partition_data <- function() {
  set.seed(1)
  trait <- stats::setNames(seq(0, 100, 10), letters[1:11])
  env <- runif(40, 0, 100)
  comm <- outer(env, trait, function(e, t) exp(-(e - t)^2 / 800))
  colnames(comm) <- names(trait)
  list(comm = comm, trait = trait, env = env, space = cbind(runif(40), runif(40)))
}

test_that("the checks of whole numbers and of numbers keep their rules", {
  expect_true(is_whole_number(3))
  expect_false(is_whole_number(2.5))
  expect_false(is_whole_number(c(1, 2)))
  expect_false(is_whole_number(Inf))
  expect_false(is_whole_number("1"))
  expect_true(is_one_number(Inf, finite = FALSE))
  expect_false(is_one_number(Inf, finite = TRUE))
  expect_false(is_one_number(NA_real_, finite = FALSE))
  expect_true(number_within(0, 0, 1, c(FALSE, FALSE)))
  expect_false(number_within(0, 0, 1, c(TRUE, FALSE)))
  expect_false(number_within(1, 0, 1, c(FALSE, TRUE)))
  expect_identical(check_count(2, "n"), 2)
  expect_null(check_count(NULL, "n", null_ok = TRUE))
  expect_error(check_count(1.5, "n"), "`n` must be one whole number, 1 or more", fixed = TRUE)
  expect_error(check_count(11, "n", max = 10),
               "`n` must be one whole number between 1 and 10", fixed = TRUE)
  expect_identical(check_number(Inf, "x", finite = FALSE), Inf)
  expect_error(check_number(Inf, "x"), "`x` must be one finite number", fixed = TRUE)
  expect_error(check_number(0, "x", lower = 0, upper = 1, open = c(TRUE, FALSE)),
               "`x` must be one finite number greater than 0 and up to 1", fixed = TRUE)
})

test_that("number_range() writes each kind of range", {
  closed <- c(FALSE, FALSE)
  expect_identical(number_range(-Inf, Inf, closed, TRUE), "one finite number")
  expect_identical(number_range(-Inf, Inf, closed, FALSE), "one number")
  expect_identical(number_range(0, 1, closed, TRUE), "one finite number between 0 and 1")
  expect_identical(number_range(0, 1, c(FALSE, TRUE), TRUE),
                   "one finite number from 0 and less than 1")
  expect_identical(number_range(0, Inf, closed, TRUE), "one finite number, 0 or more")
  expect_identical(number_range(0, Inf, c(TRUE, FALSE), TRUE),
                   "one finite number, greater than 0")
  expect_identical(number_range(-Inf, 1, closed, FALSE), "one number, 1 or less")
  expect_identical(number_range(-Inf, 1, c(FALSE, TRUE), TRUE),
                   "one finite number, less than 1")
})

test_that("seeds and priors of ABC are checked", {
  expect_identical(check_seed(12), 12)
  expect_error(check_seed(1.5), "`seed` must be one whole number between", fixed = TRUE)
  expect_error(check_seed(.Machine$integer.max, "start"), "`start` must be one whole number")
  expect_true(distinct_names(list(a = 1, b = 2)))
  expect_false(distinct_names(list(a = 1, a = 2)))
  expect_false(distinct_names(list(1, 2)))
  expect_true(is_increasing_pair(c(0, 1)))
  expect_false(is_increasing_pair(c(1, 0)))
  expect_false(is_increasing_pair(c(0, Inf)))
  expect_identical(check_prior(list(a = c(0, 1))), list(a = c(0, 1)))
  expect_error(check_prior(list(a = c(1, 0))), "increasing pair")
})

test_that("phylo_vector() takes the deprecated `predictor` and one variable only", {
  set.seed(1)
  tree <- ape::rcoal(10)
  env <- runif(30, 0, 100)
  optima <- stats::setNames(seq(0, 100, length.out = 10), tree$tip.label)
  comm <- outer(env, optima, function(e, o) as.numeric(abs(e - o) < 25))
  comm[rowSums(comm) == 0, 1] <- 1
  set.seed(2)
  expect_warning(old <- phylo_vector(comm, tree, predictor = env, n_perm = 9, n_sim = 5),
                 "`predictor` is deprecated, use `env`", fixed = TRUE)
  set.seed(2)
  new <- phylo_vector(comm, tree, env, n_perm = 9, n_sim = 5)
  expect_identical(old, new)
  expect_s3_class(new, "mc_phylo_vector")
  expect_error(phylo_vector(comm, tree, cbind(env, rev(env)), n_perm = 9),
               "`env` must be a single variable", fixed = TRUE)
})

test_that("the eigenvectors of phylo_vector() are selected down to the first that fails", {
  observed <- c(0.1, -0.8, 0.5)
  nulls <- list(site = cbind(rep(0.2, 9), rep(0.1, 9), rep(0.6, 9)),
                taxa = cbind(rep(0.2, 9), rep(0.1, 9), rep(0.1, 9)))
  expect_identical(selected_vectors(observed, nulls, 0.2), 2L)
  expect_length(selected_vectors(observed, nulls, 0.05), 0)
})

test_that("trait_partition() refuses traits that are missing, infinite or absent", {
  d <- partition_data()
  expect_error(trait_partition(d$comm, NULL, d$env, d$space), "`traits` is NULL", fixed = TRUE)
  expect_error(trait_partition(d$comm, replace(d$trait, 1, Inf), d$env, d$space),
               "trait_partition(): `traits` has infinite values", fixed = TRUE)
  expect_error(trait_partition(d$comm, replace(d$trait, 1, NA), d$env, d$space),
               "trait_partition(): `traits` has missing values", fixed = TRUE)
  text <- stats::setNames(letters[1:11], letters[1:11])
  expect_error(trait_partition(d$comm, text, d$env, d$space), "`traits` is not numeric",
               fixed = TRUE)
  expect_error(trait_partition(d$comm, replace(d$trait, 1:11, 1), d$env, d$space),
               "traits without variance among the species with individuals", fixed = TRUE)
  traits <- partition_traits(d$trait)
  expect_identical(dim(traits), c(11L, 1L))
  matched <- partition_species(d$comm, traits, FALSE)
  expect_named(matched, c("comm", "traits"))
  expect_identical(rownames(matched$traits), colnames(matched$comm))
})

test_that("the parts of trait_partition() agree with its result", {
  d <- partition_data()
  set.seed(3)
  out <- trait_partition(d$comm, d$trait, d$env, d$space, n_perm = 9)
  traits <- partition_traits(d$trait)
  m <- partition_model(d$comm, traits, check_env(d$env, d$comm), check_env(d$space, d$comm),
                       rep(TRUE, 40), "mean")
  expect_named(m, c("comm", "traits", "env", "space", "landscape_env", "landscape_space",
                    "site_weights", "total", "component", "response"))
  ss <- partition_ss(m$response, m$env, m$space, m$site_weights)
  expect_named(ss, c("abc", "ab", "bc"))
  expect_identical(unname(ss / m$total), out$fractions$r2[1:3])
  observed <- partition_statistics(m, m$response, m$env, m$space)
  expect_identical(unname(observed[c("abc", "ab", "bc", "a", "b", "c")]),
                   out$fractions$statistic[1:6])
  expect_identical(out$shared, shared_status("none", list(), m, 9, 0.05))
  expect_identical(out$permutations$moved, "together")
})

test_that("plot_node_test() draws a tree without branch lengths", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  out <- node_test(daphniids$comm, node_contrasts(daphniids$tree_leibold), daphniids$env,
                   n_perm = 9)
  tree <- daphniids$tree_leibold
  tree$edge.length <- NULL
  p <- expect_ggplot(plot_node_test(out, tree, test = "site", level = 0.5))
  with_lengths <- tree
  with_lengths$edge.length <- rep(1, nrow(tree$edge))
  q <- plot_node_test(out, with_lengths, test = "site", level = 0.5)
  expect_identical(length(p$layers), length(q$layers))
  for (i in seq_along(p$layers)) {
    expect_identical(ggplot2::layer_data(p, i), ggplot2::layer_data(q, i))
  }
  expect_identical(attr(p, "marks"), attr(q, "marks"))
})

test_that("the helpers of the plots of trees give their coordinates", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  layout <- tree_layout(tree)
  tips <- tree_tips(tree, layout)
  expect_named(tips, c("x", "y", "label"))
  expect_identical(tips$label, c("a", "b", "c"))
  expect_identical(tips$x, c(2, 2, 2))
  turned <- fan_turn(c(1, 3), 4)
  expect_identical(turned$angle, c(45, 45))
  expect_identical(turned$flip, c(FALSE, TRUE))
  hung <- hanging_tree(tree, layout, 9)
  expect_named(hung, c("x", "xend", "y", "yend"))
  expect_identical(nrow(hung), 8L)
  expect_identical(max(hung$y), 12.5)
  expect_true(has_edge_rates(list(tree = tree, edge_rates = 1:4)))
  expect_false(has_edge_rates(tree))
  given <- rates_of_tree(list(tree = tree, edge_rates = 1:4), NULL)
  expect_identical(given$rates, 1:4)
  expect_s3_class(given$tree, "phylo")
})
