test_that("clade_regimes takes node numbers", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,(d:0.5,e:0.5):0.5):1);")
  by_tips <- clade_regimes(tree, list(cde = c("c", "d", "e")))
  by_node <- clade_regimes(tree, list(cde = ape::getMRCA(tree, c("c", "d", "e"))))
  expect_identical(by_node, by_tips)
  expect_error(clade_regimes(tree, list(x = 99)), "one node")
})

test_that("random_regimes keeps the regimes and the sizes of the clades", {
  set.seed(16)
  tree <- treeshapesim::sim_tree_shape(50)
  size <- damocles_subtree_size(tree)
  expect_identical(size[tree$edge[, 2] <= ape::Ntip(tree)], rep(1, ape::Ntip(tree)))
  expect_equal(sum(size[tree$edge[, 1] == ape::Ntip(tree) + 1]), nrow(tree$edge))
  nodes <- ape::Ntip(tree) + c(2, 5)
  reg <- clade_regimes(tree, list(A = nodes[1], B = nodes[2]), stem = TRUE)
  for (i in 1:5) {
    r <- random_regimes(tree, reg)
    expect_identical(levels(r), levels(reg))
    expect_length(r, length(reg))
  }
  # a clade of n edges goes to subtrees of n / 2 to 2 n edges
  edge <- which.min(abs(size - 10))
  one <- clade_regimes(tree, list(A = tree$edge[edge, 2]), stem = TRUE)
  moved <- replicate(20, random_regimes(tree, one, max_overlap = 1) == "A")
  expect_true(all(colSums(moved) >= size[edge] / 2 & colSums(moved) <= 2 * size[edge]))
  # with the default, replacements share at most a quarter of their edges with the clade
  original <- one == "A"
  for (i in 1:10) {
    placed <- random_regimes(tree, one) == "A"
    expect_lte(sum(placed & original) / sum(placed | original), 0.25)
  }
  subtree <- damocles_subtree_edges(tree)
  expect_identical(colSums(subtree), size)
  expect_error(random_regimes(tree, one, max_overlap = 2), "max_overlap")
  expect_error(random_regimes(tree, reg[-1]), "one regime for each edge")
  expect_error(random_regimes(tree, reg, tolerance = 0.5), "tolerance")
  expect_error(random_regimes(list(), reg), "phylo")
})

test_that("damocles_samples fits each configuration and summarises them", {
  set.seed(17)
  tree <- treeshapesim::sim_tree_shape(30)
  crown <- ape::Ntip(tree) + 2
  reg <- clade_regimes(tree, list(A = crown))
  comm <- simulate_damocles(tree, ifelse(reg == "A", 3, 0.5), 1, n = 5)
  configurations <- list(reg, clade_regimes(tree, list(A = crown + 1)),
                         rep("x", nrow(tree$edge)))
  res <- damocles_samples(comm, tree, regimes = configurations, vary = "gamma")
  expect_s3_class(res, "mc_damocles_samples")
  expect_identical(nrow(res$samples), 3L)
  expect_identical(dim(res$gamma), c(nrow(tree$edge), 3L))
  # one regime is the constant model, with no test
  expect_true(is.na(res$samples$p_value[3]))
  expect_equal(res$samples$loglik[3], attr(res$constant, "loglik_total"), tolerance = 1e-6)
  expect_equal(res$samples$loglik[1], attr(damocles_fit(comm, tree, shared = TRUE,
                                                        regimes = reg, vary = "gamma"),
                                           "loglik_total"))
  expect_output(print(res), "3 configurations")
  s <- summary(res)
  expect_identical(rownames(s$parameters), "mu")
  expect_identical(nrow(s$edges), nrow(tree$edge))
  expect_output(print(s), "Parameters shared")
  covariates <- lapply(1:2, function(i) exp(stats::rnorm(nrow(tree$edge))))
  cv <- damocles_samples(comm[1:2, ], tree, edge_covariate = covariates, vary = "gamma",
                         cores = 2)
  expect_identical(rownames(summary(cv)$parameters), c("gamma", "mu", "beta_gamma"))
  both <- damocles_samples(comm[1:2, ], tree, regimes = configurations[1:2])
  expect_null(summary(both)$parameters)
  expect_output(print(summary(both)), "Rates of each edge")
  expect_error(damocles_samples(comm, tree), "give a list")
  expect_error(damocles_samples(comm, tree, regimes = reg), "list of configurations")
  expect_error(damocles_samples(comm, tree, regimes = configurations, cores = 0), "cores")
})

test_that("damocles_partition_test compares a partition with random ones", {
  set.seed(18)
  tree <- treeshapesim::sim_tree_shape(30)
  reg <- clade_regimes(tree, list(A = ape::Ntip(tree) + 2))
  comm <- simulate_damocles(tree, ifelse(reg == "A", 3, 0.5), 1, n = 5)
  res <- damocles_partition_test(comm, tree, reg, n_rand = 9, vary = "gamma")
  expect_s3_class(res, "mc_damocles_partition")
  expect_length(res$null, 9)
  expect_true(res$p_random > 0 && res$p_random <= 1)
  expect_identical(res$df, 1)
  expect_output(print(res), "random partitions")
  expect_error(damocles_partition_test(comm, tree, rep("x", length(reg))), "two or more")
  expect_error(damocles_partition_test(comm, tree, reg, n_rand = 0), "n_rand")
})

test_that("dr_rates follows the equal-splits definition", {
  tree <- ape::read.tree(text = "((a:1,b:1):2,(c:2,(d:1,e:1):1):1);")
  dr <- dr_rates(tree)
  expect_identical(names(dr), tree$tip.label)
  expect_equal(dr[["a"]], 1 / (1 + 2 / 2))
  expect_equal(dr[["d"]], 1 / (1 + 1 / 2 + 1 / 4))
  expect_equal(dr[["c"]], 1 / (2 + 1 / 2))
  per_edge <- dr_rates(tree, per = "edge")
  expect_length(per_edge, nrow(tree$edge))
  terminal <- tree$edge[, 2] <= ape::Ntip(tree)
  expect_equal(per_edge[terminal], unname(dr[tree$edge[terminal, 2]]))
  ab <- tree$edge[, 2] == ape::getMRCA(tree, c("a", "b"))
  expect_equal(per_edge[ab], sqrt(dr[["a"]] * dr[["b"]]))
  expect_error(dr_rates(list()), "phylo")
  expect_error(dr_rates(ape::read.tree(text = "(a,b);")), "branch lengths")
})
