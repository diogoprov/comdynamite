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
  # not a clade at the root, whose only random placements may split the tree in the same way
  nodes <- ape::Ntip(tree) + c(3, 6)
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
  expect_identical(res$samples$loglik[1], attr(damocles_fit(comm, tree, shared = TRUE,
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

test_that("damocles_samples integrates the likelihood over configurations", {
  set.seed(21)
  tree <- treeshapesim::sim_tree_shape(30)
  x <- exp(stats::rnorm(nrow(tree$edge), 0, 0.5))
  comm <- simulate_damocles(tree, 2 * x / exp(mean(log(x))), 1, n = 4)
  # with copies of one configuration, the integrated fit is the fit of that configuration
  single <- damocles_fit(comm, tree, shared = TRUE, edge_covariate = x, vary = "gamma")
  same <- damocles_samples(comm, tree, edge_covariate = list(x, x), vary = "gamma",
                           integrate = TRUE)
  expect_equal(same$integrated$loglik, attr(single, "loglik_total"), tolerance = 1e-4)
  expect_equal(unname(same$integrated$estimates["beta_gamma"]), single$beta_gamma[1],
               tolerance = 1e-2)
  expect_identical(same$integrated$df, 1)
  expect_equal(same$integrated$statistic,
               2 * (same$integrated$loglik - attr(same$constant, "loglik_total")))
  # averaging likelihoods: never above the best configuration, never below the mean
  noisy <- lapply(1:3, function(i) x * exp(stats::rnorm(length(x), 0, 0.5)))
  res <- damocles_samples(comm, tree, edge_covariate = noisy, vary = "gamma",
                          integrate = TRUE)
  expect_lte(res$integrated$loglik, max(res$samples$loglik) + 1e-6)
  expect_named(res$integrated, c("estimates", "loglik", "n_par", "delta_aic", "statistic",
                                 "df", "p_value", "convergence"))
  expect_output(print(res), "Integrated over configurations")
  expect_output(print(summary(res)), "beta_gamma")
  expect_null(damocles_samples(comm, tree, edge_covariate = noisy[1:2])$integrated)
  # regimes with the same levels are integrated; different levels are not
  crown <- ape::Ntip(tree) + 2
  reg <- lapply(c(crown, crown + 1), function(n) clade_regimes(tree, list(A = n)))
  r <- damocles_samples(comm, tree, regimes = reg, vary = "gamma", integrate = TRUE)
  expect_named(r$integrated$estimates, c("gamma_background", "gamma_A", "mu"))
  expect_error(damocles_samples(comm, tree, regimes = list(reg[[1]], rep("x", length(x))),
                                integrate = TRUE), "same regimes")
  expect_error(damocles_samples(comm, tree, regimes = reg, integrate = NA), "integrate")
  # configurations with one regime each: the constant model, with no test
  flat <- damocles_samples(comm, tree, regimes = list(rep("x", length(x)), rep("y", length(x))),
                           integrate = TRUE)
  expect_true(is.na(flat$integrated$p_value))
  expect_equal(flat$integrated$loglik, attr(flat$constant, "loglik_total"), tolerance = 1e-4)
  expect_equal(log_mean_exp(c(log(2), log(4))), log(3))
  expect_identical(log_mean_exp(c(-1e10, -Inf)), -1e10)
  expect_equal(log_mean_exp(c(log(2), -1e10)), log(1))
})

test_that("damocles_partition_test compares a partition with random ones", {
  set.seed(18)
  tree <- treeshapesim::sim_tree_shape(30)
  reg <- clade_regimes(tree, list(A = ape::Ntip(tree) + 2))
  comm <- simulate_damocles(tree, ifelse(reg == "A", 3, 0.5), 1, n = 5)
  res <- expect_warning_value(damocles_partition_test(comm, tree, reg, n_rand = 9, vary = "gamma"),
                              "cannot be below")
  expect_s3_class(res, "mc_damocles_partition")
  expect_length(res$null, res$n_distinct)
  expect_lte(res$n_distinct, 9)
  # each distinct partition counts once
  expect_equal(res$p_random, (sum(res$null >= res$statistic - 1e-8) + 1) / (res$n_distinct + 1))
  expect_identical(res$df, 1)
  expect_output(print(res), paste0("(", res$n_distinct, " distinct random partitions)"),
                fixed = TRUE)
  expect_error(damocles_partition_test(comm, tree, rep("x", length(reg))), "two or more")
  expect_error(damocles_partition_test(comm, tree, reg, n_rand = 0), "`n_rand` or `n_boot`")
  expect_error(damocles_partition_test(comm, tree, reg, n_rand = -1), "n_rand")
  expect_error(damocles_partition_test(comm, tree, reg, n_boot = 1.5), "n_boot")
  expect_true(is.na(res$p_boot))
  expect_length(res$null_boot, 0)
  # The partitions are drawn before the fits: the p-value does not depend on the cores.
  set.seed(3)
  one <- suppressWarnings(damocles_partition_test(comm, tree, reg, n_rand = 9, vary = "gamma"))
  set.seed(3)
  two <- suppressWarnings(damocles_partition_test(comm, tree, reg, n_rand = 9, vary = "gamma",
                                                  cores = 2))
  expect_identical(one$null, two$null)
})

test_that("damocles_partition_test runs a parametric bootstrap under constant rates", {
  set.seed(18)
  tree <- treeshapesim::sim_tree_shape(30)
  reg <- clade_regimes(tree, list(A = ape::Ntip(tree) + 2))
  comm <- simulate_damocles(tree, ifelse(reg == "A", 3, 0.5), 1, n = 5)
  set.seed(5)
  res <- damocles_partition_test(comm, tree, reg, n_rand = 0, n_boot = 9, vary = "gamma")
  expect_true(is.na(res$p_random))
  expect_identical(res$n_distinct, 0L)
  expect_length(res$null_boot, 9)
  expect_true(all(res$null_boot > -1e-6))
  expect_equal(res$p_boot, (sum(res$null_boot >= res$statistic - 1e-8) + 1) / 10)
  out <- capture.output(print(res))
  expect_true(any(grepl("(9 simulations under constant rates)", out, fixed = TRUE)))
  expect_false(any(grepl("random partitions)", out, fixed = TRUE)))
  # The communities are simulated before the fits: the p-value does not depend on the cores.
  set.seed(5)
  two <- damocles_partition_test(comm, tree, reg, n_rand = 0, n_boot = 9, vary = "gamma",
                                 cores = 2)
  expect_identical(res$null_boot, two$null_boot)
})

test_that("bootstrap fits that fail are left out", {
  set.seed(18)
  tree <- treeshapesim::sim_tree_shape(30)
  reg <- clade_regimes(tree, list(A = ape::Ntip(tree) + 2))
  comm <- simulate_damocles(tree, ifelse(reg == "A", 3, 0.5), 1, n = 5)
  real_sim <- simulate_damocles
  count <- new.env()
  count$calls <- 0
  local_mocked_bindings(simulate_damocles = function(...) {
    count$calls <- count$calls + 1
    out <- real_sim(...)
    if (count$calls == 1) out[] <- NA  # a set that cannot be fitted
    out
  })
  set.seed(6)
  res <- expect_warning_value(damocles_partition_test(comm, tree, reg, n_rand = 0, n_boot = 3,
                                                      vary = "gamma"),
                              "1 of 3 bootstrap fits failed")
  expect_length(res$null_boot, 2)
  local_mocked_bindings(simulate_damocles = function(...) {
    out <- real_sim(...); out[] <- NA; out
  })
  res <- expect_warning_value(damocles_partition_test(comm, tree, reg, n_rand = 0, n_boot = 2,
                                                      vary = "gamma"),
                              "2 of 2")
  expect_true(is.na(res$p_boot))
})

test_that("random partitions that split the tree as the original are excluded", {
  # A balanced tree: the clade of one daughter of the root has the other as its complement.
  tree <- ape::read.tree(text = paste0("(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,",
                                       "(g:1,h:1):1):1);"))
  # Regime A on the clade of node 10 and its stem, as a shift on that edge.
  reg <- factor(ifelse(seq_len(nrow(tree$edge)) %in% 1:7, "A", "0"))
  expect_identical(adjusted_rand(reg, reg), 1)
  swapped <- factor(ifelse(reg == "A", "0", "A"))
  expect_identical(adjusted_rand(reg, swapped), 1)
  expect_identical(adjusted_rand(rep("a", 4), rep("b", 4)), 1)  # one block each: same partition
  expect_error(random_regimes(tree, reg), "matching `regimes`")
  set.seed(4)
  comm <- simulate_damocles(tree, 1, 1, n = 4)
  res <- expect_warning_value(damocles_partition_test(comm, tree, reg, n_rand = 5, vary = "gamma"),
                              "no random partition")
  expect_true(is.na(res$p_random))
  expect_length(res$null, 0)
  # With max_overlap = 1 the complement and the clade itself are allowed: two partitions.
  res <- expect_warning_value(damocles_partition_test(comm, tree, reg, n_rand = 5, vary = "gamma",
                                                      max_overlap = 1),
                              "distinct random partitions")
  expect_gte(res$p_random, 1 / 3)
})

test_that("dr_rates follows the equal-splits definition", {
  tree <- ape::read.tree(text = "((a:1,b:1):2,(c:2,(d:1,e:1):1):1);")
  dr <- dr_rates(tree)
  expect_named(dr, tree$tip.label)
  expect_equal(dr[["a"]], 1 / (1 + 2 / 2))
  expect_equal(dr[["d"]], 1 / (1 + 1 / 2 + 1 / 4))
  expect_equal(dr[["c"]], 1 / (2 + 1 / 2))
  per_edge <- dr_rates(tree, per = "edge")
  expect_length(per_edge, nrow(tree$edge))
  terminal <- tree$edge[, 2] <= ape::Ntip(tree)
  expect_identical(per_edge[terminal], unname(dr[tree$edge[terminal, 2]]))
  ab <- tree$edge[, 2] == ape::getMRCA(tree, c("a", "b"))
  expect_equal(per_edge[ab], sqrt(dr[["a"]] * dr[["b"]]))
  expect_error(dr_rates(list()), "phylo")
  expect_error(dr_rates(ape::read.tree(text = "(a,b);")), "branch lengths")
})
