crown_tips <- function(tree) ape::extract.clade(tree, ape::Ntip(tree) + 2)$tip.label

test_that("rates per edge match constant rates when equal, and sum to one otherwise", {
  tree <- ape::read.tree(text = "((a:0.5,b:0.5):1,(c:1.2,(d:0.3,e:0.3):0.9):0.3);")
  patterns <- as.matrix(expand.grid(rep(list(0:1), 5)))
  colnames(patterns) <- tree$tip.label
  n_edge <- nrow(tree$edge)
  expect_equal(damocles_loglik(patterns, tree, rep(0.9, n_edge), rep(1.6, n_edge)),
               damocles_loglik(patterns, tree, 0.9, 1.6))
  gamma <- seq(0.2, 2, length.out = n_edge)
  mu <- rev(gamma)
  expect_equal(sum(exp(damocles_loglik(patterns, tree, gamma, mu, sympatry = 0.3))), 1)
  rooted <- tree
  rooted$root.edge <- 0.5
  expect_error(damocles_loglik(patterns, rooted, gamma, mu), "root edge")
  same_root <- gamma
  same_root[rooted$edge[, 1] == 6] <- 1
  expect_equal(sum(exp(damocles_loglik(patterns, rooted, same_root, mu * 0 + 1))), 1)
  expect_error(damocles_loglik(patterns, tree, gamma[-1], mu), "one per edge")
})

test_that("simulations with rates per edge have the probabilities of the likelihood", {
  tree <- ape::read.tree(text = "((a:0.5,b:0.5):1,(c:0.8,d:0.8):0.7);")
  patterns <- as.matrix(expand.grid(rep(list(0:1), 4)))
  colnames(patterns) <- tree$tip.label
  gamma <- c(0.2, 2, 0.5, 1, 3, 0.1)
  mu <- c(1, 0.3, 2, 0.6, 0.4, 1.5)
  expected <- exp(damocles_loglik(patterns, tree, gamma, mu))
  set.seed(10)
  sims <- simulate_damocles(tree, gamma, mu, n = 20000)
  key <- function(m) apply(m[, tree$tip.label, drop = FALSE], 1, paste, collapse = "")
  observed <- table(factor(key(sims), levels = key(patterns))) / nrow(sims)
  expect_lt(max(abs(as.vector(observed) - expected)), 0.012)
})

test_that("clade_regimes assigns edges to clades", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,(d:0.5,e:0.5):0.5):1);")
  reg <- clade_regimes(tree, list(cde = c("c", "d", "e"), de = c("d", "e")))
  expect_s3_class(reg, "factor")
  expect_identical(levels(reg), c("background", "cde", "de"))
  child <- tree$edge[, 2]
  expect_identical(as.character(reg[child == match("a", tree$tip.label)]), "background")
  expect_identical(as.character(reg[child == match("c", tree$tip.label)]), "cde")
  expect_identical(as.character(reg[child == match("d", tree$tip.label)]), "de")
  # the stem of the clade is background unless stem = TRUE
  stem <- child == ape::getMRCA(tree, c("c", "d", "e"))
  expect_identical(as.character(reg[stem]), "background")
  expect_identical(as.character(clade_regimes(tree, list(x = c("c", "d", "e")),
                                              stem = TRUE)[stem]), "x")
  expect_identical(as.character(clade_regimes(tree, list(x = "a"))[child == 1]), "x")
  expect_warning(clade_regimes(tree, list(x = c("a", "c"))), "not monophyletic")
  expect_error(clade_regimes(tree, list(x = "z")), "not in `tree`")
  expect_error(clade_regimes(tree, list("a")), "different names")
  expect_error(clade_regimes(tree, list(x = "a"), background = "x"), "background")
  expect_error(clade_regimes(list(), list(x = "a")), "phylo")
})

test_that("fit_damocles with regimes recovers a difference among clades", {
  set.seed(11)
  tree <- treeshapesim::sim_tree_shape(40)
  reg <- clade_regimes(tree, list(A = crown_tips(tree)))
  # a strong difference and many communities, so that the test rejects whatever the draw
  gamma <- ifelse(reg == "A", 4, 0.25)
  comm <- simulate_damocles(tree, gamma, 1, n = 20)
  constant <- fit_damocles(comm, tree, shared = TRUE)
  by_clade <- fit_damocles(comm, tree, shared = TRUE, regimes = reg, vary = "gamma")
  expect_identical(attr(by_clade, "parameters"), c("gamma_background", "gamma_A", "mu"))
  expect_gt(by_clade$gamma_A[1], by_clade$gamma_background[1])
  expect_gte(attr(by_clade, "loglik_total"),
             sum(damocles_loglik(comm, tree, gamma, 1)) - 1e-6)
  lrt <- damocles_lrt(constant, by_clade)
  expect_identical(lrt$df[2], 1)
  expect_lt(lrt$p_value[2], 0.05)
  expect_equal(lrt$n_par, c(2, 3))
  both <- fit_damocles(comm[1:2, ], tree, regimes = reg)
  expect_identical(attr(both, "parameters"),
                   c("gamma_background", "gamma_A", "mu_background", "mu_A"))
  expect_identical(attr(both, "n_par"), 4L)
  expect_identical(damocles_lrt(fit_damocles(comm[1:2, ], tree), both)$df[2], 4)
  # one regime is the constant model
  single <- fit_damocles(comm[1:2, ], tree, regimes = rep("x", nrow(tree$edge)))
  expect_identical(attr(single, "parameters"), c("gamma", "mu"))
  expect_error(fit_damocles(comm, tree, regimes = reg[-1]), "one regime for each edge")
  expect_error(fit_damocles(comm, tree, regimes = reg, edge_covariate = rep(1, length(reg))),
               "not both")
})

test_that("fit_damocles with an edge covariate estimates its slope", {
  set.seed(12)
  tree <- treeshapesim::sim_tree_shape(40)
  x <- exp(stats::rnorm(nrow(tree$edge)))
  z <- log(x) - mean(log(x))
  comm <- simulate_damocles(tree, 1 * exp(1.5 * z), 1, n = 10)
  fit <- fit_damocles(comm, tree, shared = TRUE, edge_covariate = x, vary = "gamma")
  expect_identical(attr(fit, "parameters"), c("gamma", "mu", "beta_gamma"))
  expect_gt(fit$beta_gamma[1], 0)
  rates <- attr(fit, "edge_rates")[[1]]
  expect_equal(rates$gamma, fit$gamma[1] * exp(fit$beta_gamma[1] * z))
  expect_equal(rates$mu, rep(fit$mu[1], nrow(tree$edge)))
  ci <- confint(fit, "beta_gamma")
  expect_identical(ci$community, "shared")
  expect_true(ci$lower < fit$beta_gamma[1] && fit$beta_gamma[1] < ci$upper)
  both <- fit_damocles(comm[1, , drop = FALSE], tree, edge_covariate = x, sympatry = NULL)
  expect_identical(attr(both, "parameters"), c("gamma", "mu", "beta_gamma", "beta_mu"))
  expect_error(fit_damocles(comm, tree, edge_covariate = -x), "positive")
  rooted <- tree
  rooted$root.edge <- 0.2
  expect_error(fit_damocles(comm, rooted, edge_covariate = x), "root edge")
  rooted_reg <- clade_regimes(rooted, list(A = crown_tips(rooted)), stem = TRUE)
  expect_error(fit_damocles(comm, rooted, regimes = rooted_reg), "root edge")
})

test_that("profile intervals cover the estimate and flag unbounded rates", {
  set.seed(13)
  tree <- treeshapesim::sim_tree_shape(30)
  comm <- simulate_damocles(tree, 1, 1, n = 2)
  fit <- fit_damocles(comm, tree)
  ci <- confint(fit, "gamma")
  expect_identical(ci$community, rownames(comm))
  expect_true(all(ci$lower <= ci$estimate & ci$estimate <= ci$upper))
  expect_true(all(ci$lower >= 0))
  wide <- confint(fit_damocles(comm, tree, shared = TRUE, sympatry = NULL), "sympatry",
                  level = 0.99)
  expect_true(wide$lower >= 0 && wide$upper <= 1)
  expect_error(confint(fit, "beta_gamma"), "not free parameters")
  expect_error(confint(fit, level = 2), "between 0 and 1")
  flat <- confint(fit_damocles(c(stats::setNames(rep(1, 30), tree$tip.label)), tree), "mu")
  expect_true(flat$lower == 0 || is.infinite(flat$upper))
})

test_that("damocles_lrt checks its inputs", {
  set.seed(14)
  tree <- treeshapesim::sim_tree_shape(20)
  comm <- simulate_damocles(tree, 1, 1, n = 3)
  separate <- fit_damocles(comm, tree)
  shared <- fit_damocles(comm, tree, shared = TRUE)
  expect_identical(damocles_lrt(shared, separate)$df[2], 4)
  expect_error(damocles_lrt(separate, shared), "more free parameters")
  expect_error(damocles_lrt(shared, fit_damocles(comm[1:2, ], tree)), "same communities")
  expect_error(damocles_lrt(shared, data.frame()), "fits of")
})

test_that("damocles_test simulates with regimes and with rates given per edge", {
  set.seed(15)
  tree <- treeshapesim::sim_tree_shape(30)
  reg <- clade_regimes(tree, list(A = crown_tips(tree)))
  comm <- simulate_damocles(tree, ifelse(reg == "A", 2, 0.5), 1, n = 3)
  res <- damocles_test(comm, tree, n_sim = 19, shared = TRUE, regimes = reg, vary = "gamma")
  expect_true(all(c("gamma_background", "gamma_A", "mu") %in% names(res)))
  expect_output(print(res), "gamma_A")
  given <- damocles_test(comm, tree, n_sim = 19, gamma = list(ifelse(reg == "A", 2, 0.5)),
                         mu = 1)
  expect_true(all(is.na(given$gamma)))
  expect_equal(given$mu, rep(1, 3))
  expect_error(damocles_test(comm, tree, gamma = 1, mu = 1, regimes = reg), "estimated rates")
})
