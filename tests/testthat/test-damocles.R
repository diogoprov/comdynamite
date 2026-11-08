all_patterns <- function(tips) {
  m <- as.matrix(expand.grid(rep(list(0:1), length(tips))))
  dimnames(m) <- list(NULL, tips)
  m
}

test_that("the likelihood of two sister tips follows the speciation rule", {
  tree <- ape::read.tree(text = "(a:1.3,b:1.3);")
  gamma <- 0.7
  mu <- 0.4
  s <- gamma + mu
  decay <- exp(-s * 1.3)
  p01 <- gamma * (1 - decay) / s
  p11 <- (gamma + mu * decay) / s
  # a present ancestor leaves one present daughter; an absent one, two absent daughters
  present <- 0.5 * (p11 * (1 - p01) + p01 * (1 - p11))
  absent <- p01 * (1 - p01)
  expect_equal(damocles_loglik(c(a = 1, b = 0), tree, gamma, mu, root = "present"),
               log(present), ignore_attr = TRUE)
  expect_equal(damocles_loglik(c(a = 1, b = 0), tree, gamma, mu, root = "absent"),
               log(absent), ignore_attr = TRUE)
  expect_equal(damocles_loglik(c(a = 1, b = 0), tree, gamma, mu),
               log((present + absent) / 2), ignore_attr = TRUE)
})

test_that("likelihoods of all communities sum to one, also with a root edge", {
  tree <- ape::read.tree(text = "((a:0.5,b:0.5):1,(c:1.2,(d:0.3,e:0.3):0.9):0.3);")
  patterns <- all_patterns(tree$tip.label)
  for (root in c("both", "absent", "present")) {
    expect_equal(sum(exp(damocles_loglik(patterns, tree, 0.9, 1.6, root = root))), 1)
  }
  tree$root.edge <- 0.8
  expect_equal(sum(exp(damocles_loglik(patterns, tree, 0.2, 3))), 1)
  expect_equal(sum(exp(damocles_loglik(patterns, tree, 0, 0))), 1)
})

test_that("simulated communities have the probabilities of the likelihood", {
  tree <- ape::read.tree(text = "((a:0.5,b:0.5):1,(c:0.8,d:0.8):0.7);")
  patterns <- all_patterns(tree$tip.label)
  expected <- exp(damocles_loglik(patterns, tree, 0.6, 0.9))
  set.seed(1)
  sims <- simulate_damocles(tree, 0.6, 0.9, n = 20000)
  key <- function(m) apply(m[, tree$tip.label, drop = FALSE], 1, paste, collapse = "")
  observed <- table(factor(key(sims), levels = key(patterns))) / nrow(sims)
  expect_lt(max(abs(as.vector(observed) - expected)), 0.012)
})

test_that("limits of the rates give the expected communities", {
  set.seed(2)
  tree <- treeshapesim::sim_tree_shape(15)
  expect_true(all(simulate_damocles(tree, 0, 2, n = 50, root = "absent") == 0))
  # without colonisation or extinction the presence of the root passes to one tip
  expect_true(all(rowSums(simulate_damocles(tree, 0, 0, n = 50, root = "present")) == 1))
  sims <- simulate_damocles(tree, 1, 1, n = 3)
  expect_identical(dim(sims), c(3L, 15L))
  expect_identical(colnames(sims), tree$tip.label)
  expect_true(all(sims %in% 0:1))
  # one present tip is impossible without colonisation from an absent root
  one <- stats::setNames(c(1, rep(0, 14)), tree$tip.label)
  expect_identical(unname(damocles_loglik(one, tree, 0, 1, root = "absent")), -Inf)
})

test_that("fit_damocles reaches at least the likelihood of the true rates", {
  set.seed(3)
  tree <- treeshapesim::sim_tree_shape(80)
  comm <- simulate_damocles(tree, 1, 0.5, n = 2)
  fit <- fit_damocles(comm, tree)
  expect_s3_class(fit, "mc_damocles_fit")
  expect_identical(rownames(fit), rownames(comm))
  expect_equal(fit$richness, unname(rowSums(comm)))
  expect_true(all(fit$loglik >= damocles_loglik(comm, tree, 1, 0.5) - 1e-6))
  expect_true(all(fit$convergence == 0))
  expect_output(print(fit), "maximum likelihood")
  one <- fit_damocles(comm, tree, start = c(1, 1))
  expect_equal(one$loglik, fit$loglik, tolerance = 1e-4)
})

test_that("damocles_test compares each community with both null models", {
  set.seed(4)
  tree <- treeshapesim::sim_tree_shape(30)
  comm <- simulate_damocles(tree, 1, 1, n = 2)
  rownames(comm) <- c("north", "south")
  res <- damocles_test(comm, tree, n_sim = 49, gamma = 1, mu = 1)
  expect_s3_class(res, "mc_damocles_test")
  expect_identical(rownames(res), c("north", "south"))
  expect_true(all(is.na(res$loglik)))
  expect_true(all(res$p_mpd_random > 0 & res$p_mpd_random <= 1))
  expect_true(all(res$n_valid <= 49))
  expect_output(print(res), "DAMOCLES null model")
  fitted <- damocles_test(comm[1, , drop = FALSE], tree, n_sim = 19)
  expect_true(is.finite(fitted$loglik))
  expect_error(damocles_test(comm, tree, gamma = 1), "both")
  expect_error(damocles_test(comm, tree, n_sim = 1), "n_sim")
  expect_error(damocles_test(comm, tree, gamma = 1:3, mu = 1), "one per community")
})

test_that("inputs are checked", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  expect_error(simulate_damocles(tree, -1, 1), "zero or more")
  expect_error(simulate_damocles(tree, 1, 1, n = 0), "`n`")
  expect_error(simulate_damocles(ape::read.tree(text = "(a:1,b:1,c:1);"), 1, 1),
               "dichotomous")
  expect_error(simulate_damocles(ape::read.tree(text = "((a,b),c);"), 1, 1), "branch lengths")
  expect_error(damocles_loglik(c(1, 0, 1, 0), tree, 1, 1), "species names")
  expect_warning(ll <- damocles_loglik(c(a = 1, b = 0, c = 1), tree, 1, 1), "counted as absent")
  expect_equal(ll, damocles_loglik(c(a = 1, b = 0, c = 1, d = 0), tree, 1, 1))
  expect_warning(damocles_loglik(c(a = 1, b = 0, c = 1, d = 0, z = 1), tree, 1, 1),
                 "dropped")
  expect_error(fit_damocles(c(a = 1, b = 0, c = 1, d = 0), tree, start = c(1, -1)), "start")
})
