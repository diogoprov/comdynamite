# a hand-built metacommunity: species a stays in site1, b leaves and comes back, c never occurs
fake_sim <- function(alpha = diag(3)) {
  h <- array(0L, c(2, 3, 3), dimnames = list(c("site1", "site2"), c("a", "b", "c"),
                                             c("step0", "step10", "step20")))
  h["site1", "a", ] <- 5L
  h["site1", "b", ] <- c(3L, 0L, 3L)
  h["site2", "a", ] <- c(4L, 4L, 0L)
  dimnames(alpha) <- list(c("a", "b", "c"), c("a", "b", "c"))
  structure(list(history = h, alpha = alpha), class = "mc_sim")
}

test_that("partner_turnover() counts entries and exits per site and step", {
  x <- fake_sim()
  # a: site1 over 20 steps, b leaves and returns (2 changes); site2 over 10 steps, no change
  expect_equal(partner_turnover(x), c(a = 2 / 30, b = 0, c = 0))
  # without competition between a and b there are no partners to count
  expect_equal(partner_turnover(x, "competitors"), c(a = 0, b = 0, c = 0))
  alpha <- diag(3)
  alpha[2, 1] <- 0.1
  expect_equal(partner_turnover(fake_sim(alpha), "competitors"), c(a = 2 / 30, b = 0, c = 0))
  # without step names each interval counts as one step
  y <- x
  dimnames(y$history) <- c(dimnames(y$history)[1:2], list(NULL))
  expect_equal(partner_turnover(y), c(a = 2 / 3, b = 0, c = 0))
})

test_that("turnover_rates() raises only the tip branches", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  x <- fake_sim()
  r <- turnover_rates(x, tree, sigma2 = 2, k = 3)
  tip_edge <- match(1:3, tree$edge[, 2])
  expect_equal(r[tip_edge], 2 * (1 + 3 * c(2 / 30, 0, 0)))
  expect_equal(r[-tip_edge], 2)
  expect_equal(attr(r, "turnover"), partner_turnover(x))
  # a window of half of the branch of a halves the increase
  r_half <- turnover_rates(x, tree, k = 3, window = 0.5)
  expect_equal(r_half[tip_edge[1]], 1 + 3 * (2 / 30) * 0.5)
  # windows longer than a branch are cut to it
  expect_equal(turnover_rates(x, tree, k = 3, window = 5)[tip_edge],
               turnover_rates(x, tree, k = 3)[tip_edge])
  # the rates go to trait_rates()
  out <- trait_rates(tree, edge_rates = r)
  expect_equal(as.numeric(out$edge_rate), as.numeric(r))
})

test_that("turnover_rates() works on a simulation", {
  skip_if_not_installed("treeshapesim")
  set.seed(1)
  tr <- treeshapesim::sim_tree_shape(8)
  sim <- simulate_metacommunity(landscape_implicit(10), species_traits(tr),
                                competition_matrix(tr, "none"), n_steps = 200,
                                record_steps = seq(100, 200, 20), extirpation = 0.01)
  r <- turnover_rates(sim, k = 10)
  expect_length(r, nrow(tr$edge))
  expect_true(all(r >= 1))
})

test_that("partner_turnover() and turnover_rates() check their input", {
  x <- fake_sim()
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  expect_error(partner_turnover(list()), "mc_sim")
  y <- x
  y$history <- NULL
  expect_error(partner_turnover(y), "record_steps")
  y <- x
  y$history <- y$history[, , 1, drop = FALSE]
  expect_error(partner_turnover(y), "two recorded")
  expect_error(turnover_rates(x), "phylo")
  expect_error(turnover_rates(x, tree, sigma2 = 0), "sigma2")
  expect_error(turnover_rates(x, tree, k = -1), "`k`")
  expect_error(turnover_rates(x, tree, window = 0), "window")
  expect_error(turnover_rates(x, ape::read.tree(text = "((a:1,b:1):1,d:2);")), "differ")
})
