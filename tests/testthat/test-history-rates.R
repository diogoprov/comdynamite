# a history with every branch in state "x", then changes set by hand
fixed_history <- function(text) {
  tree <- ape::read.tree(text = text)
  trait_mk(tree, q = 0, states = c("x", "y"), root = "x")
}
set_map <- function(h, e, m) {
  h$maps[[e]] <- m
  h
}

test_that("state rates are means over time of the rate of each state", {
  h <- fixed_history("((a:1,b:1):1,c:2);")
  expect_equal(history_rates(h, "state", state_rates = c(x = 2, y = 5)), rep(2, 4))
  e <- which(h$edge[, 2] == 1)
  h <- set_map(h, e, c(x = 0.25, y = 0.75))
  rates <- history_rates(h, "state", state_rates = c(x = 2, y = 5))
  expect_equal(rates[e], 0.25 * 2 + 0.75 * 5)
  expect_equal(rates[-e], rep(2, 3))
  expect_error(history_rates(h, "state"), "named by each state")
  expect_error(history_rates(h, "state", state_rates = c(x = 1, y = -1)), "positive")
})

test_that("a burst follows each change and decays along descending branches", {
  h <- fixed_history("((a:1,b:1):1,c:2);")
  expect_equal(history_rates(h, "burst", sigma2 = 3), rep(3, 4))
  inner <- which(h$edge[, 2] == 5)
  h <- set_map(h, inner, c(x = 0.5, y = 0.5))
  rates <- history_rates(h, "burst", burst = 4, decay = 0.2)
  excess <- function(u, span) 4 * 0.2 * (exp(-u / 0.2) - exp(-(u + span) / 0.2))
  expect_equal(rates[inner], 0.5 + 0.5 + excess(0, 0.5))
  # the tips below start 0.5 after the change
  below <- which(h$edge[, 2] %in% c(1, 2))
  expect_equal(rates[below], rep(1 + excess(0.5, 1), 2))
  expect_equal(rates[h$edge[, 2] == 3], 1)
  # a later change resets the burst
  h <- set_map(h, below[1], c(y = 0.5, x = 0.5))
  expect_equal(history_rates(h, "burst", burst = 4, decay = 0.2)[below[1]],
               1 + excess(0.5, 0.5) + excess(0, 0.5))
  expect_error(history_rates(h, "burst", decay = 0), "positive")
})

test_that("sympatry counts the lineages that share a state at each moment", {
  h <- fixed_history("(a:1,b:1);")
  expect_equal(history_rates(h, "sympatry", r = -0.5), rep(exp(-1), 2))
  h <- set_map(h, 2, c(y = 1))
  expect_equal(history_rates(h, "sympatry", r = -0.5, sigma2 = 2), rep(2 * exp(-0.5), 2))
  # half the time together, half apart
  h <- set_map(h, 2, c(x = 0.5, y = 0.5))
  expect_equal(history_rates(h, "sympatry", r = -0.5),
               rep(0.5 * exp(-1) + 0.5 * exp(-0.5), 2))
  three <- fixed_history("((a:1,b:1):1,c:2);")
  rates <- history_rates(three, "sympatry", r = -0.5)
  # c is with one lineage, then with two
  expect_equal(rates[three$edge[, 2] == 3], (exp(-1) + exp(-1.5)) / 2)
  expect_error(history_rates(three, "sympatry", r = NA), "one number")
})

test_that("history rates feed trait_rates and are checked", {
  set.seed(1)
  h <- trait_mk(ape::rcoal(12), q = 1, states = c("x", "y"))
  rates <- history_rates(h, "state", state_rates = c(x = 1, y = 4))
  x <- trait_rates(h, edge_rates = rates)
  expect_equal(x$edge_rate, rates)
  expect_equal(x$scaled$edge.length, h$edge.length * rates)
  expect_null(x$shifts)
  expect_error(trait_rates(h, edge_rates = 1:3), "one positive number")
  expect_error(history_rates(ape::rcoal(4)), "simmap")
  expect_error(history_rates(h, sigma2 = 0), "positive number")
  zero <- fixed_history("((a:0,b:1):1,c:2);")
  expect_true(all(is.finite(history_rates(zero, "sympatry"))))
  expect_equal(history_rates(zero, "state", state_rates = c(x = 3, y = 1))[1:2], c(3, 3))
})
