test_that("trait_rates multiplies branch lengths by the rates and keeps them", {
  set.seed(1)
  tree <- ape::rcoal(10)
  for (model in c("shift", "gbm", "iid")) {
    x <- trait_rates(tree, model, root = 5)
    expect_named(x, c("tips", "nodes", "edge_rates", "edge_rate", "node_rate", "shifts", "tree", "scaled"))
    expect_named(x$tips, tree$tip.label)
    expect_named(x$nodes, as.character(11:19))
    expect_length(x$edge_rates, nrow(tree$edge))
    expect_true(all(x$edge_rates > 0))
    expect_equal(x$scaled$edge.length, tree$edge.length * x$edge_rates)
    expect_null(attr(x$tips, "nodes"))
  }
  expect_null(trait_rates(tree, "iid")$shifts)
  expect_null(trait_rates(tree, "shift")$node_rate)
  # without variation among branches, every rate is sigma2
  expect_equal(trait_rates(tree, "gbm", sigma2 = 2, rate_sd = 0)$edge_rates, rep(2, 18))
  expect_equal(trait_rates(tree, "iid", sigma2 = 2, rate_sd = 0)$edge_rates, rep(2, 18))
  expect_equal(trait_rates(tree, "shift", sigma2 = 2, shifts = 0)$edge_rates, rep(2, 18))
})

test_that("shifts are inherited by descendants unless they change one branch", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  inner <- which(tree$edge[, 2] == 6)  # the branch leading to (a, b)
  below <- which(tree$edge[, 2] %in% c(1, 2))
  x <- trait_rates(tree, "shift", shift_edges = inner, shift_size = 10)
  expect_equal(x$edge_rates[c(inner, below)], rep(10, 3))
  expect_equal(x$edge_rates[-c(inner, below)], rep(1, 3))
  expect_equal(x$shifts, data.frame(edge = inner, node = 6L, position = 0, multiplier = 10,
                                    rate = 10))
  only <- trait_rates(tree, "shift", shift_edges = inner, shift_size = 10, shift_branch = TRUE)
  expect_equal(only$edge_rates[inner], 10)
  expect_equal(only$edge_rates[-inner], rep(1, 5))
  # nested shifts multiply
  nested <- trait_rates(tree, "shift", shift_edges = c(inner, below[1]), shift_size = c(10, 0.5))
  expect_equal(nested$edge_rates[below], c(5, 10))
  set.seed(2)
  drawn <- trait_rates(tree, "shift", shifts = 2)
  expect_equal(nrow(drawn$shifts), 2)
  expect_error(trait_rates(tree, "shift", shift_edges = 9), "row numbers")
  expect_error(trait_rates(tree, "shift", shift_edges = 1, shift_size = -1), "positive multiplier")
  expect_error(trait_rates(tree, "shift", shift_edges = 1, shift_size = function(n) 1:2),
               "one positive multiplier")
  expect_error(trait_rates(tree, "shift", shift_rate = -1), "zero or more")
  expect_error(trait_rates(tree, sigma2 = 0), "positive number")
  expect_error(trait_rates(tree, rate_sd = -1), "zero or more")
})

test_that("the log rate of gbm is Brownian motion and tips follow the rates", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  set.seed(3)
  logs <- replicate(4000, log(trait_rates(tree, "gbm", rate_sd = 0.5)$node_rate[["1"]]))
  # tip a is two time units from the root
  expect_equal(stats::var(logs), 0.5^2 * 2, tolerance = 0.08)
  g <- trait_rates(tree, "gbm")
  ends <- log(g$node_rate[tree$edge])
  expect_equal(log(g$edge_rates), rowMeans(matrix(ends, ncol = 2)))
  inner <- which(tree$edge[, 2] == 6)
  tips <- replicate(4000, trait_rates(tree, "shift", shift_edges = inner, shift_size = 9)$tips)
  # variance of a tip is the sum of rate times length along its path: 9 + 9 for a, 2 for c
  expect_equal(apply(tips, 1, stats::var)[c("a", "c")], c(a = 18, c = 2), tolerance = 0.08)
})

test_that("trait_mk varies rates among branches", {
  set.seed(4)
  tree <- ape::rcoal(15)
  set.seed(5)
  doubled <- trait_mk(tree, q = 1, root_state = "0", edge_rates = rep(2, 28))
  set.seed(5)
  faster <- trait_mk(tree, q = 2, root_state = "0")
  expect_equal(doubled$maps, faster$maps)
  expect_equal(doubled$edge_rates, rep(2, 28))
  expect_equal(doubled$Q, faster$Q / 2)
  set.seed(6)
  g <- trait_mk(tree, q = 1, alpha = 0.5)
  expect_length(g$edge_rates, 28)
  # edge_rates follow the edges into cladewise order
  shuffled <- ape::reorder.phylo(tree, "postorder")
  rates <- seq_len(28)
  m <- trait_mk(shuffled, edge_rates = rates)
  expect_equal(m$edge_rates, rates[match(m$edge[, 2], shuffled$edge[, 2])])
  stopped <- trait_mk(tree, q = 5, root_state = "1", edge_rates = rep(0, 28))
  expect_true(all(stopped$states == "1"))
  expect_error(trait_mk(tree, alpha = 1, edge_rates = rates), "not both")
  expect_error(trait_mk(tree, alpha = 0), "positive number")
  expect_error(trait_mk(tree, edge_rates = 1:3), "per row")
  # mean of the gamma multipliers is 1
  pair <- ape::read.tree(text = "(a:1,b:1);")
  set.seed(1)
  changes <- replicate(8000, length(trait_mk(pair, q = 0.7, root_state = "0", alpha = 0.5)$maps[[1]]) - 1)
  # Poisson with a gamma(0.5, 0.5) multiplier: mean 0.7, variance 0.7 + 0.7^2 / 0.5
  expect_equal(c(mean(changes), stats::var(changes)), c(0.7, 1.68), tolerance = 0.08)
})

test_that("hidden rate classes give an observed and a hidden history", {
  set.seed(8)
  tree <- ape::rcoal(20)
  h <- trait_mk(tree, q = 1, states = c("a", "b"), rate_classes = c(slow = 0.1, fast = 3),
                class_switch = 0.5)
  expect_s3_class(h, "simmap")
  expect_s3_class(h$hidden, "simmap")
  expect_equal(colnames(h$hidden$mapped.edge), c("a|slow", "b|slow", "a|fast", "b|fast"))
  expect_equal(colnames(h$mapped.edge), c("a", "b"))
  expect_equal(vapply(h$maps, sum, numeric(1)), h$edge.length)
  # consecutive segments of the observed history differ in state
  expect_true(all(vapply(h$maps, function(m) all(names(m)[-1] != names(m)[-length(m)]),
                         logical(1))))
  expect_equal(unname(h$states), unname(sub("\\|.*", "", h$hidden$states)))
  expect_equal(unname(h$classes), unname(sub(".*\\|", "", h$hidden$states)))
  expect_equal(unname(h$mapped.edge[, "a"]), unname(rowSums(h$hidden$mapped.edge[, c(1, 3)])))
  expect_equal(h$hidden$Q["a|slow", "b|slow"], 0.1)
  expect_equal(h$hidden$Q["a|slow", "a|fast"], 0.5)
  expect_equal(h$hidden$Q["a|slow", "b|fast"], 0)
  unnamed <- trait_mk(tree, rate_classes = c(1, 2))
  expect_true(all(unnamed$classes %in% c("R1", "R2")))
  expect_error(trait_mk(tree, rate_classes = 1), "two or more")
  expect_error(trait_mk(tree, rate_classes = c(a = 1, a = 2)), "different names")
  expect_error(trait_mk(tree, rate_classes = c(1, 2), class_switch = -1), "zero or more")
  skip_if_not_installed("phytools")
  expect_equal(sum(phytools::countSimmap(h)$N), sum(vapply(h$maps, length, integer(1)) - 1))
})

test_that("shifts can follow a Poisson process and start regimes of their own", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  set.seed(9)
  counts <- replicate(3000, nrow(trait_rates(tree, "shift", shift_rate = 0.25)$shifts))
  # total length 6: Poisson with mean 1.5
  expect_equal(c(mean(counts), stats::var(counts)), c(1.5, 1.5), tolerance = 0.08)
  set.seed(10)
  x <- trait_rates(tree, "shift", shift_rate = 1, shift_size = function(n) stats::rexp(n))
  expect_true(all(x$shifts$position >= 0 & x$shifts$position <= 1))
  expect_equal(x$scaled$edge.length, tree$edge.length * x$edge_rates)
  # a shift halfway along a branch: the rate of the branch is the mean over time
  inner <- which(tree$edge[, 2] == 6)
  below <- which(tree$edge[, 2] %in% c(1, 2))
  pair <- ape::read.tree(text = "(a:1,b:1);")
  set.seed(11)
  repeat {
    p <- trait_rates(pair, "shift", shift_rate = 1, shift_size = 3)
    once <- as.integer(names(which(table(p$shifts$edge) == 1)))
    if (length(once) > 0) break
  }
  shift <- p$shifts[p$shifts$edge == once[1], ]
  expect_equal(p$edge_rates[once[1]], shift$position + 3 * (1 - shift$position))
  # relative = FALSE: nested shifts set the rate from sigma2, not from the parent regime
  own <- trait_rates(tree, "shift", sigma2 = 2, shift_edges = c(inner, below[1]),
                     shift_size = c(10, 0.5), relative = FALSE)
  expect_equal(own$edge_rates[c(inner, below)], c(20, 1, 20))
  expect_equal(own$shifts$rate, c(20, 1)[order(c(inner, below[1]))])
  only <- trait_rates(tree, "shift", shift_edges = inner, shift_size = 10, relative = FALSE,
                      shift_branch = TRUE)
  expect_equal(only$edge_rates[below], c(1, 1))
})
