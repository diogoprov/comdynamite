# rOpenSci statistical standards tested in this file (see the srr tags in R/)
#' @srrstats {RE7.1, RE7.1a} Tested in this file, as described in R/srr-stats-regression.R.
NULL

# noiseless data: the response is an exact linear function of the predictors

flex_data <- function(noise) {
  set.seed(1)
  env <- data.frame(site = paste0("s", 1:8), moisture = c(3, 1, 4, 1, 5, 9, 2, 6))
  ind <- expand.grid(species = c("a", "b", "c"), site = env$site, stringsAsFactors = FALSE)
  ind$trait <- c(a = 10, b = 20, c = 30)[ind$species] +
    env$moisture[match(ind$site, env$site)] + noise * stats::rnorm(nrow(ind))
  list(ind = ind, env = env)
}

test_that("warn_perfect_fit() flags fits without residual variation", {
  expect_silent(warn_perfect_fit(c(0.5, NA), c("x", "y"), "these"))
  expect_warning(exact <- warn_perfect_fit(c(0.5, 1 - 1e-12), c("x", "y"), "these"),
                 "these fit exactly.*: y\\. Their F statistics")
  expect_identical(exact, c(FALSE, TRUE))
})

test_that("trait_flex() warns on a noiseless model and recovers it", {
  d <- flex_data(0)
  expect_warning(res <- trait_flex(d$ind, formula = ~moisture, env = d$env),
                 "fit exactly.*trait \\(intraspecific\\), trait \\(total\\)")
  total <- res$decomposition[res$decomposition$component == "total", ]
  expect_equal(total$residual, 0, tolerance = 1e-8)
  expect_equal(total$explained, total$total)
  noisy <- flex_data(1)
  expect_no_warning(trait_flex(noisy$ind, formula = ~moisture, env = noisy$env))
})

test_that("drift_test() warns when the variances are exactly proportional", {
  set.seed(2)
  q <- 4
  values <- c(8, 4, 2, 1)
  vectors <- qr.Q(qr(matrix(stats::rnorm(q * q), q)))
  W <- vectors %*% diag(values) %*% t(vectors)
  z <- qr.Q(qr(scale(matrix(stats::rnorm(8 * q), 8), scale = FALSE)))
  means <- z %*% diag(sqrt(values)) %*% t(vectors)
  expect_warning(res <- drift_test(means, 1:8, W = W, correct = FALSE),
                 "log\\(between variance\\) and log\\(eigenvalue\\) fit exactly")
  expect_equal(res$slope, 1, tolerance = 1e-8)
  expect_equal(res$se, 0, tolerance = 1e-8)
})

test_that("a noiseless fit takes no longer than a noisy one", {
  exact <- flex_data(0)
  noisy <- flex_data(1)
  time_exact <- system.time(for (i in 1:20) {
    suppressWarnings(trait_flex(exact$ind, formula = ~moisture, env = exact$env))
  })[["elapsed"]]
  time_noisy <- system.time(for (i in 1:20) {
    trait_flex(noisy$ind, formula = ~moisture, env = noisy$env)
  })[["elapsed"]]
  expect_lt(time_exact, 3 * time_noisy + 0.5)
})

test_that("node_test() warns when the node response is an exact function of `env`", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  comm <- rbind(c(1, 0, 1), c(0, 1, 1), c(1, 1, 1), c(1, 0, 1), c(0, 1, 1), c(1, 1, 0),
                c(1, 0, 0), c(0, 1, 0))
  dimnames(comm) <- list(paste0("s", 1:8), c("a", "b", "c"))
  # a and b occur at the same number of sites, so the response of node 5 is b - a
  exact <- comm[, "b"] - comm[, "a"]
  expect_warning(out <- node_test(comm, node_contrasts(tree), exact, n_perm = 9),
                 "fit exactly.*: node 5\\.")
  expect_equal(out$nodes$r2[out$nodes$node == "5"], 1)
  set.seed(2)
  expect_no_warning(node_test(comm, node_contrasts(tree), exact + stats::rnorm(8, 0, 0.3),
                              n_perm = 9))
})

test_that("pcps_test() warns when an axis is an exact function of `env`", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  comm <- rbind(c(1, 0, 1, 0), c(0, 1, 1, 1), c(1, 1, 0, 0), c(1, 0, 0, 1), c(0, 1, 1, 0),
                c(1, 1, 1, 0), c(0, 0, 1, 1))
  dimnames(comm) <- list(paste0("s", 1:7), letters[1:4])
  axis <- pcps_test(comm, tree, 1:7, n_perm = 9)$scores[, "pcps1"]
  expect_warning(out <- pcps_test(comm, tree, axis, n_perm = 9), "fit exactly.*: pcps1\\.")
  expect_equal(out$tests$r2[out$tests$test == "pcps1"], 1)
  set.seed(2)
  expect_no_warning(pcps_test(comm, tree, axis + stats::rnorm(7, 0, 0.3), n_perm = 9))
})
