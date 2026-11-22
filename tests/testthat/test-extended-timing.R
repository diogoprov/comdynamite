# rOpenSci statistical standards tested in this file (see the srr tags in R/)
#' @srrstats {BS7.3} Tested in this file, as described in R/srr-stats-bayesian.R.
NULL

# Extended tests of how computing time scales with the size of the input (rOpenSci standard
# BS7.3), for abc_fit() and trait_gradient_bayes(). They run only with
# COMDYNAMITE_EXTENDED_TESTS=true (helper-extended.R). Timings depend on the machine, so the
# bounds are generous: they catch a change of order (for example quadratic growth where it
# should be close to linear), not small changes of speed. The scaling of the permutation tests
# is documented in data-raw/benchmark_regression.R and in the README.
#
# Set COMDYNAMITE_REPORT=true to print the observed slope and ratio.

timing_report <- function(...) {
  if (identical(Sys.getenv("COMDYNAMITE_REPORT"), "true")) message(...)
}

# shortest of `times` runs, which is less sensitive to other load on the machine than one run
fastest <- function(f, times = 3) {
  min(vapply(seq_len(times), function(i) system.time(f())[["elapsed"]], numeric(1)))
}

test_that("time of abc_fit() grows no faster than about n log n with the reference table", {
  skip_if_not_extended()
  skip_if_not_installed("abcrf")
  # A random forest grows each tree in about n log n operations for n rows, so the slope of
  # log time on log n should be close to 1. Reference tables of 500, 1000 and 2000 rows, one
  # parameter and four statistics (two informative, two noise), 50 trees, one core. The
  # slope must be below 1.6, between n log n (about 1.1 over this range) and n^2.
  set.seed(4201)
  sizes <- c(500, 1000, 2000)
  make_table <- function(n) {
    table <- data.frame(a = stats::runif(n))
    table$s1 <- table$a + stats::rnorm(n, 0, 0.05)
    table$s2 <- table$a^2 + stats::rnorm(n, 0, 0.05)
    table$s3 <- stats::rnorm(n)
    table$s4 <- stats::runif(n)
    table
  }
  tables <- lapply(sizes, make_table)
  invisible(abc_fit(tables[[1]], parameters = "a", ntree = 10, ncores = 1))
  seconds <- vapply(tables, function(table) {
    fastest(function() abc_fit(table, parameters = "a", ntree = 50, ncores = 1))
  }, numeric(1))
  slope <- unname(stats::coef(stats::lm(log(seconds) ~ log(sizes)))[2])
  timing_report("abc_fit: seconds ", paste(signif(seconds, 3), collapse = ", "),
                " for ", paste(sizes, collapse = ", "), " rows; slope ", signif(slope, 3))
  expect_true(all(is.finite(seconds) & seconds > 0))
  expect_lt(slope, 1.6)
})

test_that("time of trait_gradient_bayes() grows moderately with the number of species", {
  skip_if_not_extended()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "CmdStan is not installed")
  skip_if_not_installed("treeshapesim")
  # The likelihood uses the Cholesky factor of the phylogenetic covariance of the species, so
  # each gradient costs at most of order n^3 for n species, and doubling the species should
  # multiply the time by 8 or less (less in practice, since the populations and the sampler
  # overhead grow linearly). Fits on 15 and 30 species, 15 sites, one chain of 200 + 200
  # draws, after a first fit that compiles the model if needed. The ratio of times must be
  # below 20.
  make_data <- function(n_species, n_sites = 15) {
    tree <- treeshapesim::sim_tree_shape(n_species)
    tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
    optimum <- ape::rTraitCont(tree, sigma = 0.4)
    m <- 2 * ape::rTraitCont(tree, sigma = 1)
    env <- seq(min(optimum) - 0.1, max(optimum) + 0.1, length.out = n_sites)
    rows <- lapply(seq_len(n_species), function(i) {
      present <- which(stats::runif(n_sites) < exp(-(env - optimum[i])^2 / (2 * 0.3^2)))
      if (length(present) < 3) present <- order(abs(env - optimum[i]))[1:3]
      data.frame(site = paste0("s", present), species = tree$tip.label[i],
                 trait = m[i] + stats::rnorm(length(present), 0, 0.3))
    })
    list(tree = tree, tg = trait_gradient(do.call(rbind, rows)))
  }
  fit <- function(d) {
    suppressWarnings(suppressMessages(
      trait_gradient_bayes(d$tg, d$tree, chains = 1, iter_warmup = 200, iter_sampling = 200,
                           seed = 1, refresh = 0, show_messages = FALSE)
    ))
  }
  set.seed(4202)
  small <- make_data(15)
  large <- make_data(30)
  invisible(fit(small))
  seconds <- c(small = system.time(fit(small))[["elapsed"]],
               large = system.time(fit(large))[["elapsed"]])
  ratio <- seconds[["large"]] / seconds[["small"]]
  timing_report("trait_gradient_bayes: seconds ", paste(signif(seconds, 3), collapse = ", "),
                " for 15 and 30 species; ratio ", signif(ratio, 3))
  expect_true(all(is.finite(seconds) & seconds > 0))
  expect_lt(ratio, 20)
})
