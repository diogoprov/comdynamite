# Extracted from test-plot-rates.R:68

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "comdynamite", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("ggplot2")
set.seed(6)
tree <- ape::rtree(10)
comm <- matrix(rpois(200, 2), 20, dimnames = list(NULL, tree$tip.label))
p <- plot(pcps_test(comm, tree, seq_len(20), n_perm = 9))
expect_equal(p$labels$colour, "env1")
