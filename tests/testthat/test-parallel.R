sim_small <- function() {
  set.seed(1)
  tree <- ape::rcoal(8)
  comm <- matrix(rpois(20 * 8, 2) + 1, 20, dimnames = list(paste0("s", 1:20), tree$tip.label))
  list(tree = tree, comm = comm, env = cbind(x = seq_len(20)))
}

test_that("permutations under a future plan repeat with a seed, for any number of workers", {
  skip_if_not_installed("future.apply")
  skip_on_os("windows")
  d <- sim_small()
  run <- function(workers) {
    future::plan(future::multicore, workers = workers)
    set.seed(3)
    pcps_test(d$comm, d$tree, d$env, n_perm = 19)$tests
  }
  on.exit(future::plan(future::sequential), add = TRUE)
  two <- run(2)
  expect_identical(two, run(2))
  expect_identical(two, run(3))
  future::plan(future::sequential)
  expect_equal(future_workers(), 1)
  expect_equal(future_workers(installed = FALSE), 1)
})

test_that("tibbles enter and results leave as tables", {
  skip_if_not_installed("tibble")
  d <- sim_small()
  comm <- tibble::as_tibble(d$comm)
  env <- tibble::as_tibble(d$env)
  rownames(d$comm) <- NULL
  set.seed(4)
  tidy <- pcc_test(comm, d$tree, env, n_perm = 9)
  set.seed(4)
  expect_equal(as.data.frame(tidy), as.data.frame(pcc_test(d$comm, d$tree, d$env, n_perm = 9)))
  expect_s3_class(tibble::as_tibble(as.data.frame(tidy)), "tbl_df")
  expect_s3_class(tibble::as_tibble(pcps_test(comm, d$tree, env, n_perm = 3)), "tbl_df")
  expect_s3_class(phylo_beta_test(comm, d$tree, env, n_perm = 3), "mc_phylo_beta")
})
