test_that("observe_detection keeps everything with detection 1 and nothing with 0", {
  comm <- matrix(c(5, 0, 2, 10, 1, 0), 2, dimnames = list(c("s1", "s2"), c("a", "b", "c")))
  data <- mc_data(comm, series = array(c(comm, comm * 2), c(2, 3, 2),
                                       dimnames = c(dimnames(comm), list(NULL))))
  all_seen <- observe_detection(data, 1, visits = 2)
  expect_equal(all_seen$comm, comm)
  expect_equal(all_seen$series, data$series, ignore_attr = TRUE)
  expect_true(all(observe_detection(data, 0)$comm == 0))
  presence <- observe_detection(data, 1, record = "presence")
  expect_equal(presence$comm, (comm > 0) * 1)
  expect_equal(attr(presence, "detection")$record, "presence")
})

test_that("observe_detection loses rare species and more visits recover them", {
  set.seed(2)
  comm <- matrix(rep(c(1, 50), each = 2000), 2000, dimnames = list(NULL, c("rare", "common")))
  comm <- cbind(comm[, 1], comm[, 2])
  colnames(comm) <- c("rare", "common")
  data <- mc_data(comm)
  one <- observe_detection(data, 0.2, visits = 1, record = "presence")
  three <- observe_detection(data, 0.2, visits = 3, record = "presence")
  expect_equal(mean(one$comm[, "rare"]), 0.2, tolerance = 0.1)
  expect_equal(mean(three$comm[, "rare"]), 1 - 0.8^3, tolerance = 0.1)
  expect_gt(mean(one$comm[, "common"]), 0.99)
  counts <- observe_detection(data, 0.2, visits = 1)
  expect_equal(mean(counts$comm[, "common"]), 10, tolerance = 0.05)
  by_species <- observe_detection(data, c(common = 0, rare = 1), visits = 1)
  expect_equal(unname(colSums(by_species$comm)), c(2000, 0))
})

test_that("observe_detection works on simulations and checks its input", {
  skip_if_not_installed("treeshapesim")
  set.seed(3)
  tree <- treeshapesim::sim_tree_shape(8)
  sim <- simulate_metacommunity(landscape_grid(side = 3), species_traits(tree, 0.5),
                                competition_matrix(tree, "none"), dispersal = 0.05,
                                n_steps = 30, record_steps = c(20, 30))
  seen <- observe_detection(sim, 0.3)
  expect_s3_class(seen, "mc_data")
  expect_true(all(seen$comm <= sim$abundance))
  expect_true(all(seen$series <= sim$history))
  expect_type(summary_stats(seen, n_null = 9), "double")
  data <- as_mc_data(sim)
  expect_error(observe_detection(data, 2))
  expect_error(observe_detection(data, c(t1 = 0.5)), "one value per species")
  data$comm[1, 1] <- 0.5
  expect_error(observe_detection(data, 0.5), "counts")
})
