test_that("river_network() grows a tree draining to one outlet", {
  set.seed(1)
  net <- river_network(40)
  expect_equal(nrow(net), 40)
  expect_true(is.na(net$downstream[1]))
  expect_true(all(net$downstream[-1] < seq_len(40)[-1]))
  expect_true(all(table(net$downstream) <= 2))
  expect_equal(net$magnitude[1], sum(net$headwater))
  expect_equal(net$distance[-1], net$distance[net$downstream[-1]] + net$length[-1])
  # distances along the channel are never shorter than straight lines
  expect_true(all(sqrt(net$x^2 + net$y^2) <= net$distance + 1e-9))
  chain <- river_network(10, branching = 0, fixed_length = TRUE)
  expect_equal(chain$distance, 0:9)
  expect_equal(sum(chain$headwater), 1)
  expect_true(all(chain$order == 1))
  set.seed(2)
  bushy <- river_network(60, branching = 1, max_branches = 3)
  expect_true(any(table(bushy$downstream) == 3))
  expect_gt(sum(bushy$headwater), sum(river_network(60, branching = 0.1)$headwater))
})

test_that("Strahler order and magnitude follow their definitions", {
  # outlet 1, confluence at 2 of two branches (3 and 4), 4 a confluence of 5 and 6
  net <- data.frame(downstream = c(NA, 1, 2, 2, 4, 4), length = c(0, 1, 1, 1, 1, 1))
  tab <- comdynamite:::river_table(net$downstream, net$length)
  expect_equal(tab$magnitude, c(3, 3, 1, 2, 1, 1))
  expect_equal(tab$order, c(2L, 2L, 1L, 2L, 1L, 1L))
  expect_equal(tab$distance, c(0, 1, 2, 2, 3, 3))
})

test_that("landscape_river() splits distances into downstream and upstream parts", {
  net <- data.frame(downstream = c(NA, 1, 2, 2, 4, 4), length = c(0, 1, 2, 1, 1, 3))
  land <- landscape_river(net, scale_down = 1, scale_up = 1)
  d <- land$downstream
  # from 3 to 5: 2 down to site 2, then 1 + 1 up
  expect_equal(d[3, 5], 2)
  expect_equal(d[5, 3], 2)
  expect_equal(land$dist[3, 5], 4)
  expect_equal(d[1, 6], 0)
  expect_equal(d[6, 1], 5)
  expect_true(isSymmetric(land$dist))
  expect_equal(unname(rowSums(land$dispersal)), rep(1, 6))
  # with equal scales the kernel depends only on the distance along the channel
  k <- exp(-land$dist[1, -1])
  expect_equal(land$dispersal[1, -1], k / sum(k))
  # drift: from the middle site, more emigrants go down than up
  drift <- landscape_river(net, scale_down = 2, scale_up = 0.5)
  expect_gt(drift$dispersal[2, 1], drift$dispersal[2, 4])
  expect_equal(drift$env, net$length * 0 + c(0, 20, 60, 40, 60, 100))
  expect_s3_class(drift, "mc_landscape")
  expect_output(print(drift), "river network \\(3 headwaters\\)")
})

test_that("landscape_river() takes coordinates, env and overland dispersal", {
  net <- data.frame(downstream = c(NA, 1, 1), length = c(0, 1, 1), x = c(0, -1, 1),
                    y = c(0, 1, 1))
  land <- landscape_river(net, env = c(5, 6, 7), overland = 0.5, overland_scale = 1)
  expect_equal(land$env, c(5, 6, 7))
  expect_equal(land$coords$x, c(0, -1, 1))
  straight <- dispersal_kernel(as.matrix(dist(land$coords)), 1)
  channel <- landscape_river(net)$dispersal
  expect_equal(land$dispersal, unname(0.5 * channel + 0.5 * straight))
  set.seed(3)
  sim <- simulate_metacommunity(landscape_river(river_network(20)),
                                data.frame(optimum = 50, sigma = 30, r_max = 3,
                                           row.names = "sp1"),
                                matrix(0.05, dimnames = list("sp1", "sp1")), n_steps = 20)
  expect_equal(nrow(sim$abundance), 20)
})

test_that("bad networks give errors", {
  expect_error(landscape_river(data.frame(downstream = c(NA, NA), length = 1)), "one outlet")
  expect_error(landscape_river(data.frame(downstream = c(NA, 3, 2), length = 1)), "loops")
  expect_error(landscape_river(data.frame(downstream = c(NA, 2), length = 1)), "another site")
  expect_error(landscape_river(data.frame(downstream = c(NA, 1), length = c(0, -1))),
               "positive")
  expect_error(landscape_river(data.frame(downstream = NA, length = 0)), "at least 2")
})
