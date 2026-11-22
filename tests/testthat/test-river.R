test_that("river_network() grows a tree draining to one outlet", {
  set.seed(1)
  net <- river_network(40)
  expect_identical(nrow(net), 40L)
  expect_true(is.na(net$downstream[1]))
  expect_true(all(net$downstream[-1] < seq_len(40)[-1]))
  expect_true(all(table(net$downstream) <= 2))
  expect_equal(net$magnitude[1], sum(net$headwater))
  expect_identical(net$distance[-1], net$distance[net$downstream[-1]] + net$length[-1])
  # distances along the channel are never shorter than straight lines
  expect_true(all(sqrt(net$x^2 + net$y^2) <= net$distance + 1e-9))
  chain <- river_network(10, branching = 0, fixed_length = TRUE)
  expect_equal(chain$distance, 0:9)
  expect_identical(sum(chain$headwater), 1L)
  expect_true(all(chain$order == 1))
  set.seed(2)
  bushy <- river_network(60, branching = 1, max_branches = 3)
  expect_true(any(table(bushy$downstream) == 3))
  expect_gt(sum(bushy$headwater), sum(river_network(60, branching = 0.1)$headwater))
})

test_that("Strahler order and magnitude follow their definitions", {
  # outlet 1, confluence at 2 of two branches (3 and 4), 4 a confluence of 5 and 6
  net <- data.frame(downstream = c(NA, 1, 2, 2, 4, 4), length = c(0, 1, 1, 1, 1, 1))
  tab <- river_table(net$downstream, net$length)
  expect_identical(tab$magnitude, c(3, 3, 1, 2, 1, 1))
  expect_identical(tab$order, c(2L, 2L, 1L, 2L, 1L, 1L))
  expect_identical(tab$distance, c(0, 1, 2, 2, 3, 3))
  expect_identical(tab$upstream_length, c(5, 5, 1, 3, 1, 1))
})

test_that("landscape_river() splits distances into downstream and upstream parts", {
  net <- data.frame(downstream = c(NA, 1, 2, 2, 4, 4), length = c(0, 1, 2, 1, 1, 3))
  land <- landscape_river(net, scale_down = 1, scale_up = 1)
  d <- land$downstream
  # from 3 to 5: 2 down to site 2, then 1 + 1 up
  expect_identical(d[3, 5], 2)
  expect_identical(d[5, 3], 2)
  expect_identical(land$dist[3, 5], 4)
  expect_identical(d[1, 6], 0)
  expect_identical(d[6, 1], 5)
  expect_true(isSymmetric(land$dist))
  expect_equal(unname(rowSums(land$dispersal)), rep(1, 6))
  # with equal scales the kernel depends only on the distance along the channel
  k <- exp(-land$dist[1, -1])
  expect_equal(land$dispersal[1, -1], k / sum(k))
  # drift: from the middle site, more emigrants go down than up
  drift <- landscape_river(net, scale_down = 2, scale_up = 0.5)
  expect_gt(drift$dispersal[2, 1], drift$dispersal[2, 4])
  expect_identical(drift$env, net$length * 0 + c(0, 20, 60, 40, 60, 100))
  expect_s3_class(drift, "mc_landscape")
  expect_output(print(drift), "river network \\(3 headwaters\\)")
})

test_that("landscape_river() takes coordinates, env and overland dispersal", {
  net <- data.frame(downstream = c(NA, 1, 1), length = c(0, 1, 1), x = c(0, -1, 1),
                    y = c(0, 1, 1))
  land <- landscape_river(net, env = c(5, 6, 7), overland = 0.5, overland_scale = 1)
  expect_identical(land$env, c(5, 6, 7))
  expect_identical(land$coords$x, c(0, -1, 1))
  straight <- dispersal_kernel(as.matrix(dist(land$coords)), 1)
  channel <- landscape_river(net)$dispersal
  expect_equal(land$dispersal, unname(0.5 * channel + 0.5 * straight))
  set.seed(3)
  sim <- simulate_metacommunity(landscape_river(river_network(20)),
                                data.frame(optimum = 50, sigma = 30, r_max = 3,
                                           row.names = "sp1"),
                                matrix(0.05, dimnames = list("sp1", "sp1")), n_steps = 20)
  expect_identical(nrow(sim$abundance), 20L)
})

test_that("bad networks give errors", {
  expect_error(landscape_river(data.frame(downstream = c(NA, NA), length = 1)), "one outlet")
  expect_error(landscape_river(data.frame(downstream = c(NA, 3, 2), length = 1)), "loops")
  expect_error(landscape_river(data.frame(downstream = c(NA, 2), length = 1)), "another site")
  expect_error(landscape_river(data.frame(downstream = c(NA, 1), length = c(0, -1))),
               "positive")
  expect_error(landscape_river(data.frame(downstream = NA, length = 0)), "at least 2")
})

test_that("river_covariance() builds valid tail-up, tail-down and Euclidean components", {
  set.seed(4)
  for (b in c(0.2, 0.6, 1)) {
    land <- landscape_river(river_network(60, branching = b))
    for (w in c("magnitude", "upstream_length")) {
      v <- river_covariance(land, tail_up = c(2, 3), weights = w)
      expect_true(isSymmetric(v))
      expect_identical(diag(v), rep(2, 60))
      expect_gt(min(eigen(v, symmetric = TRUE, only.values = TRUE)$values), -1e-10)
      unconnected <- land$downstream > 0 & t(land$downstream) > 0
      expect_true(all(v[unconnected] == 0))
    }
    v <- river_covariance(land, tail_down = c(1, 3), euclidean = c(0.5, 2), nugget = 0.1)
    expect_gt(min(eigen(v, symmetric = TRUE, only.values = TRUE)$values), 0)
    expect_equal(diag(v), rep(1.6, 60))
  }
  # with Shreve weights the additive function is the magnitude over that of the outlet
  af <- river_additive(land$network, "magnitude")
  expect_equal(af, land$network$magnitude / land$network$magnitude[1])
  # hand example: 1 <- 2 <- (3, 4); 4 <- (5, 6)
  net <- data.frame(downstream = c(NA, 1, 2, 2, 4, 4), length = c(0, 1, 2, 1, 1, 3))
  small <- landscape_river(net)
  v <- river_covariance(small, tail_up = c(sill = 1, range = 2))
  expect_equal(v[5, 1], exp(-3 / 2) * sqrt(1 / 3))
  expect_identical(v[3, 4], 0)
  expect_equal(v[4, 2], exp(-1 / 2) * sqrt(2 / 3))
  v2 <- river_covariance(small, tail_up = c(range = 2, sill = 1), weights = c(1, 1, 1, 3, 1, 1))
  expect_equal(v2[4, 2], exp(-1 / 2) * sqrt(3 / 4))
  expect_equal(river_covariance(small, tail_down = c(1, 2))[3, 5], exp(-4 / 2))
})

test_that("env_river() draws fields and series with the covariance asked for", {
  set.seed(5)
  land <- landscape_river(river_network(15))
  v <- river_covariance(land, tail_up = c(1, 4), nugget = 0.2)
  draws <- replicate(4000, env_river(land, tail_up = c(1, 4), nugget = 0.2))
  expect_lt(max(abs(stats::cov(t(draws)) - v)), 0.12)
  e <- env_river(land, mean = 10, tail_down = c(1, 2))
  expect_length(e, 15)
  s <- env_river(land, n_steps = 3000, tail_up = c(1, 4), temporal_scale = 10)
  expect_identical(dim(s), c(15L, 3000L))
  lag1 <- mean(vapply(1:15, function(i) stats::cor(s[i, -1], s[i, -3000]), numeric(1)))
  expect_equal(lag1, exp(-1 / 10), tolerance = 0.05)
  still <- env_river(land, n_steps = 5, euclidean = c(1, 2))
  expect_true(all(still == still[, 1]))
  free <- env_river(land, n_steps = 2000, nugget = 1, temporal_scale = 0)
  expect_lt(abs(stats::cor(free[1, -1], free[1, -2000])), 0.1)
  sim <- simulate_metacommunity(land, data.frame(optimum = 0, sigma = 2, r_max = 3,
                                                 row.names = "sp1"),
                                matrix(0.05, dimnames = list("sp1", "sp1")), n_steps = 20,
                                env_time = env_river(land, n_steps = 20, tail_up = c(1, 4)))
  expect_identical(nrow(sim$abundance), 15L)
})

test_that("river_covariance() checks its arguments", {
  land <- landscape_river(data.frame(downstream = c(NA, 1, 1), length = c(0, 1, 1)))
  expect_error(river_covariance(landscape_grid(3)), "landscape_river")
  expect_error(river_covariance(land), "at least one")
  expect_error(river_covariance(land, tail_up = 1), "two positive")
  expect_error(river_covariance(land, tail_up = c(1, 1), weights = c(1, -1, 1)), "weights")
  expect_error(river_covariance(land, tail_up = c(1, 1), weights = "area"))
})
