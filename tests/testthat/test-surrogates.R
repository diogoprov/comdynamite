field_on <- function(coords, range, seed) {
  set.seed(seed)
  d <- as.matrix(dist(coords))
  unname(drop(t(chol(exp(-d / range) + diag(1e-8, nrow(d)))) %*% rnorm(nrow(d))))
}
grid <- as.matrix(expand.grid(1:12, 1:12))

test_that("parametric surrogates are permutations of whole rows", {
  x <- field_on(grid, 3, 1)
  env <- cbind(a = x, b = x + rnorm(144, sd = 0.5))
  draws <- env_surrogates(env, grid, n = 4, method = "parametric")
  rows <- attr(draws, "rows")
  expect_length(draws, 4)
  expect_equal(dim(rows), c(4L, 144L))
  expect_true(all(apply(rows, 1, function(r) all(sort(r) == 1:144))))
  expect_equal(draws[[2]], env[rows[2, ], ])
  expect_equal(colnames(draws[[1]]), c("a", "b"))
  one <- env_surrogates(x, grid, n = 1, method = "parametric")
  expect_equal(sort(one[[1]][, 1]), sort(x))
  expect_equal(dim(attr(one, "rows")), c(1L, 144L))
})

test_that("the fitted range follows the range of the field", {
  d <- as.matrix(dist(grid))
  short <- median(vapply(1:5, function(s) fit_exponential_range(field_on(grid, 1, s), d), 1))
  long <- median(vapply(1:5, function(s) fit_exponential_range(field_on(grid, 6, s), d), 1))
  expect_lt(short, long)
})

test_that("variogram surrogates keep the autocorrelation and lose the values", {
  x <- field_on(grid, 4, 2)
  draws <- env_surrogates(x, grid, n = 20)
  setup <- variogram_setup(as.matrix(dist(grid)))
  target <- binned_variogram(x, setup)
  mean_vgram <- rowMeans(vapply(draws, function(e) binned_variogram(e[, 1], setup),
                                numeric(length(target))))
  permuted <- binned_variogram(sample(x), setup)
  # at short distances the surrogates are much closer to the observed variogram than a
  # permutation, which has no autocorrelation
  expect_lt(abs(mean_vgram[1] - target[1]), abs(permuted[1] - target[1]) / 2)
  expect_lt(abs(mean(vapply(draws, function(e) cor(e[, 1], x), 1))), 0.2)
  two <- env_surrogates(cbind(p = x, q = -x), grid, n = 2)
  expect_equal(colnames(two[[1]]), c("p", "q"))
})

test_that("env_surrogates() checks its input", {
  x <- field_on(grid, 3, 1)
  expect_error(env_surrogates(x, grid[-1, ]), "nrow")
  expect_error(env_surrogates(x, rbind(grid[-1, ], grid[1, , drop = FALSE])[c(144, 1:143), ][
    c(1, 1:143), ]), "same coordinates")
  expect_error(env_surrogates(x, grid, method = "torus"), "should be one of")
})

test_that("the tests take surrogates as their null of sites", {
  set.seed(5)
  tree <- ape::rcoal(8)
  coords <- as.matrix(expand.grid(1:6, 1:6))
  env <- cbind(x = field_on(coords, 2, 3))
  rate <- exp(-outer(env[, 1], seq(-2, 2, length.out = 8), "-")^2)
  comm <- matrix(rpois(36 * 8, 4 * rate) + 1, 36,
                 dimnames = list(paste0("s", 1:36), tree$tip.label))
  out <- pcps_test(comm, tree, env, n_perm = 9, null = "variogram", coords = coords)
  expect_equal(out$permutations$null, "variogram")
  expect_output(print(out), "9 variogram surrogates of the predictors")
  beta <- phylo_beta_test(comm, tree, env, n_perm = 9, null = "parametric", coords = coords)
  expect_equal(beta$permutations$null, "parametric")
  expect_gte(beta$p_site, 1 / 10)
  nodes <- node_contrasts(tree)
  node <- node_test(comm, nodes, env, n_perm = 9, null = "parametric", coords = coords)
  expect_equal(node$permutations$null, "parametric")
  expect_true(all(node$nodes$p_site >= 1 / 10, na.rm = TRUE))
  data <- mc_data(comm, tree = tree, env = data.frame(env), coords = coords)
  from_data <- pcps_test(data, n_perm = 9, null = "parametric")
  expect_equal(from_data$permutations$null, "parametric")
  expect_error(pcps_test(comm, tree, env, n_perm = 9, null = "variogram"), "needs `coords`")
  expect_error(pcps_test(comm, tree, env, n_perm = 9, null = "variogram", coords = coords[-1, ]),
               "one row per site")
  expect_error(pcps_test(comm, tree, env, n_perm = 9, null = "parametric", coords = coords,
                         site_permutations = torus_permutations(6)), "either")
})

test_that("landscapes get surrogate environments with their own distances", {
  set.seed(7)
  land <- landscape_grid(8, env_range = c(0, 100))
  nulls <- env_surrogates(land, n = 3, method = "parametric")
  expect_length(nulls, 3)
  expect_s3_class(nulls[[1]], "mc_landscape")
  expect_equal(sort(nulls[[1]]$env), sort(land$env))
  expect_equal(nulls[[2]]$env, land$env[attr(nulls, "rows")[2, ]])
  expect_identical(nulls[[1]]$dispersal, land$dispersal)
  axes <- landscape_grid(8, n_axes = 2, axis_cor = 0.5)
  moved <- env_surrogates(axes, n = 2)
  expect_equal(dim(moved[[1]]$env_axes), dim(axes$env_axes))
  expect_equal(moved[[1]]$env, unname(moved[[1]]$env_axes[, 1]))
  points <- landscape_points(as.matrix(expand.grid(1:6, 1:6)), env = runif(36), disp_scale = 2)
  points$dist <- NULL
  expect_length(env_surrogates(points, n = 1), 1)
  expect_error(env_surrogates(landscape_implicit(10), n = 1), "neither distances")
})

test_that("an env_time matrix moves whole series with the parametric surrogate", {
  set.seed(8)
  coords <- as.matrix(expand.grid(1:6, 1:6))
  series <- env_space_time(coords, n_steps = 30, burn_in = 5, spatial_scale = 3,
                           temporal_scale = 10)
  series[, 1] <- 0.5
  draws <- env_surrogates(series, coords, n = 2, method = "parametric")
  expect_equal(draws[[1]], series[attr(draws, "rows")[1, ], ])
  smooth <- env_surrogates(series[, 1:4], coords, n = 1)
  expect_equal(smooth[[1]][, 1], rep(0.5, 36))
  expect_error(env_surrogates(matrix(1, 36, 2), coords), "no variation")
  expect_error(env_surrogates(c(1:35, NA), coords), "missing values")
  expect_error(env_surrogates(data.frame(a = letters[1:36]), coords), "dummy")
  expect_error(env_surrogates(matrix(letters[1:36]), coords), "not numeric")
})

test_that("surrogates use great-circle distances and stay within strata", {
  set.seed(9)
  lonlat <- cbind(-48 + runif(30, 0, 0.5), -21 + runif(30, 0, 0.5))
  x <- lonlat[, 1] + rnorm(30, sd = 0.05)
  expect_length(env_surrogates(x, lonlat, n = 2, lonlat = TRUE), 2)
  expect_equal(great_circle(rbind(c(0, 0), c(0, 1)))[1, 2], 6371 * pi / 180)
  city <- rep(c("a", "b", "c"), each = 10)
  env <- cbind(x = x + 10 * (city == "b"))
  env[city == "c", ] <- 3
  draws <- env_surrogates(env, lonlat, n = 3, method = "parametric", lonlat = TRUE,
                          strata = city)
  rows <- attr(draws, "rows")
  expect_true(all(city[rows[1, ]] == city))
  expect_equal(draws[[2]][, 1], env[rows[2, ], 1])
  expect_equal(draws[[1]][city == "c", 1], rep(3, 10))
  smooth <- env_surrogates(env, lonlat, n = 2, lonlat = TRUE, strata = city)
  expect_lt(abs(mean(smooth[[1]][city == "b", 1]) - mean(env[city == "b", 1])), 1)
  expect_error(env_surrogates(env, lonlat, strata = city[-1]), "one value per site")
  expect_error(env_surrogates(env, lonlat, strata = rep(1:6, each = 5)[c(1:29, 1)]), "at least 5")
})

test_that("the tests take a list of surrogates", {
  set.seed(5)
  tree <- ape::rcoal(8)
  coords <- as.matrix(expand.grid(1:6, 1:6))
  env <- cbind(x = field_on(coords, 2, 3))
  rate <- exp(-outer(env[, 1], seq(-2, 2, length.out = 8), "-")^2)
  comm <- matrix(rpois(36 * 8, 4 * rate) + 1, 36,
                 dimnames = list(paste0("s", 1:36), tree$tip.label))
  strata <- rep(1:2, each = 18)
  draws <- env_surrogates(env, coords, n = 7, method = "parametric", strata = strata)
  out <- pcps_test(comm, tree, env, n_perm = 9, null = draws)
  expect_equal(out$permutations$null, "surrogates")
  expect_equal(out$permutations$sites, 7)
  expect_output(print(out), "7 surrogates of the predictors")
  beta <- phylo_beta_test(comm, tree, env, n_perm = 9, null = draws)
  expect_equal(beta$permutations$sites, 7)
  node <- node_test(comm, node_contrasts(tree), env, n_perm = 9, null = draws)
  expect_equal(node$permutations$sites, 7)
  expect_error(pcps_test(comm, tree, env, n_perm = 9, null = list(1:3)), "numeric matrices")
  expect_error(pcps_test(comm, tree, env, n_perm = 9, null = lapply(draws, `[`, 1:30, ,
                                                                     drop = FALSE)),
               "one row per site")
  expect_error(pcps_test(comm, tree, env, n_perm = 9, null = lapply(draws, cbind, 1)),
               "one column per predictor")
  expect_error(pcps_test(comm, tree, env, n_perm = 9, null = draws,
                         site_permutations = torus_permutations(6)), "either")
})
