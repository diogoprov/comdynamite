# argument checks of the simulators, null models and food webs

test_that("food web plots check their scalar arguments", {
  skip_if_not_installed("ggplot2")
  x <- food_web_dynamics(seq(0, 100, by = 10), dt = 1)
  expect_error(plot(x, time = NA), "`time`")
  expect_error(plot(x, time = "a"), "`time`")
  expect_error(plot(x, log_abundance = NA), "`log_abundance`")
  expect_error(plot(x, log_abundance = c(TRUE, FALSE)), "`log_abundance`")
  expect_error(plot(x, linewidth = 0), "`linewidth`")
  expect_error(plot(x, linewidth = Inf), "`linewidth`")
  skip_if_not_installed("gifski")
  expect_error(animate_food_web(x, file = 1), "`file`")
  expect_error(animate_food_web(x, frames = 2.5), "`frames`")
  expect_error(animate_food_web(x, frames = NA), "`frames`")
  expect_error(animate_food_web(x, frames = 3, delay = "a"), "`delay`")
  expect_error(animate_food_web(x, frames = 3, width = 0), "`width`")
  expect_error(animate_food_web(x, frames = 3, height = c(1, 2)), "`height`")
})

test_that("mc_dynamics() checks its arguments", {
  r <- matrix(2, 3, 2)
  a <- diag(0.05, 2)
  expect_error(mc_dynamics(r, a, 0, n_steps = 0), "`n_steps`")
  expect_error(mc_dynamics(r, a, 0, n_steps = NA), "`n_steps`")
  expect_error(mc_dynamics(r, a, 0, n_steps = 2, seed_mean = -1), "`seed_mean`")
  expect_error(mc_dynamics(r, a, 0, n_steps = 2, seed_mean = Inf), "`seed_mean`")
  expect_error(mc_dynamics(r, a, 0, n_steps = 2, seed_steps = NA), "`seed_steps`")
  expect_error(mc_dynamics(r, a, 0, n_steps = 2, record_steps = "a"), "`record_steps`")
  expect_error(mc_dynamics(r, a, 1.5, n_steps = 2), "`dispersal`")
  expect_error(mc_dynamics(r, a, c(0.1, 0.1, 0.1), n_steps = 2), "`dispersal`")
  expect_error(mc_dynamics(r, a, NA, n_steps = 2), "`dispersal`")
  expect_error(mc_dynamics(r, a, 0, n_steps = 2, extirpation = -0.1), "`extirpation`")
  expect_error(mc_dynamics(r, diag(3), 0, n_steps = 2), "`alpha`")
  expect_error(mc_dynamics(r, a, 0, n_steps = 2, disp_matrix = diag(2)), "`disp_matrix`")
  expect_error(mc_dynamics("a", a, 0, n_steps = 2), "`r`")
})

test_that("simulate_metacommunity() checks its arguments", {
  tr <- data.frame(optimum = c(0.2, 0.8), sigma = 0.3, r_max = 3, row.names = c("a", "b"))
  ls <- landscape_implicit(4)
  expect_error(simulate_metacommunity(list(), tr, diag(0.05, 2)), "`landscape`")
  expect_error(simulate_metacommunity(ls, tr, diag(0.05, 2), n_steps = -1), "`n_steps`")
  expect_error(simulate_metacommunity(ls, tr, diag(0.05, 2), n_steps = 5,
                                      env_time = matrix(0.5, 3, 5)), "`env_time`")
})

test_that("static simulations check their arguments", {
  a <- data.frame(optimum = 1, sigma = 1, r_max = 1, row.names = "a")
  expect_error(simulate_static(a, 1:3, incidence = NA), "`incidence`")
  expect_error(simulate_static(a, 1:3, incidence = TRUE, sampling = "expected"), "`incidence")
  expect_error(simulate_static(a, c(1, NA)), "`env`")
  expect_error(simulate_static(a, c(1, Inf)), "`env`")
  expect_error(simulate_static(a, "x"), "`env`")
  expect_error(simulate_static(list(a, a), list(1:3)), "`traits`")
  expect_error(simulate_static(list(a, a), list(1:3, 1:2)), "`env`")
  expect_error(simulate_static(transform(a, r_max = -1), 1:3), "`traits`")
  tree <- ape::compute.brlen(ape::stree(4, "balanced"))
  g <- 1:5
  expect_error(simulate_logistic(tree, g, intercept_mean = NA), "`intercept_mean`")
  expect_error(simulate_logistic(tree, g, intercept_sd = -1), "`intercept_sd`")
  expect_error(simulate_logistic(tree, g, slope_sd = c(1, 2)), "`slope_sd`")
  expect_error(simulate_logistic(tree, g, lambda = 2), "`lambda`")
  expect_error(simulate_logistic(tree, g, reduction = 0.5), "`reduction`")
  expect_error(simulate_logistic(tree, g, reduction = c(1, 0.5)), "`reduction`")
  expect_error(simulate_logistic(tree, c(1, NA)), "`gradients`")
  expect_error(simulate_logistic(tree, list(1:3, 1:4)), "`gradients`")
  expect_error(logistic_gradients(8), "`scenario`")
  expect_error(logistic_gradients(c(1, 2)), "`scenario`")
  expect_error(logistic_gradients(1, n_sites = "a"), "`n_sites`")
})

test_that("competition_matrix() checks its arguments", {
  expect_error(competition_matrix(3, inter_max = -1), "`inter_max`")
  expect_error(competition_matrix(3, inter_max = NA), "`inter_max`")
  expect_error(competition_matrix(3, scale = 0), "`scale`")
  expect_error(competition_matrix(3, scale = Inf), "`scale`")
  expect_error(competition_matrix(2.5), "`species`")
  expect_error(competition_matrix("a"), "`species`")
})

test_that("assemble_community() checks its arguments", {
  x <- c(a = 1, b = 2, c = 4, d = 7)
  expect_error(assemble_community(x, 5), "`size`")
  expect_error(assemble_community(x, 2.5), "`size`")
  expect_error(assemble_community(x, NA), "`size`")
  expect_error(assemble_community(x, 2, min_distance = -1), "`min_distance`")
  expect_error(assemble_community(x, 2, min_distance = c(1, 2)), "`min_distance`")
  expect_error(assemble_community(x, 2, "filtering", optimum = NA_real_), "`optimum`")
  expect_error(assemble_community(x, 2, "filtering", optimum = TRUE), "`optimum`")
  expect_error(assemble_community(c(x, e = Inf), 2), "`traits`")
})

test_that("patch_permutations() checks its scalar arguments", {
  coords <- expand.grid(x = 1:3, y = 1:3)
  event <- cbind(a = c(1, 1, 0, 0, 0, 0, 0, -1, -1))
  expect_error(patch_permutations(coords, event, n_perm = 0), "`n_perm`")
  expect_error(patch_permutations(coords, event, n_perm = "a"), "`n_perm`")
  expect_error(patch_permutations(coords, event, gap = 0), "`gap`")
  expect_error(patch_permutations(coords, event, gap = NA), "`gap`")
  coords$x[1] <- Inf
  expect_error(patch_permutations(coords, event), "`coords`")
})

test_that("observe_detection() checks its arguments", {
  sim <- structure(list(abundance = matrix(1:4, 2, dimnames = list(c("s1", "s2"),
                                                                    c("a", "b")))),
                   class = "mc_sim")
  expect_error(observe_detection(sim, 2), "`detection`")
  expect_error(observe_detection(sim, NA), "`detection`")
  expect_error(observe_detection(sim, "a"), "`detection`")
  expect_error(observe_detection(sim, 0.5, visits = 0), "`visits`")
  expect_error(observe_detection(sim, 0.5, visits = c(1, 2)), "`visits`")
})

test_that("env_surrogates() checks its scalar arguments", {
  coords <- as.matrix(expand.grid(1:4, 1:4))
  env <- coords[, 1] + coords[, 2]
  expect_error(env_surrogates(env, coords, n = 0), "`n`")
  expect_error(env_surrogates(env, coords, n = 1.5), "`n`")
  expect_error(env_surrogates(env, coords, n = 2, lonlat = NA), "`lonlat`")
  expect_error(env_surrogates(c(env[-1], Inf), coords, n = 2), "`env`")
  expect_error(env_surrogates(env[1:4], coords[1:4, ], n = 2), "at least 5")
  expect_error(env_surrogates(rep(1 + 1e-17, 16), coords, n = 2), "no variation")
})

test_that("persistence() checks `window`", {
  x <- array(rbinom(24, 1, 0.5), c(3, 2, 4))
  expect_error(persistence(x, window = 0), "`window`")
  expect_error(persistence(x, window = c(1, 2)), "`window`")
  expect_error(persistence(x, window = NA), "`window`")
  expect_equal(nrow(persistence(x, window = 2)), 2)
})
