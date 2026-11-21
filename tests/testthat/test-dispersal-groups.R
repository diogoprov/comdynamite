two_species <- function() {
  data.frame(optimum = c(50, 50), sigma = 1000, r_max = 2, row.names = c("sp1", "sp2"))
}
three_sites <- function() landscape_points(cbind(0:2, 0), env = rep(50, 3), disp_scale = 1)
# site 1 sends only to site 2 in "near", only to site 3 in "far"
near <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
far <- rbind(c(0, 0, 1), c(1, 0, 1), c(1, 0, 0))
start <- matrix(c(300L, 0L, 0L), 3, 2, dimnames = list(NULL, c("sp1", "sp2")))

test_that("the defaults of the new arguments reproduce earlier runs to the last bit", {
  set.seed(1)
  net <- river_network(15)
  land <- landscape_river(net)
  tr <- data.frame(optimum = c(20, 50, 80, 50), sigma = 30, r_max = 4,
                   row.names = paste0("sp", 1:4))
  alpha <- matrix(0.05, 4, 4, dimnames = list(rownames(tr), rownames(tr)))
  run <- function(..., immigration = "poisson") {
    set.seed(2)
    simulate_metacommunity(land, tr, alpha, dispersal = 0.1, n_steps = 120, disturbance = 0.05,
                           immigration = immigration, ...)$abundance
  }
  base <- run()
  expect_identical(run(disturbance_sensitivity = 1), base)
  expect_identical(run(dispersal_site = matrix(1, 15, 120)), base)
  same <- run(kernels = list(land$dispersal), kernel_step = matrix(1, 120, 2),
              kernel_group = c(1, 1, 2, 2))
  expect_identical(same, base)
  multi <- run(kernels = list(land$dispersal), kernel_step = rep(1, 120),
               immigration = "multinomial")
  multi_groups <- run(kernels = list(land$dispersal), kernel_step = matrix(1, 120, 2),
                      kernel_group = c(1, 1, 2, 2), immigration = "multinomial")
  expect_identical(multi_groups, multi)
})

test_that("kernel_group sends each group of species along its own routes", {
  for (immigration in c("poisson", "multinomial")) {
    set.seed(3)
    sim <- simulate_metacommunity(three_sites(), two_species(), diag(0.01, 2), dispersal = 0.2,
                                  n_steps = 1, initial = start, seed_steps = numeric(0),
                                  kernels = list(near = near, far = far),
                                  kernel_step = cbind(walkers = "near", fliers = "far"),
                                  kernel_group = c("walkers", "fliers"),
                                  immigration = immigration)
    expect_gt(sim$abundance[2, "sp1"], 0)
    expect_equal(sim$abundance[3, "sp1"], 0)
    expect_equal(sim$abundance[2, "sp2"], 0)
    expect_gt(sim$abundance[3, "sp2"], 0)
    expect_equal(sim$kernel_group, c("walkers", "fliers"))
  }
  # groups by index, and a change of routes through time
  set.seed(4)
  steps <- cbind(c(1, 2), c(2, 2))
  sim <- simulate_metacommunity(three_sites(), two_species(), diag(0.01, 2), dispersal = 0.2,
                                n_steps = 2, initial = start, seed_steps = numeric(0),
                                kernels = list(near, far), kernel_step = steps,
                                kernel_group = c(1, 2), record_steps = 1)
  expect_equal(sim$history[3, "sp1", 1], 0)
  expect_equal(sim$history[2, "sp2", 1], 0)
})

test_that("kernel_group works with pool_from and univariate evolution", {
  set.seed(5)
  open <- matrix(c(1, 1, 1, 1, 1, 50), 3, 2, dimnames = list(NULL, c("sp1", "sp2")))
  sim <- simulate_metacommunity(three_sites(), two_species(), diag(0.01, 2), dispersal = 0.2,
                                n_steps = 3, initial = start, seed_steps = numeric(0),
                                kernels = list(near, far), kernel_step = cbind(rep(1, 3),
                                                                               rep(2, 3)),
                                kernel_group = 1:2, pool_from = open)
  # sp2 flies only to site 3 from site 1, which is closed to it until step 50
  expect_equal(sim$abundance[2, "sp2"], 0)
  expect_equal(sim$abundance[3, "sp2"], 0)
  set.seed(5)
  multi <- simulate_metacommunity(three_sites(), two_species(), diag(0.01, 2), dispersal = 0.2,
                                  n_steps = 3, initial = start, seed_steps = numeric(0),
                                  kernels = list(near, far),
                                  kernel_step = cbind(rep(1, 3), rep(2, 3)),
                                  kernel_group = 1:2, pool_from = open,
                                  immigration = "multinomial")
  expect_equal(multi$abundance[3, "sp2"], 0)
  expect_gt(multi$abundance[2, "sp1"], 0)
  set.seed(6)
  ev <- simulate_metacommunity(three_sites(), two_species(), diag(0.01, 2), dispersal = 0.2,
                               n_steps = 1, initial = start, seed_steps = numeric(0),
                               kernels = list(near, far),
                               kernel_step = cbind(1, 2), kernel_group = 1:2,
                               evolution = evolution_control(10))
  expect_true(is.finite(ev$evolution$genetic[3, "sp2"]))
  expect_true(is.na(ev$evolution$genetic[3, "sp1"]))
})

test_that("kernel_group and kernel_step are checked", {
  run <- function(...) {
    simulate_metacommunity(three_sites(), two_species(), diag(0.01, 2), dispersal = 0.2,
                           n_steps = 2, initial = start, seed_steps = numeric(0), ...)
  }
  expect_error(run(kernels = list(near), kernel_step = matrix(1, 2, 2)), "goes with")
  expect_error(run(kernels = list(near), kernel_step = rep(1, 2), kernel_group = 1:2),
               "goes with")
  expect_error(run(kernel_group = 1:2), "needs `kernels`")
  expect_error(run(kernels = list(near), kernel_step = matrix(1, 2, 2), kernel_group = 1),
               "one group per species")
  expect_error(run(kernels = list(near), kernel_step = matrix(1, 2, 2), kernel_group = c(1, 3)),
               "column")
  expect_error(run(kernels = list(near), kernel_step = matrix(1, 2, 2),
                   kernel_group = c("a", "b")), "column names")
  expect_error(run(kernels = list(near), kernel_step = matrix(1, 1, 2), kernel_group = 1:2),
               "per step")
})

test_that("dispersal_site changes emigration site by site", {
  init <- matrix(c(300L, 300L, 0L), 3, 1, dimnames = list(NULL, "sp1"))
  tr <- two_species()[1, , drop = FALSE]
  # site 1 sends only to site 2 and site 2 only to site 3; site 2 does not emigrate
  route <- rbind(c(0, 1, 0), c(0, 0, 1), c(1, 0, 0))
  set.seed(7)
  sim <- simulate_metacommunity(three_sites(), tr, matrix(0.01), dispersal = 0.3, n_steps = 1,
                                initial = init, seed_steps = numeric(0),
                                kernels = list(route), kernel_step = 1,
                                dispersal_site = matrix(c(1, 0, 1), 3, 1))
  expect_equal(sim$abundance[3, 1], 0)
  expect_equal(sim$dispersal_site, matrix(c(1, 0, 1), 3, 1))
  # multipliers above 1 raise emigration, capped at 1
  set.seed(8)
  all_out <- simulate_metacommunity(three_sites(), tr, matrix(0.01), dispersal = 0.3,
                                    n_steps = 1, initial = init, seed_steps = numeric(0),
                                    kernels = list(route), kernel_step = 1,
                                    dispersal_site = matrix(c(10, 0, 0), 3, 1),
                                    immigration = "multinomial")
  expect_equal(all_out$abundance[1, 1], 0)
  expect_error(simulate_metacommunity(three_sites(), tr, matrix(0.01), n_steps = 2,
                                      dispersal_site = matrix(1, 2, 2)), "one row per site")
  expect_error(simulate_metacommunity(three_sites(), tr, matrix(0.01), n_steps = 2,
                                      dispersal_site = matrix(-1, 3, 2)), "multipliers")
})

test_that("disturbance_sensitivity lets resistant species keep their populations", {
  init <- matrix(200L, 3, 2, dimnames = list(NULL, c("sp1", "sp2")))
  set.seed(9)
  sim <- simulate_metacommunity(three_sites(), two_species(), diag(0.01, 2), dispersal = 0,
                                n_steps = 1, initial = init, seed_steps = numeric(0),
                                disturbance = 1, disturbance_sensitivity = c(1, 0))
  expect_true(all(sim$abundance[, "sp1"] == 0))
  expect_true(all(sim$abundance[, "sp2"] > 0))
  expect_equal(sim$disturbance_sensitivity, c(1, 0))
  # only disturbed sites lose populations
  dist <- matrix(c(1, 0, 0), 3, 1)
  set.seed(10)
  part <- simulate_metacommunity(three_sites(), two_species(), diag(0.01, 2), dispersal = 0,
                                 n_steps = 1, initial = init, seed_steps = numeric(0),
                                 disturbance = dist, disturbance_sensitivity = 0.5)
  expect_true(all(part$abundance[2:3, ] > 0))
  run <- function(s) {
    simulate_metacommunity(three_sites(), two_species(), diag(0.01, 2), n_steps = 2,
                           disturbance = 0.1, disturbance_sensitivity = s)
  }
  expect_error(run(2), "between 0 and 1")
  expect_error(run(c(0.1, 0.2, 0.3)), "one per species")
  expect_error(run(NA_real_), "between 0 and 1")
})

test_that("the evolving engine takes dispersal_site and disturbance_sensitivity", {
  set.seed(11)
  ls <- landscape_grid(5, disp_scale = 0.5)
  anc <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
  run <- function(...) {
    set.seed(12)
    simulate_evolving_metacommunity(ls, anc, speciation_control("point", rate = 0.002),
                                    dispersal = 0.05, n_steps = 80, disturbance = 0.02, ...)
  }
  base <- run()
  same <- run(disturbance_sensitivity = 1, dispersal_site = matrix(1, 25, 80))
  expect_identical(same$abundance, base$abundance)
  tough <- run(disturbance_sensitivity = 0.2)
  expect_gt(nrow(tough$events), 0)
  expect_length(tough$disturbance_sensitivity, ncol(tough$abundance))
  expect_true(all(tough$disturbance_sensitivity == 0.2))
})

test_that("the remaining checks of the new arguments", {
  # groups as a factor
  set.seed(13)
  sim <- simulate_metacommunity(three_sites(), two_species(), diag(0.01, 2), dispersal = 0.2,
                                n_steps = 1, initial = start, seed_steps = numeric(0),
                                kernels = list(near = near, far = far),
                                kernel_step = cbind(walkers = "near", fliers = "far"),
                                kernel_group = factor(c("walkers", "fliers")))
  expect_equal(sim$abundance[2, "sp2"], 0)
  # one site cannot take a dispersal pulse
  expect_error(mc_dynamics(matrix(5, 1, 1), matrix(0.05), dispersal = 0, n_steps = 2,
                           dispersal_site = matrix(1, 1, 2)), "two sites")
})

test_that("kernel_group with immigrant selection and several traits: same matrix, same run", {
  set.seed(14)
  ls <- landscape_grid(4, disp_scale = 0.7)
  other <- dispersal_kernel(as.matrix(stats::dist(expand.grid(1:4, 1:4))), 2)
  tr <- data.frame(optimum = c(30, 50, 70), sigma = 25, r_max = 4,
                   row.names = paste0("sp", 1:3))
  alpha <- matrix(0.03, 3, 3, dimnames = list(rownames(tr), rownames(tr)))
  flat <- data.frame(optimum = rep(0, 3), sigma = Inf, r_max = 1, row.names = rownames(tr))
  env <- cbind(a = ls$env, b = 0)
  steps <- rep(1:2, 20)
  run <- function(evolution, several = FALSE, groups = FALSE) {
    set.seed(15)
    sim <- simulate_metacommunity(ls, if (several) list(a = tr, b = flat) else tr, alpha,
                                  env = if (several) env, dispersal = 0.1, n_steps = 40,
                                  evolution = evolution, kernels = list(ls$dispersal, other),
                                  kernel_step = if (groups) cbind(steps, steps) else steps,
                                  kernel_group = if (groups) c(1, 2, 1))
    sim[c("abundance", "evolution", "immigrants", "migrant_load")]
  }
  G <- matrix(c(8, 2, 2, 3), 2)
  settings <- list(list(evolution_control(8, immigrant_selection = TRUE), FALSE),
                   list(evolution_control(8, mutation_var = 0.01, immigrant_selection = TRUE),
                        FALSE),
                   list(evolution_control(G), TRUE),
                   list(evolution_control(G, immigrant_selection = TRUE), TRUE))
  for (s in settings) {
    expect_identical(run(s[[1]], s[[2]], groups = TRUE), run(s[[1]], s[[2]]))
  }
})

test_that("kernel_group routes immigrants and their genes with selection and several traits", {
  flat <- data.frame(optimum = c(0, 0), sigma = Inf, r_max = 1, row.names = c("sp1", "sp2"))
  G <- diag(c(10, 1))
  settings <- list(list(evolution_control(10, immigrant_selection = TRUE), FALSE),
                   list(evolution_control(G), TRUE),
                   list(evolution_control(G, immigrant_selection = TRUE), TRUE))
  for (s in settings) {
    set.seed(16)
    sim <- simulate_metacommunity(three_sites(),
                                  if (s[[2]]) list(a = two_species(), b = flat) else
                                    two_species(),
                                  diag(0.01, 2), env = if (s[[2]]) cbind(a = rep(50, 3), b = 0),
                                  dispersal = 0.2, n_steps = 1, initial = start,
                                  seed_steps = numeric(0), kernels = list(near, far),
                                  kernel_step = cbind(1, 2), kernel_group = 1:2,
                                  evolution = s[[1]])
    expect_gt(sim$abundance[2, "sp1"], 0)
    expect_equal(sim$abundance[3, "sp1"], 0)
    expect_equal(sim$abundance[2, "sp2"], 0)
    expect_gt(sim$abundance[3, "sp2"], 0)
    z <- if (s[[2]]) sim$evolution$genetic[, , 1] else sim$evolution$genetic
    expect_true(is.finite(z[2, "sp1"]) && is.na(z[3, "sp1"]))
    expect_true(is.finite(z[3, "sp2"]) && is.na(z[2, "sp2"]))
  }
})

test_that("the evolving engine takes kernel_group: same matrix, same run", {
  set.seed(17)
  ls <- landscape_grid(5, disp_scale = 0.5)
  other <- dispersal_kernel(as.matrix(stats::dist(expand.grid(1:5, 1:5))), 2)
  tr <- species_traits(ape::rcoal(3))
  env <- cbind(a = ls$env, b = rev(ls$env))
  steps <- rep(1:2, 40)
  run <- function(speciation, several = FALSE, groups = FALSE) {
    set.seed(18)
    simulate_evolving_metacommunity(
      ls, if (several) list(tr, tr) else tr, speciation, env = if (several) env,
      evolution = evolution_control(if (several) diag(5, 2) else 5,
                                    immigrant_selection = TRUE),
      dispersal = 0.02, n_steps = 80, kernels = list(ls$dispersal, other),
      kernel_step = if (groups) cbind(a = steps, b = steps) else steps,
      kernel_group = if (groups) c("a", "b", "a"))
  }
  same <- function(x, y) {
    keep <- setdiff(names(x), c("kernel_step", "kernel_group"))
    expect_identical(x[keep], y[keep])
  }
  modes <- list(speciation_control("geographic", threshold = 30),
                speciation_control("ecological", threshold = 5, min_generations = 5),
                speciation_control("fission", rate = 0.05),
                speciation_control("point", rate = 0.002))
  for (sp in modes) same(run(sp, groups = TRUE), run(sp))
  for (sp in modes[1:3]) same(run(sp, several = TRUE, groups = TRUE), run(sp, several = TRUE))
  # without groups, no groups in the output
  expect_null(run(modes[[3]])$kernel_group)
  # new species stay in the group of their mother
  sim <- run(modes[[1]], groups = TRUE)
  expect_gt(nrow(sim$events), 0)
  expect_length(sim$kernel_group, ncol(sim$abundance))
  first <- stats::setNames(c("a", "b", "a"), sim$species$species[1:3])
  root <- function(sp) {
    while (!sp %in% names(first)) sp <- sim$species$parent[sim$species$species == sp]
    sp
  }
  expect_equal(sim$kernel_group,
               unname(first[vapply(colnames(sim$abundance), root, character(1))]))
  expect_gt(length(unique(sim$kernel_group)), 1)
})

test_that("the groups of the evolving engine follow their own routes", {
  anc <- species_traits(ape::rcoal(2), optima = c(50, 50), sigma = c(1000, 1000), r_max = c(2, 2),
                        rescale_optima = FALSE)
  init <- matrix(c(300L, 0L, 0L), 3, 2, dimnames = list(NULL, rownames(anc)))
  still <- speciation_control("fission", rate = 0)
  run <- function(speciation, ...) {
    simulate_evolving_metacommunity(three_sites(), anc, speciation, dispersal = 0.2,
                                    n_steps = 1, initial = init, competition = "none",
                                    kernels = list(near, far), kernel_step = cbind(1, 2), ...)
  }
  for (immigration in c("poisson", "multinomial")) {
    set.seed(19)
    sim <- run(still, kernel_group = 1:2, immigration = immigration)
    expect_gt(sim$abundance[2, 1], 0)
    expect_equal(sim$abundance[3, 1], 0)
    expect_equal(sim$abundance[2, 2], 0)
    expect_gt(sim$abundance[3, 2], 0)
    expect_identical(sim$kernel_group, 1:2)
  }
  expect_error(run(still, kernel_group = 1), "one group per species")
  expect_error(run(still, kernel_group = c(1, 3)), "column")
  expect_error(run(still, kernel_group = c("a", "b")), "column names")
  expect_error(simulate_evolving_metacommunity(three_sites(), anc, still, n_steps = 1,
                                               kernels = list(near), kernel_step = 1,
                                               kernel_group = 1:2), "goes with")
})

test_that("the gene flow that keeps populations together follows the kernel of each group", {
  # two alive species (the second and third ever) in groups 2 and 1
  routes <- structure(list(near, far), group = c(2L, 1L))
  st <- list2env(list(em = cbind(c(10, 0, 0), c(10, 0, 0)), k = 2, alive = 2:3,
                      kernel_group = c(1L, 2L, 1L), ev = NULL))
  zero <- matrix(0, 3, 3)
  flow <- evolving_flow(st, list(NULL, zero, zero), routes)
  expect_equal(flow[[2]], 10 * rbind(far[1, ], 0, 0))
  expect_equal(flow[[3]], 10 * rbind(near[1, ], 0, 0))
  # one kernel for all, as before
  expect_equal(evolving_flow(st, list(NULL, zero, zero), near)[[2]], 10 * rbind(near[1, ], 0, 0))
  expect_identical(species_kernel(near, 2), near)
  expect_identical(species_kernel(routes, 1), far)
  expect_null(species_kernel(NULL, 1, 1L))
  # a new population in site 3 with no recorded flow takes the divergences of the population
  # that reaches it by the kernel of the group of its species: site 1 by `one`, site 2 by `two`
  one <- rbind(c(0, 0, 1), c(1, 0, 0), c(1, 0, 0))
  two <- rbind(c(0, 1, 0), c(0, 0, 1), c(1, 0, 0))
  d <- matrix(c(0, 7, NA, 7, 0, NA, NA, NA, NA), 3)
  sp <- speciation_control("geographic", threshold = 100, divergence = "gradual")
  for (g in 1:2) {
    st <- list2env(list(k = 1, alive = 2L, n = matrix(c(100, 100, 5), 3, 1),
                        disp_matrix = structure(list(one, two), group = g),
                        kernel_group = c(3L - g, g), ev = NULL))
    upd <- evolving_divergence(st, list(NULL, d), list(NULL, zero), list(NULL, zero), sp)
    expect_equal(upd$divergence[[2]][3, 1:2], if (g == 1) c(1, 8) else c(8, 1))
  }
})
