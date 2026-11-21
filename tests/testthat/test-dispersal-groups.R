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
  # groups need one route for every species, which immigrant selection does not have
  expect_error(simulate_metacommunity(three_sites(), two_species(), diag(0.01, 2),
                                      dispersal = 0.2, n_steps = 1, initial = start,
                                      seed_steps = numeric(0), kernels = list(near, far),
                                      kernel_step = cbind(1, 2), kernel_group = 1:2,
                                      evolution = evolution_control(10,
                                                                    immigrant_selection = TRUE)),
               "immigrant selection")
})
