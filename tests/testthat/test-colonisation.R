churn_sim <- function(seed = 1, n_steps = 120, record = TRUE) {
  set.seed(seed)
  simulate_metacommunity(landscape_implicit(15), species_traits(8),
                         competition_matrix(8, "random"), n_steps = n_steps,
                         dispersal = 0.01, extirpation = 0.02, seed_steps = 1:5,
                         record_steps = if (record) seq_len(n_steps))
}

test_that("recording the colonisation history leaves abundances as before", {
  # sums from comdynamite 0.6.3.9000 before the history was recorded, same seeds
  set.seed(1)
  a <- simulate_metacommunity(landscape_implicit(12), species_traits(8),
                              competition_matrix(8, "none"), n_steps = 60, seed_steps = 1:5,
                              record_steps = c(10, 60))
  expect_equal(sum(a$abundance), 2011)
  set.seed(2)
  b <- simulate_metacommunity(landscape_grid(5), species_traits(6),
                              competition_matrix(6, "random"), n_steps = 80,
                              extirpation = 0.01, dispersal = 0.05)
  expect_equal(sum(b$abundance), 2387)
  set.seed(3)
  c3 <- simulate_metacommunity(landscape_grid(5), species_traits(6),
                               competition_matrix(6, "random"), n_steps = 80,
                               extirpation = 0.01, dispersal = 0.05,
                               immigration = "multinomial")
  expect_equal(sum(c3$abundance), 1701)
  set.seed(4)
  ind <- sample_individuals(a, n = 5, local_adaptation = 0.3, population_sd = 0.5,
                            displacement = 1)
  expect_equal(sum(ind$trait), 12729.16, tolerance = 1e-6)
})

test_that("sample_individuals() defaults ignore the history and draw nothing more", {
  sim <- churn_sim(record = FALSE)
  set.seed(5)
  old <- sample_individuals(sim, local_adaptation = 0.4, population_sd = 1)
  set.seed(5)
  new <- sample_individuals(sim, local_adaptation = 0.4, population_sd = 1, founder_sd = 0,
                            monopolization = 0, gene_flow = 3)
  expect_identical(old, new)
  expect_false("arrival" %in% names(old))
  # a sample after the arrival-dependent draws continues the same stream
  set.seed(5)
  mono <- sample_individuals(sim, local_adaptation = 0.4, population_sd = 1,
                             monopolization = 0.1)
  expect_equal(mono$trait - mono$population_mean, old$trait - old$population_mean)
})

test_that("arrival matches the recorded abundances, including recolonisation", {
  sim <- churn_sim()
  present <- sim$abundance > 0
  expect_identical(!is.na(sim$arrival), present)
  expect_identical(!is.na(sim$propagule), present)
  expect_true(all(sim$propagule[present] >= 1))
  expect_true(all(sim$immigrants[present] >= 0))
  expect_true(all(sim$migrant_load[present] >= 0))
  recolonised <- 0
  for (cell in which(present)) {
    row <- (cell - 1) %% nrow(present) + 1
    col <- (cell - 1) %/% nrow(present) + 1
    path <- sim$history[row, col, ]
    absent <- which(path == 0)
    expected <- if (length(absent)) max(absent) + 1 else 1
    expect_equal(sim$arrival[cell], expected)
    if (length(absent) && any(path[seq_len(min(absent))] > 0)) recolonised <- recolonised + 1
  }
  expect_gt(recolonised, 0)
  expect_equal(sim$n_steps, 120)
})

test_that("populations from `initial` arrive at step 1 with their founders", {
  set.seed(2)
  init <- matrix(c(5, 0, 3, 7, 0, 2), 3, 2, dimnames = list(NULL, c("a", "b")))
  out <- mc_dynamics(matrix(2, 3, 2, dimnames = list(NULL, c("a", "b"))),
                     diag(0.1, 2), dispersal = 0, n_steps = 30, seed_steps = numeric(0),
                     initial = init)
  present <- out$abundance > 0
  expect_true(all(out$arrival[present] == 1))
  expect_equal(out$propagule[present], init[present])
  expect_true(all(out$immigrants[present] == 0))
  expect_true(all(out$migrant_load[present] == 0))
  expect_true(all(is.na(out$arrival[!present])))
})

test_that("immigrants and seeds after arrival add to immigrants and migrant load", {
  set.seed(3)
  out <- mc_dynamics(matrix(3, 6, 2), diag(0.05, 2), dispersal = 0.2, n_steps = 40,
                     seed_steps = 1:10)
  present <- out$abundance > 0
  expect_true(any(out$immigrants[present] > 0))
  expect_true(all((out$immigrants[present] > 0) == (out$migrant_load[present] > 0)))
})

test_that("colonisation_history() lists the populations present at the last step", {
  sim <- churn_sim(record = FALSE)
  h <- colonisation_history(sim)
  expect_named(h, c("site", "species", "abundance", "arrival", "age", "propagule",
                    "immigrants", "migrant_load"))
  expect_equal(nrow(h), sum(sim$abundance > 0))
  expect_equal(h$age, sim$n_steps - h$arrival)
  expect_equal(h$abundance, sim$abundance[cbind(h$site, h$species)])
  expect_false(is.unsorted(match(h$site, rownames(sim$abundance))))
  expect_error(colonisation_history(list()), "simulate_metacommunity")
  old <- sim
  old$arrival <- NULL
  expect_error(colonisation_history(old), "no colonisation history")
})

test_that("founder deviations scale with founders and shrink with gene flow", {
  sim <- churn_sim(record = FALSE)
  draw <- function(flow) {
    set.seed(9)
    sample_individuals(sim, n = 1, intraspecific_sd = 0, founder_sd = 2, gene_flow = flow)
  }
  isolated <- draw(0)
  mixed <- draw(1)
  strong <- draw(5)
  h <- colonisation_history(sim)
  key <- paste(isolated$site, isolated$species)
  load <- h$migrant_load[match(key, paste(h$site, h$species))]
  dev0 <- isolated$population_mean - isolated$optimum
  expect_equal(mixed$population_mean - mixed$optimum, dev0 * exp(-load))
  expect_equal(strong$population_mean - strong$optimum, dev0 * exp(-5 * load))
  expect_lt(stats::sd(strong$population_mean - strong$optimum), stats::sd(dev0))
  expect_true("arrival" %in% names(isolated))
  # standardised by founders, the deviations are standard normal
  founders <- h$propagule[match(key, paste(h$site, h$species))]
  expect_gt(stats::sd(dev0 * sqrt(founders) / 2), 0.5)
  expect_lt(stats::sd(dev0 * sqrt(founders) / 2), 1.6)
})

test_that("monopolization closes the gap to the environment with the age of a population", {
  sim <- churn_sim(record = FALSE)
  h <- colonisation_history(sim)
  closed_share <- function(ind) (ind$population_mean - ind$optimum) / (ind$env - ind$optimum)
  hist_of <- function(ind) h[match(paste(ind$site, ind$species), paste(h$site, h$species)), ]
  isolated <- sample_individuals(sim, n = 1, intraspecific_sd = 0, monopolization = 0.05,
                                 gene_flow = 0)
  expect_equal(closed_share(isolated), 1 - exp(-0.05 * hist_of(isolated)$age))
  mixed <- sample_individuals(sim, n = 1, intraspecific_sd = 0, monopolization = 0.05)
  hm <- hist_of(mixed)
  rate <- 0.05 * hm$age
  expected <- ifelse(hm$age == 0, 0,
                     rate / (rate + hm$migrant_load) * (1 - exp(-rate - hm$migrant_load)))
  expect_equal(closed_share(mixed), expected)
  expect_true(all(closed_share(mixed) <= closed_share(isolated) + 1e-12))
  both <- sample_individuals(sim, n = 1, intraspecific_sd = 0, monopolization = 0.05,
                             gene_flow = 0, local_adaptation = 0.5)
  expect_equal(closed_share(both), 0.5 + 0.5 * closed_share(isolated))
})

test_that("the gap between early and late colonists grows with the rate of adaptation", {
  sim <- churn_sim(seed = 4, n_steps = 300, record = FALSE)
  gap <- function(rate) {
    ind <- sample_individuals(sim, n = 1, intraspecific_sd = 0, monopolization = rate)
    closed <- (ind$population_mean - ind$optimum) / (ind$env - ind$optimum)
    early <- ind$arrival <= stats::median(ind$arrival)
    mean(closed[early]) - mean(closed[!early])
  }
  gaps <- vapply(c(0.001, 0.003, 0.01), gap, numeric(1))
  expect_true(all(diff(gaps) > 0))
  expect_gt(gaps[1], 0)
})

test_that("per-species values and validation of the arrival-dependent arguments", {
  sim <- churn_sim(record = FALSE)
  sp <- colnames(sim$abundance)
  rates <- stats::setNames(rep(c(0, 0.1), 4), sp)
  ind <- sample_individuals(sim, n = 1, intraspecific_sd = 0, monopolization = rates)
  fixed <- ind$species %in% sp[rates == 0]
  expect_equal(ind$population_mean[fixed], ind$optimum[fixed])
  expect_error(sample_individuals(sim, founder_sd = -1), "must not be negative")
  expect_error(sample_individuals(sim, monopolization = -0.1), "must not be negative")
  expect_error(sample_individuals(sim, monopolization = 0.1, gene_flow = -1),
               "must not be negative")
  expect_error(sample_individuals(sim, founder_sd = "a"), "must be numbers")
  expect_error(sample_individuals(sim, monopolization = c(0.1, 0.2)), "one per species")
  recorded <- churn_sim(n_steps = 30)
  expect_error(sample_individuals(recorded, step = 10, founder_sd = 1), "last step")
  expect_silent(sample_individuals(recorded, step = 30, founder_sd = 1))
  old <- sim
  old$arrival <- NULL
  expect_error(sample_individuals(old, monopolization = 0.1), "no colonisation history")
  expect_silent(sample_individuals(old))
})
