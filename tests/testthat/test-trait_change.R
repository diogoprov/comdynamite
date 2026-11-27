tc_sim <- function(g = 1, b = 0.3, rate = 0.5, record = c(20, 40, 60)) {
  set.seed(8)
  ls <- landscape_grid(4)
  tr <- ape::rcoal(4)
  env_time <- outer(ls$env, rate * seq_len(60), "+")
  simulate_metacommunity(ls, species_traits(tr), competition_matrix(tr, "phylo"),
                         n_steps = 60, record_steps = record, env_time = env_time,
                         evolution = evolution_control(g, plasticity = b))
}

test_that("the community partition matches a hand calculation", {
  d <- data.frame(site = "a", species = rep(c("x", "y", "w"), 2), step = rep(1:2, each = 3),
                  abundance = c(10, 10, 0, 5, 15, 5), trait = c(1, 3, NA, 1.5, 3, 4))
  res <- trait_change(d)
  com <- res$community
  # relative abundances 0.5, 0.5, 0 and 0.2, 0.6, 0.2
  expect_equal(com$mean_from, 2)
  expect_equal(com$mean_to, 0.2 * 1.5 + 0.6 * 3 + 0.2 * 4)
  expect_equal(com$sorting, -0.3 * 1 + 0.1 * 3)
  expect_equal(com$turnover, 0.2 * 4)
  expect_equal(com$within, 0.5 * 0.5)
  expect_equal(com$covariance, -0.3 * 0.5)
  expect_equal(com$total, com$sorting + com$turnover + com$within + com$covariance)
  expect_true(all(is.na(com[c("genetic", "plastic", "cov_genetic", "cov_plastic")])))
  expect_null(res$growth)
  expect_false(res$split)
})

test_that("a genetic column splits the change within species", {
  d <- data.frame(site = rep(c("a", "b"), each = 4), species = rep(c("x", "y"), 4),
                  step = rep(rep(1:2, each = 2), 2), abundance = c(4, 6, 3, 7, 5, 5, 8, 2),
                  trait = c(1, 2, 2, 2.5, 0, 1, 0.5, 2), genetic = c(1, 2, 1.5, 2.2, 0, 1, 0, 1.4))
  com <- trait_change(d, genetic = "genetic")$community
  expect_equal(com$within, com$genetic + com$plastic)
  expect_equal(com$covariance, com$cov_genetic + com$cov_plastic)
  expect_equal(com$total, com$sorting + com$turnover + com$within + com$covariance)
  s <- summary(trait_change(d, genetic = "genetic"))
  expect_identical(s$community$sites, 2L)
  expect_equal(s$community$total, mean(com$total))
  expect_output(print(trait_change(d, genetic = "genetic")), "Community-weighted mean")
})

test_that("empty sites give NA and are left out of the summary", {
  d <- data.frame(site = rep(c("a", "b"), each = 2), species = "x", step = rep(1:2, 2),
                  abundance = c(1, 2, 0, 3), trait = c(1, 2, NA, 2))
  res <- trait_change(d)
  expect_true(all(is.na(res$community[res$community$site == "b", -(1:3)])))
  expect_identical(summary(res)$community$sites, 1L)
})

test_that("the partition of a simulation adds up, at both levels", {
  sim <- tc_sim()
  res <- trait_change(sim)
  com <- res$community
  expect_identical(nrow(com), 2L * 16L)
  ok <- !is.na(com$total)
  expect_equal(com$total[ok], (com$sorting + com$turnover + com$within + com$covariance)[ok],
               tolerance = 1e-10)
  expect_equal(com$within[ok], (com$genetic + com$plastic)[ok], tolerance = 1e-10)
  # with all populations persisting, the plastic part is b times the change in environment
  shared <- com$turnover[ok] == 0
  expect_equal(com$plastic[ok][shared], rep(0.3 * 10, sum(shared)), tolerance = 1e-10)
  gr <- res$growth
  expect_gt(nrow(gr), 0)
  expect_equal(gr$total, gr$genetic + gr$plastic + gr$environment, tolerance = 1e-10)
  # the total is the change in log growth rate computed from the recorded state
  ev <- sim$evolution
  r <- gr[1, ]
  niche <- sim$traits[r$species, ]
  x_of <- function(s) {
    lab <- paste0("step", s)
    p <- ev$history[r$site, r$species, lab] +
      0.3 * (ev$environment_history[r$site, lab] - ev$reference)
    e <- ev$environment_history[r$site, lab]
    log(evolution_growth(matrix(p), e, niche$sigma, niche$r_max, 1))
  }
  expect_equal(r$total, c(x_of(r$to) - x_of(r$from)), tolerance = 1e-10)
  expect_s3_class(res, "trait_change")
  expect_output(print(res), "Growth rate")
  expect_identical(trait_change(sim, steps = c(20, 60))$steps, c(20, 60))
  expect_null(trait_change(sim, growth = FALSE)$growth)
})

test_that("each part vanishes when its process is off", {
  com_b0 <- trait_change(tc_sim(b = 0))
  expect_equal(com_b0$community$plastic[!is.na(com_b0$community$plastic)], rep(0, 32 -
    sum(is.na(com_b0$community$plastic))))
  expect_true(all(com_b0$growth$plastic == 0))
  com_g0 <- trait_change(tc_sim(g = 0))
  expect_true(all(abs(com_g0$community$genetic) < 1e-12, na.rm = TRUE))
  expect_true(all(abs(com_g0$growth$genetic) < 1e-12))
  static <- trait_change(tc_sim(rate = 0))
  expect_true(all(abs(static$community$plastic) < 1e-12, na.rm = TRUE))
  expect_true(all(abs(static$growth$environment) < 1e-12))
})

test_that("several traits are taken one at a time, without growth", {
  set.seed(8)
  ls <- landscape_grid(4)
  tr <- ape::rcoal(4)
  traits <- species_traits(tr)
  a <- competition_matrix(tr, "phylo")
  two <- simulate_metacommunity(ls, list(traits, traits), a, n_steps = 20,
                                env = cbind(ls$env, ls$env), record_steps = c(10, 20),
                                evolution = evolution_control(1, plasticity = 0.2))
  res <- trait_change(two, trait = 2)
  expect_null(res$growth)
  ok <- !is.na(res$community$total)
  expect_equal(res$community$within[ok], (res$community$genetic + res$community$plastic)[ok])
  expect_error(trait_change(two, trait = "nope"), "name one of the traits")
  expect_error(trait_change(two, trait = 2, growth = TRUE), "one evolving trait")
})

test_that("plots are ggplot objects", {
  skip_if_not_installed("ggplot2")
  res <- trait_change(tc_sim())
  expect_s3_class(plot(res), "ggplot")
  expect_s3_class(plot(res, "growth"), "ggplot")
  d <- data.frame(site = "a", species = rep(c("x", "y"), 2), step = rep(1:2, each = 2),
                  abundance = c(10, 10, 5, 15), trait = c(1, 3, 1.5, 3))
  expect_s3_class(plot(trait_change(d)), "ggplot")
  expect_error(plot(trait_change(d), "growth"), "growth = TRUE")
  expect_error(plot(res, colours = "red"), "needs 5 colours")
})

test_that("inputs are checked", {
  d <- data.frame(site = "a", species = rep(c("x", "y"), 2), step = rep(1:2, each = 2),
                  abundance = c(10, 10, 5, 15), trait = c(1, 3, 1.5, 3))
  expect_error(trait_change(d, growth = TRUE), "not a data frame")
  expect_error(trait_change(d, growth = NA), "growth")
  expect_error(trait_change(d, steps = 1), "at least two")
  expect_error(trait_change(d, steps = c(2, 1)), "increasing")
  expect_error(trait_change(d, steps = c(1, 3)), "not recorded")
  expect_error(trait_change(rbind(d, d)), "one row per")
  bad <- d
  bad$abundance[1] <- -1
  expect_error(trait_change(bad), "non-negative")
  bad <- d
  bad$trait[1] <- NA
  expect_error(trait_change(bad), "finite where")
  bad <- d
  bad$trait <- as.character(bad$trait)
  expect_error(trait_change(bad), "numeric")
  bad <- d
  bad$step <- as.character(bad$step)
  expect_error(trait_change(bad), "numbers")
  expect_error(trait_change(list(a = 1)), "simulation with `evolution`")
  sim <- tc_sim(record = NULL)
  expect_error(trait_change(sim), "record_steps")
  sim <- tc_sim()
  sim$traits <- NULL
  expect_null(trait_change(sim)$growth)
  expect_error(trait_change(sim, growth = TRUE), "sigma")
})

test_that("evolving metacommunities, with species named in the settings, are taken", {
  set.seed(1)
  ls <- landscape_grid(5, disp_scale = 0.5)
  anc <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
  sim <- simulate_evolving_metacommunity(ls, anc, speciation_control("geographic", 15),
                                         evolution = evolution_control(10, ne_ratio = 0.5,
                                                                       plasticity = 0.2),
                                         dispersal = 0.01, n_steps = 150,
                                         record_steps = c(100, 125, 150))
  res <- trait_change(sim)
  expect_gt(nrow(res$growth), 0)
  expect_equal(res$growth$total, res$growth$genetic + res$growth$plastic +
                 res$growth$environment, tolerance = 1e-10)
  # a static environment: no plastic change and no direct effect of the environment
  expect_true(all(res$growth$plastic == 0))
  expect_true(all(res$growth$environment == 0))
  sim$traits <- sim$traits[-1, ]
  expect_error(trait_change(sim, growth = TRUE), "lack species")
})

test_that("an interval without shared populations gives an empty growth table", {
  s <- list(abundance = matrix(c(1, 0), 1, dimnames = list("a", c("x", "y"))),
            genetic = matrix(c(0, NA), 1, dimnames = list("a", c("x", "y"))),
            plastic = matrix(c(0, NA), 1, dimnames = list("a", c("x", "y"))), environment = 0)
  s1 <- s
  s1$abundance[] <- c(0, 1)
  s1$genetic[] <- s1$plastic[] <- c(NA, 0)
  niche <- list(sigma = c(x = 1, y = 1), r_max = c(x = 2, y = 2), v = c(x = 0, y = 0))
  expect_identical(nrow(change_growth(s, s1, 1, 2, niche)), 0L)
})
