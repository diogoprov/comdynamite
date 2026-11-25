pop_data <- function() {
  # three species, four sites in a line, two individuals per population and survey
  grid <- expand.grid(individual = 1:2, site = paste0("s", 1:4), species = c("a", "b", "c"),
                      step = c(1, 2), stringsAsFactors = FALSE)
  env <- c(s1 = 0, s2 = 1, s3 = 2, s4 = 3)
  grid$env <- env[grid$site]
  slope <- c(a = 1, b = 2, c = 0.5)[grid$species]
  offset <- c(s1 = 0.1, s2 = -0.2, s3 = 0.3, s4 = -0.1)[grid$site]
  grid$trait <- slope * grid$env + offset + ifelse(grid$individual == 1, -0.5, 0.5) +
    0.2 * (grid$step == 2)
  grid$abundance <- c(s1 = 10, s2 = 20, s3 = 40, s4 = 80)[grid$site] *
    c(a = 1, b = 2, c = 1)[grid$species]
  grid
}
pop_coords <- matrix(c(0, 1, 2, 3, 0, 0, 0, 0), 4, dimnames = list(paste0("s", 1:4), NULL))

test_that("population_stats() returns the same names, with NA where inputs are missing", {
  ind <- pop_data()
  full <- population_stats(ind, abundance = "abundance", step = "step", coords = pop_coords)
  bare <- population_stats(ind[ind$step == 2, ])
  expect_named(full, c("slope", "r2", "sd", "resid_sd", "within", "within_slope", "qst",
                       "near", "far", "dist_cor", "size_cor", "temporal", "temporal_slope",
                       "temporal_step", "temporal_step_slope",
                       "slope_min", "slope_max"))
  expect_named(bare, names(full))
  expect_true(all(is.na(bare[c("within_slope", "near", "far", "dist_cor", "size_cor",
                               "temporal", "temporal_slope")])))
  expect_false(anyNA(full[c("slope", "r2", "within", "qst", "near", "temporal")]))
})

test_that("population_stats() computes regression, within and QST-like statistics by hand", {
  ind <- pop_data()
  last <- ind[ind$step == 2 & ind$species == "a", ]
  means <- tapply(last$trait, last$site, mean)
  vars <- tapply(last$trait, last$site, stats::var)
  fit <- stats::lm(means ~ c(0, 1, 2, 3))
  v_res <- stats::var(stats::residuals(fit))
  by_sp <- population_stats(ind, abundance = "abundance", step = "step", coords = pop_coords,
                            by_species = TRUE)
  a <- by_sp[by_sp$species == "a", ]
  expect_equal(a$slope, unname(stats::coef(fit)[2]))
  expect_equal(a$within, log(mean(vars)))
  expect_equal(a$qst, v_res / (v_res + 2 * mean(vars)))
  expect_identical(by_sp$weight, c(150, 300, 150))
  expect_equal(population_stats(ind, abundance = "abundance")[["slope"]],
               sum(by_sp$slope * by_sp$weight) / 600)
  expect_identical(population_stats(ind, abundance = "abundance", step = "step")[["slope_max"]],
                   max(by_sp$slope))
})

test_that("population_stats() measures the change between the first and the last survey", {
  ind <- pop_data()
  # every mean shifts by 0.2; the sampling variance of the change is 2 * 0.5 / 2 per population
  out <- population_stats(ind, abundance = "abundance", step = "step")
  expect_equal(out[["temporal"]], (0.2^2 - 0.5) / 0.5)
  expect_equal(out[["temporal_slope"]], 0)
  # surveys in between are ignored
  mid <- ind[ind$step == 1, ]
  mid$step <- 1.5
  mid$trait <- mid$trait + 100
  keep <- setdiff(names(population_stats(ind)), c("temporal_step", "temporal_step_slope"))
  expect_identical(population_stats(rbind(ind, mid), step = "step")[keep],
                   population_stats(ind, step = "step")[keep])
})

test_that("population_stats() measures the change between consecutive surveys", {
  ind <- pop_data()
  # two surveys one step apart: the same change as `temporal`, per step
  out <- population_stats(ind, abundance = "abundance", step = "step")
  expect_equal(out[["temporal_step"]], out[["temporal"]])
  expect_equal(out[["temporal_step_slope"]], 0)
  # a third survey two steps later, with every mean shifted by 0.4 more: the changes per step
  # are 0.2^2 - 0.5 and (0.4^2 - 0.5) / 2, relative to the variance within, 0.5
  late <- ind[ind$step == 2, ]
  late$step <- 4
  late$trait <- late$trait + 0.4
  out3 <- population_stats(rbind(ind, late), abundance = "abundance", step = "step")
  expect_equal(out3[["temporal_step"]], mean(c(0.2^2 - 0.5, (0.4^2 - 0.5) / 2)) / 0.5)
  expect_equal(out3[["temporal"]], (0.6^2 - 0.5) / 0.5)
  # a step that is not numeric gives NA
  ind$step <- as.character(ind$step)
  expect_true(is.na(population_stats(ind, step = "step")[["temporal_step"]]))
})

test_that("population_stats() leaves out species with too few populations", {
  ind <- pop_data()
  ind <- ind[!(ind$species == "c" & ind$site %in% c("s3", "s4")), ]
  by_sp <- population_stats(ind, by_species = TRUE)
  expect_identical(by_sp$weight[by_sp$species == "c"], 0)
  expect_true(all(is.na(unlist(by_sp[by_sp$species == "c", -(1:2)]))))
  # single individuals give no variance, so the population is not used
  one <- ind[ind$step == 2 & !(ind$species == "a" & ind$individual == 2 & ind$site == "s1"), ]
  expect_identical(population_stats(one[one$species == "a" & one$site != "s2", ],
                                    by_species = TRUE)$weight, 0)
})

test_that("population_stats() gives distance statistics that reflect spatial structure", {
  ind <- pop_data()
  out <- population_stats(ind, coords = pop_coords, near = 1, far = 2)
  expect_true(is.finite(out[["near"]]) && is.finite(out[["far"]]))
  expect_true(all(is.na(population_stats(ind, coords = pop_coords[1:3, ])[
    c("near", "far", "dist_cor")])))
})

test_that("population_stats() works on sample_individuals() output", {
  set.seed(1)
  land <- landscape_grid(side = 5, env_bandwidth = 2)
  traits <- data.frame(optimum = c(30, 60), sigma = 15, r_max = 5, row.names = c("a", "b"))
  alpha <- matrix(0.01, 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
  diag(alpha) <- 0.02
  sim <- simulate_metacommunity(land, traits, alpha, dispersal = 0.05, n_steps = 60,
                                record_steps = c(50, 60),
                                evolution = evolution_control(4, ne_ratio = 0.5, env_var = 1))
  ind <- sample_individuals(sim, n = 8, source = "evolution", step = c(50, 60))
  coords <- land$coords
  rownames(coords) <- rownames(sim$abundance)
  out <- population_stats(ind, abundance = "abundance", step = "step", coords = coords)
  expect_type(out, "double")
  expect_true(is.finite(out[["within"]]))
  expect_lt(abs(exp(out[["within"]]) - 5), 2)
})

test_that("population_stats() checks its arguments", {
  ind <- pop_data()
  expect_error(population_stats(as.list(ind)), "data frame")
  expect_error(population_stats(ind, trait = "height"), "no column height")
  ind_chr <- ind
  ind_chr$env <- as.character(ind_chr$env)
  expect_error(population_stats(ind_chr), "`env` must name a numeric column")
  expect_error(population_stats(ind, near = 3, far = 1), "far")
  expect_error(population_stats(ind, coords = unname(pop_coords)), "row names")
  expect_error(population_stats(ind, by_species = NA), "by_species")
})
