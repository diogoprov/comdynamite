# data with known local adaptation (one slope per variable) and displacement
make_shifts <- function(n_sites = 25, n_species = 8, n_ind = 6, h = 0.5, h2 = 0, c = 0,
                        w = 0.5, within_sd = 0.2, species_h = NULL, seed = 1) {
  set.seed(seed)
  theta <- seq(0, 3, length.out = n_species)
  env <- stats::rnorm(n_sites)
  env2 <- stats::rnorm(n_sites)
  present <- matrix(stats::runif(n_sites * n_species) < 0.6, n_sites, n_species)
  abundance <- matrix(stats::rpois(n_sites * n_species, 10) + 1, n_sites, n_species)
  slope <- if (is.null(species_h)) rep(h, n_species) else species_h
  rows <- list()
  for (i in seq_len(n_sites)) {
    j <- which(present[i, ])
    if (length(j) == 0) next
    target <- theta[j] + slope[j] * env[i] + h2 * env2[i]
    m <- if (c > 0 && length(j) > 1) displace(target, abundance[i, j], c, w) else target
    rows[[i]] <- data.frame(site = paste0("s", i), species = paste0("sp", rep(j, each = n_ind)),
                            env = env[i], env2 = env2[i],
                            abundance = rep(abundance[i, j], each = n_ind),
                            trait = rep(m, each = n_ind) + stats::rnorm(length(j) * n_ind, 0,
                                                                         within_sd))
  }
  out <- do.call(rbind, rows)
  attr(out, "community") <- structure(abundance * present,
                                      dimnames = list(paste0("s", seq_len(n_sites)),
                                                      paste0("sp", seq_len(n_species))))
  out
}

test_that("population_shifts() keeps its results and class with one variable", {
  ind <- make_shifts(c = 2)
  fit <- population_shifts(ind, width = 0.5)
  expect_s3_class(fit, "population_shifts")
  expect_named(fit, c("local_adaptation", "displacement", "width", "r_squared",
                      "n_populations", "n_species"))
  expect_equal(fit$local_adaptation, 0.5, tolerance = 0.1)
  expect_equal(fit$displacement, 2, tolerance = 0.25)
  expect_named(attr(fit, "profile"), c("width", "rss", "r_squared", "displacement",
                                       "local_adaptation"))
})

test_that("the grid of widths starts at the spread within populations", {
  ind <- make_shifts(within_sd = 0.6, c = 2)
  fit <- suppressWarnings(population_shifts(ind))
  prof <- attr(fit, "profile")
  pops <- split(ind$trait, paste(ind$site, ind$species))
  pooled <- sqrt(sum(vapply(pops, function(v) sum((v - mean(v))^2), 0)) /
                   (nrow(ind) - length(pops)))
  expect_equal(min(prof$width), pooled)
  expect_equal(attr(fit, "floor"), pooled)
  expect_equal(nrow(prof), 22)
  narrow <- suppressWarnings(population_shifts(make_shifts(within_sd = 0.01, c = 2)))
  expect_null(attr(narrow, "floor"))
})

test_that("an unidentified width gives a warning", {
  ind <- make_shifts(c = 0, seed = 3)
  expect_warning(population_shifts(ind), "poorly identified")
  expect_no_warning(population_shifts(ind, width = 0.5))
  strong <- make_shifts(c = 4, w = 0.6, n_ind = 20, seed = 2)
  expect_no_warning(population_shifts(strong))
})

test_that("several environmental variables get their own slopes", {
  ind <- make_shifts(h = 0.5, h2 = -0.3, c = 2)
  fit <- population_shifts(ind, env = c("env", "env2"), width = 0.5)
  expect_named(fit, c("local_adaptation_env", "local_adaptation_env2", "displacement",
                      "width", "r_squared", "n_populations", "n_species"))
  expect_equal(fit$local_adaptation_env, 0.5, tolerance = 0.1)
  expect_equal(fit$local_adaptation_env2, -0.3, tolerance = 0.15)
  expect_true(all(c("local_adaptation_env", "local_adaptation_env2") %in%
                    names(attr(fit, "profile"))))
})

test_that("species_slopes gives each species its own slope", {
  slopes <- seq(0, 0.8, length.out = 8)
  ind <- make_shifts(species_h = slopes, n_sites = 40, n_ind = 10)
  fit <- population_shifts(ind, species_slopes = TRUE, width = 0.5)
  per <- attr(fit, "species_slopes")
  expect_named(per, c("species", "variable", "slope"))
  est <- per$slope[match(paste0("sp", 1:8), per$species)]
  expect_gt(stats::cor(est, slopes), 0.95)
  expect_equal(fit$local_adaptation, mean(per$slope))
})

test_that("community adds competitors that were not measured", {
  ind <- make_shifts(c = 2)
  community <- attr(ind, "community")
  # with the measured populations only (and abundances at least the number measured),
  # community and abundance give the same fit
  ind$abundance <- pmax(ind$abundance, 6)
  only_measured <- community * 0
  pops <- unique(ind[c("site", "species", "abundance")])
  only_measured[cbind(pops$site, pops$species)] <- pops$abundance
  set.seed(1)
  a <- population_shifts(ind, width = 0.5)
  set.seed(1)
  b <- population_shifts(ind, width = 0.5, community = only_measured, abundance = NULL)
  expect_equal(a$displacement, b$displacement)
  # drop the measurements of three species at half the sites: community still counts them,
  # which reduces the attenuation of c
  est <- vapply(1:8, function(s) {
    full <- make_shifts(c = 2, n_ind = 10, within_sd = 0.1, seed = 100 + s)
    gone <- full$species %in% c("sp3", "sp5", "sp7") &
      as.integer(sub("s", "", full$site)) %% 2 == 0
    set.seed(s)
    without <- population_shifts(full[!gone, ], width = 0.5)$displacement
    set.seed(s)
    with <- population_shifts(full[!gone, ], width = 0.5, community = attr(full, "community"),
                              abundance = NULL)$displacement
    c(without, with)
  }, numeric(2))
  expect_gt(mean(est[2, ]) - mean(est[1, ]), 0.15)
  # measured populations absent from the census take the number measured
  zeroed <- community
  zeroed[pops$site[1], pops$species[1]] <- 0
  expect_no_error(population_shifts(ind, width = 0.5, community = zeroed))
  # never measured species are dropped with a warning
  extra <- cbind(community, ghost = 5)
  expect_warning(population_shifts(ind, width = 0.5, community = extra), "ghost")
  expect_error(population_shifts(ind, community = community[-1, ]), "no row for site")
  expect_error(population_shifts(ind, community = community[, -1]), "no column for species")
  bad <- community
  bad[1, 1] <- -1
  expect_error(population_shifts(ind, community = bad), "non-negative")
  expect_error(population_shifts(ind, community = unname(community)), "row names")
  expect_error(population_shifts(ind, abundance = NULL), "abundance")
})

test_that("the bootstrap tests and corrects displacement", {
  strong <- make_shifts(c = 3, w = 0.5, n_ind = 4, seed = 5)
  set.seed(2)
  fit <- population_shifts(strong, width = 0.5, nboot = 39)
  expect_true(all(c("displacement_corrected", "displacement_low", "displacement_high",
                    "p_displacement", "local_adaptation_low", "local_adaptation_high") %in%
                    names(fit)))
  expect_equal(fit$p_displacement, 1 / 40)
  expect_lt(fit$displacement, 3)
  expect_lt(abs(fit$displacement_corrected - 3), abs(fit$displacement - 3))
  expect_lt(fit$displacement_low, fit$displacement_high)
  expect_lt(fit$local_adaptation_low, 0.5)
  expect_gt(fit$local_adaptation_high, 0.5)
  boot <- attr(fit, "bootstrap")
  expect_equal(as.vector(table(boot$set)[c("null", "fitted")]), c(39, 39))
  null <- make_shifts(c = 0, seed = 6)
  set.seed(3)
  fit0 <- population_shifts(null, width = 0.5, nboot = 39)
  expect_gt(fit0$p_displacement, 0.05)
  # the bootstrap also runs with an estimated width, several variables and species slopes
  set.seed(4)
  many <- suppressWarnings(population_shifts(make_shifts(c = 2, h2 = 0.2, seed = 7),
                                             env = c("env", "env2"), species_slopes = TRUE,
                                             nboot = 5))
  expect_true(all(c("local_adaptation_env_low", "local_adaptation_env2_high") %in%
                    names(many)))
  community <- attr(make_shifts(c = 2, seed = 8), "community")
  set.seed(5)
  comm <- population_shifts(make_shifts(c = 2, seed = 8), width = 0.5, community = community,
                            abundance = NULL, nboot = 5)
  expect_false(is.na(comm$p_displacement))
})

test_that("the bootstrap gives NA for displacement without competitors", {
  ind <- make_shifts(seed = 9)
  ind$site <- paste(ind$site, ind$species)
  ind$env <- ind$env + stats::rnorm(nrow(ind))[match(ind$site, unique(ind$site))]
  set.seed(1)
  fit <- population_shifts(ind, nboot = 5)
  expect_true(is.na(fit$displacement))
  expect_true(is.na(fit$p_displacement))
  expect_false(is.na(fit$local_adaptation_low))
})

test_that("plot() draws the profile of a fit", {
  skip_if_not_installed("ggplot2")
  ind <- make_shifts(c = 3, w = 0.5, within_sd = 0.4, seed = 10)
  set.seed(1)
  fit <- suppressWarnings(population_shifts(ind, nboot = 9))
  p <- plot(fit)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_equal(nlevels(b$data[[2]]$PANEL), 3)
  expect_true(any(vapply(b$data, function(l) "xintercept" %in% names(l) &&
                           any(abs(l$xintercept - log10(fit$width)) < 1e-8), logical(1))))
  expect_error(plot(population_shifts(ind, width = 0.5)), "fixed")
})

test_that("unmeasured competitors follow the gradient with the fitted slopes", {
  ind <- make_shifts(h = 1, c = 2, n_ind = 8, seed = 11)
  community <- attr(ind, "community")
  gone <- ind$species == "sp4" & as.integer(sub("s", "", ind$site)) %% 2 == 0
  kept <- ind[!gone, ]
  set.seed(1)
  fit <- population_shifts(kept, width = 0.5, community = community, abundance = NULL)
  set.seed(1)
  d <- data.frame(trait = kept$trait, site = kept$site, species = kept$species,
                  .env1 = kept$env)
  prep <- shifts_prepare(d, ".env1", community)
  first <- shifts_fit(prep, 0.5, FALSE, "local_adaptation")$table$local_adaptation
  full <- shifts_fit_full(prep, 0.5, FALSE, "local_adaptation")
  expect_equal(full$table$displacement, fit$displacement)
  pops <- unique(d[c("site", "species", ".env1")])
  sp4 <- pops[pops$species == "sp4", ]
  members <- full$prep$members
  site_names <- unique(pops$site)
  for (k in seq_along(members)) {
    m <- members[[k]]
    if ("sp4" %in% m$ext_species) {
      e <- pops$.env1[pops$site == site_names[k]][1]
      j <- match("sp4", m$ext_species)
      expect_equal(m$ext_mean[j] - m$ext_base[j],
                   first * (e - mean(sp4$.env1)))
    }
  }
  expect_true(any(vapply(members, function(m) "sp4" %in% m$ext_species, logical(1))))
})

test_that("a community with one species keeps its names", {
  ind <- make_shifts(n_species = 1, n_sites = 20, seed = 12)
  community <- attr(ind, "community")
  expect_equal(ncol(community), 1)
  fit <- population_shifts(ind, width = 0.5, community = community, abundance = NULL)
  expect_true(is.na(fit$displacement))
})

test_that("the bootstrap centres the p-value and corrects from the generating value", {
  null <- make_shifts(c = 0, seed = 13)
  set.seed(4)
  fit <- population_shifts(null, width = 0.5, nboot = 19)
  boot <- attr(fit, "bootstrap")
  c_null <- boot$displacement[boot$set == "null"]
  c_fitted <- boot$displacement[boot$set == "fitted"]
  centre <- mean(c_null)
  expect_equal(fit$p_displacement,
               (1 + sum(abs(c_null - centre) >= abs(fit$displacement - centre))) / 20)
  c0 <- max(fit$displacement, 0)
  expect_equal(fit$displacement_corrected, fit$displacement - (mean(c_fitted) - c0))
  q <- stats::quantile(c_fitted, c(0.975, 0.025), names = FALSE)
  expect_equal(c(fit$displacement_low, fit$displacement_high), fit$displacement - (q - c0))
})
