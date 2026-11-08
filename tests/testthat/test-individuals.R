small_sim <- function(seed = 1, record_steps = NULL) {
  set.seed(seed)
  simulate_metacommunity(landscape_implicit(12), species_traits(8),
                         competition_matrix(8, "none"), n_steps = 60, seed_steps = 1:5,
                         record_steps = record_steps)
}

test_that("sample_individuals() draws n individuals per population, or all of a small one", {
  sim <- small_sim()
  ind <- sample_individuals(sim, n = 5)
  counts <- table(paste(ind$site, ind$species))
  expected <- pmin(sim$abundance[sim$abundance > 0], 5)
  expect_equal(sort(as.vector(counts)), sort(as.vector(expected)))
  expect_named(ind, c("site", "species", "individual", "abundance", "env", "optimum",
                      "population_mean", "trait"))
  expect_true(all(ind$abundance >= ind$individual))
  expect_false(is.unsorted(match(ind$site, rownames(sim$abundance))))
})

test_that("sample_individuals() thins individuals binomially with fraction", {
  sim <- small_sim()
  ind <- sample_individuals(sim, fraction = 0.5)
  counts <- table(paste(ind$site, ind$species))
  expect_true(abs(nrow(ind) / sum(sim$abundance) - 0.5) < 0.1)
  expect_true(all(counts <= ind$abundance[match(names(counts),
                                                paste(ind$site, ind$species))]))
})

test_that("population means follow the optimum, the environment and local adaptation", {
  sim <- small_sim()
  none <- sample_individuals(sim, intraspecific_sd = 0)
  expect_equal(none$trait, none$optimum)
  full <- sample_individuals(sim, local_adaptation = 1, intraspecific_sd = 0)
  expect_equal(full$trait, full$env)
  half <- sample_individuals(sim, local_adaptation = 0.5)
  expect_equal(half$population_mean, (half$optimum + half$env) / 2)
  sp <- rownames(sim$traits)
  per_sp <- sample_individuals(sim, local_adaptation = stats::setNames(rep(c(0, 1), 4), sp),
                               intraspecific_sd = 0)
  adapted <- per_sp$species %in% sp[c(2, 4, 6, 8)]
  expect_equal(per_sp$trait[adapted], per_sp$env[adapted])
  expect_equal(per_sp$trait[!adapted], per_sp$optimum[!adapted])
})

test_that("sample_individuals() uses recorded steps and a given environment", {
  sim <- small_sim(record_steps = c(30, 60))
  at30 <- sample_individuals(sim, step = 30, n = 1e6, intraspecific_sd = 0)
  expect_equal(sum(at30$abundance[at30$individual == 1]), sum(sim$history[, , "step30"]))
  expect_error(sample_individuals(sim, step = 31), "not recorded")
  expect_error(sample_individuals(small_sim(), step = 30), "record_steps")
  env <- seq(0, 1, length.out = 12)
  ind <- sample_individuals(sim, env = env, local_adaptation = 1, intraspecific_sd = 0)
  expect_equal(ind$trait, env[match(ind$site, rownames(sim$abundance))])
  expect_error(sample_individuals(sim, env = 1:3), "one number per site")
})

test_that("sample_individuals() checks its arguments", {
  sim <- small_sim()
  expect_error(sample_individuals(sim, local_adaptation = 2), "between 0 and 1")
  expect_error(sample_individuals(sim, intraspecific_sd = -1), "negative")
  expect_error(sample_individuals(sim, population_sd = 1:3), "one per species")
  expect_error(sample_individuals(sim, local_adaptation = c(sp1 = 0.5)), "no value for")
  expect_error(sample_individuals(list()), "mc_sim")
  expect_error(sample_individuals(sim, local_adaptation = "a"), "must be numbers")
  expect_error(sample_individuals(sim, intraspecific_sd = NA_real_), "must be numbers")
})

test_that("unnamed per-species values follow the order of the species", {
  sim <- small_sim()
  sp <- rownames(sim$traits)
  ind <- sample_individuals(sim, local_adaptation = rep(c(0, 1), 4), intraspecific_sd = 0)
  adapted <- ind$species %in% sp[c(2, 4, 6, 8)]
  expect_equal(ind$trait[adapted], ind$env[adapted])
  expect_equal(ind$trait[!adapted], ind$optimum[!adapted])
})

# a hand-made data set with known variances
toy <- data.frame(site = rep(c("a", "b"), each = 6),
                  species = rep(rep(c("x", "y"), each = 3), 2),
                  trait = c(1, 2, 3, 7, 8, 9, 2, 3, 4, 2, 3, 4))

test_that("tstats() gives the variance ratios of Violle et al. (2012)", {
  res <- tstats(toy, nperm = 0)
  expect_named(res, c("site", "trait", "n_populations", "T_IP.IC", "T_IC.IR", "T_PC.PR"))
  expect_equal(res$site, c("a", "b"))
  expect_equal(res$n_populations, c(2, 2))
  z <- toy$trait
  expect_equal(res$T_IP.IC, c(1 / var(z[1:6]), 1 / var(z[7:12])))
  expect_equal(res$T_IC.IR, c(var(z[1:6]), var(z[7:12])) / var(z))
  means <- c(2, 8, 3, 3)
  expect_equal(res$T_PC.PR, c(var(means[1:2]), var(means[3:4])) / var(means))
})

test_that("tstats() gives NA where a ratio is undefined and handles several traits", {
  one_species <- rbind(toy, data.frame(site = "c", species = "x", trait = c(5, 6)))
  one_species$other <- one_species$trait * 2
  res <- tstats(one_species, traits = c("trait", "other"), nperm = 0)
  expect_equal(nrow(res), 6)
  expect_true(is.na(res$T_PC.PR[res$site == "c"][1]))
  expect_equal(res$T_IP.IC[res$trait == "trait"], res$T_IP.IC[res$trait == "other"])
  expect_error(tstats(toy, traits = "height"), "no column height")
  expect_error(tstats(transform(toy, trait = as.character(trait))), "numeric")
})

test_that("tstats() matches cati::Tstats(), also with site names that are prefixes", {
  skip_if_not_installed("cati")
  sim <- small_sim()
  set.seed(2)
  ind <- sample_individuals(sim, local_adaptation = 0.5, population_sd = 2,
                            intraspecific_sd = 3)
  base <- tstats(ind, nperm = 0)
  res <- tstats(ind, method = "cati", nperm = NULL)
  expect_s3_class(res, "Tstats")
  expect_true(any(base$site == "site1") && any(base$site == "site10"))
  for (stat in c("T_IP.IC", "T_IC.IR", "T_PC.PR")) {
    expect_equal(unname(res$Tstats[[stat]][base$site, 1]), base[[stat]])
  }
  with_null <- tstats(ind, method = "cati", nperm = 9)
  expect_equal(dim(with_null$Tstats$T_IP.IC_nm)[3], length(unique(ind$site)))
})

test_that("tstats() explains how to get cati", {
  local_mocked_bindings(cati_available = function() FALSE)
  expect_error(tstats(toy, method = "cati"), "install.packages")
})

test_that("local adaptation raises T_IP.IC and lowers T_PC.PR", {
  sim <- small_sim()
  set.seed(3)
  ratios <- function(h) {
    res <- tstats(sample_individuals(sim, local_adaptation = h, intraspecific_sd = 2),
                  nperm = 0)
    colMeans(res[c("T_IP.IC", "T_IC.IR", "T_PC.PR")], na.rm = TRUE)
  }
  none <- ratios(0)
  strong <- ratios(0.9)
  expect_gt(strong[["T_IP.IC"]], none[["T_IP.IC"]])
  expect_lt(strong[["T_PC.PR"]], none[["T_PC.PR"]])
  expect_lt(strong[["T_IC.IR"]], none[["T_IC.IR"]])
})

test_that("tstats() leaves out missing values trait by trait", {
  with_na <- toy
  with_na$other <- with_na$trait
  with_na$other[1] <- NA
  res <- tstats(with_na, traits = c("trait", "other"), nperm = 0)
  expect_equal(res$T_IP.IC[res$trait == "other"],
               tstats(toy[-1, ], nperm = 0)$T_IP.IC)
})

test_that("tstats() adds effect sizes and p-values from the null models", {
  set.seed(1)
  res <- tstats(toy, nperm = 49)
  expect_named(res, c("site", "trait", "n_populations", "T_IP.IC", "T_IC.IR", "T_PC.PR",
                      "ses_T_IP.IC", "p_T_IP.IC", "ses_T_IC.IR", "p_T_IC.IR",
                      "ses_T_PC.PR", "p_T_PC.PR"))
  null <- attr(res, "null")
  expect_equal(dim(null), c(49, 2, 3))
  # permuting within a site keeps its variance, and the two species of site a never mix less
  # than observed, so its T_IP.IC is the lowest possible
  expect_true(all(null[, 1, "T_IP.IC"] >= res$T_IP.IC[1] - 1e-12))
  expect_equal(res$p_T_IP.IC[1], (1 + sum(null[, 1, "T_IP.IC"] <= res$T_IP.IC[1])) / 50)
  expect_equal(res$ses_T_IP.IC[1],
               (res$T_IP.IC[1] - mean(null[, 1, "T_IP.IC"])) / sd(null[, 1, "T_IP.IC"]))
  greater <- tstats(toy, nperm = 49, alternative = "greater")
  expect_true(all(greater$p_T_IP.IC > 0 & greater$p_T_IP.IC <= 1))
  two <- tstats(toy, nperm = 49, alternative = "two.sided")
  expect_true(all(two$p_T_IP.IC <= 1))
  expect_identical(names(tstats(toy, nperm = NULL)), names(tstats(toy, nperm = 0)))
  expect_error(tstats(toy, nperm = -1))
})

test_that("null models of tstats() keep what they should", {
  set.seed(2)
  d <- data.frame(site = rep(c("a", "b", "c"), each = 20),
                  species = rep(rep(c("x", "y"), each = 10), 3))
  d$trait <- rnorm(60, mean = rep(c(0, 10, 20), each = 20) + (d$species == "y") * 3)
  g <- tstats_groups(d$site, d$species, unique(d$site))
  null <- tstats_null(d$trait, g, 30)
  # individuals drawn from the whole region: sites lose their own range, T_IC.IR near 1
  expect_gt(mean(null[, , "T_IC.IR"]), 0.7)
  # population means drawn from the whole region: T_PC.PR near 1 on average
  expect_gt(mean(null[, , "T_PC.PR"], na.rm = TRUE), 0.5)
  # random data: no structure, effect sizes small on average
  d$trait <- rnorm(60)
  res <- tstats(d, nperm = 99)
  expect_lt(abs(mean(res$ses_T_IP.IC)), 1.5)
  expect_lt(abs(mean(res$ses_T_IC.IR)), 1.5)
})

test_that("character displacement pushes co-occurring populations apart", {
  pushed <- displace(c(0, 1), c(10, 10), strength = 2, width = 2)
  expect_lt(pushed[1], 0)
  expect_gt(pushed[2], 1)
  expect_equal(pushed[1] + pushed[2], 1)
  # abundant competitors push harder
  uneven <- displace(c(0, 1), c(1, 99), strength = 2, width = 2)
  expect_lt(uneven[1], pushed[1])
  # the equilibrium satisfies its equation
  mu <- displace(c(0, 1, 5), c(5, 3, 2), strength = 3, width = 2)
  gap <- outer(mu, mu, "-")
  rhs <- c(0, 1, 5) + 3 * as.vector((gap * exp(-gap^2 / 8)) %*% (c(5, 3, 2) / 10))
  expect_equal(mu, rhs, tolerance = 1e-6)
  expect_equal(displace(c(2, 2), c(1, 1), 5, 1), c(2, 2))
  expect_warning(displace(c(0, 1), c(1, 1), 2, 2, max_iter = 1), "did not converge")
})

test_that("sample_individuals() applies character displacement within sites", {
  sim <- small_sim()
  set.seed(4)
  base <- sample_individuals(sim, intraspecific_sd = 0)
  set.seed(4)
  apart <- sample_individuals(sim, intraspecific_sd = 0, displacement = 3)
  spread <- function(ind) {
    pops <- unique(ind[c("site", "population_mean")])
    mean(tapply(pops$population_mean, pops$site, var), na.rm = TRUE)
  }
  expect_gt(spread(apart), spread(base))
  set.seed(4)
  narrow <- sample_individuals(sim, intraspecific_sd = 0, displacement = 3,
                               displacement_width = 0.01)
  expect_equal(narrow$trait, base$trait, tolerance = 1e-6)
  expect_error(sample_individuals(sim, displacement = -1))
  expect_error(sample_individuals(sim, displacement = 1, displacement_width = 0))
  set.seed(5)
  t_none <- tstats(sample_individuals(sim, intraspecific_sd = 2), nperm = 0)
  set.seed(5)
  t_apart <- tstats(sample_individuals(sim, intraspecific_sd = 2, displacement = 5), nperm = 0)
  expect_lt(mean(t_apart$T_IP.IC, na.rm = TRUE), mean(t_none$T_IP.IC, na.rm = TRUE))
})

test_that("tstats() stacks the null values of several traits", {
  two <- transform(toy, other = trait^2)
  set.seed(6)
  res <- tstats(two, traits = c("trait", "other"), nperm = 19)
  expect_equal(nrow(res), 4)
  expect_equal(dim(attr(res, "null")), c(19, 4, 3))
  expect_false(anyNA(res$p_T_IP.IC))
})

test_that("character displacement leaves sites with a single species alone", {
  sim <- small_sim()
  sim$abundance[1, ] <- 0
  sim$abundance[1, 1] <- 20
  set.seed(7)
  base <- sample_individuals(sim, intraspecific_sd = 0)
  set.seed(7)
  apart <- sample_individuals(sim, intraspecific_sd = 0, displacement = 3)
  first <- rownames(sim$abundance)[1]
  expect_equal(apart$trait[apart$site == first], base$trait[base$site == first])
})

test_that("the species null of T_IP.IC moves populations only within species", {
  sim <- small_sim()
  set.seed(2)
  ind <- sample_individuals(sim, n = 6, intraspecific_sd = 1)
  res <- tstats(ind, nperm = 20, null_ip = "species", alternative = "two.sided")
  base <- tstats(ind, nperm = 0)
  expect_equal(res[names(base)], base)
  # when every species occurs at one site, nothing can move and the null equals the observation
  single <- ind
  single$species <- paste(single$site, single$species)
  kept <- tstats(single, nperm = 5, null_ip = "species")
  expect_equal(attr(kept, "null")[, , "T_IP.IC"],
               matrix(kept$T_IP.IC, 5, nrow(kept), byrow = TRUE))
  expect_true(all(is.finite(res$ses_T_IP.IC) | is.na(res$ses_T_IP.IC)))
  expect_true(all(res$p_T_IP.IC <= 1, na.rm = TRUE))
  expect_error(tstats(ind, null_ip = "region"))
})

test_that("the species null of T_IP.IC tells convergence from divergence", {
  set.seed(3)
  sim <- simulate_metacommunity(landscape_implicit(15), species_traits(10),
                                competition_matrix(10, "none"), n_steps = 100,
                                seed_steps = 1:5)
  ses <- function(...) {
    ind <- sample_individuals(sim, n = 10, intraspecific_sd = 1, ...)
    mean(tstats(ind, nperm = 39, null_ip = "species")$ses_T_IP.IC, na.rm = TRUE)
  }
  expect_gt(ses(local_adaptation = 0.9), 2)
  expect_lt(ses(displacement = 10), 0)
})

test_that("group_sum() gives zero to groups without members", {
  expect_equal(group_sum(c(1, 2, 3), c(1, 3, 3), 4), c(1, 0, 5, 0))
  expect_equal(group_mean(c(1, 2, 3), c(1, 3, 3), 3), c(1, NaN, 2.5))
})

test_that("population_shifts() recovers local adaptation and displacement", {
  set.seed(4)
  sim <- simulate_metacommunity(landscape_implicit(15), species_traits(10),
                                competition_matrix(10, "none"), n_steps = 100,
                                seed_steps = 1:5)
  w <- mean(sim$traits$sigma)
  # without displacement the width is not identified, so it is fixed here
  local <- population_shifts(sample_individuals(sim, n = 200, local_adaptation = 0.6,
                                                intraspecific_sd = 0.5), width = w)
  expect_equal(local$local_adaptation, 0.6, tolerance = 0.05)
  expect_lt(abs(local$displacement), 0.5)
  both <- population_shifts(sample_individuals(sim, n = 200, local_adaptation = 0.4,
                                               displacement = 3, intraspecific_sd = 0.5),
                            width = w)
  expect_equal(both$local_adaptation, 0.4, tolerance = 0.1)
  expect_equal(both$displacement, 3, tolerance = 0.2)
  expect_gt(both$r_squared, 0.8)
  expect_equal(attr(both, "profile")$width, w)
  expect_named(both, c("local_adaptation", "displacement", "width", "r_squared",
                       "n_populations", "n_species"))
  free <- population_shifts(sample_individuals(sim, n = 200, local_adaptation = 0.4,
                                               displacement = 3, intraspecific_sd = 0.5))
  expect_equal(nrow(attr(free, "profile")), 22)
  expect_equal(free$width, attr(free, "profile")$width[which.min(attr(free, "profile")$rss)])
})

test_that("population_shifts() checks its input", {
  sim <- small_sim()
  ind <- sample_individuals(sim, n = 4)
  expect_error(population_shifts(ind[setdiff(names(ind), "env")]), "no column env")
  ind$env <- as.character(ind$env)
  expect_error(population_shifts(ind), "numeric")
  one <- sample_individuals(sim, n = 1)
  expect_error(population_shifts(one), "two populations")
  # with no competitors anywhere, displacement cannot be estimated
  solo <- sample_individuals(sim, n = 5)
  solo$site <- paste(solo$site, solo$species)
  solo$env <- solo$env + rnorm(nrow(solo))[match(solo$site, unique(solo$site))]
  solo$species <- "a"
  expect_true(is.na(population_shifts(solo)$displacement))
})

test_that("tstats() returns NA for a trait without values", {
  sim <- small_sim()
  ind <- sample_individuals(sim, n = 4)
  ind$empty <- NA_real_
  res <- tstats(ind, traits = c("trait", "empty"), nperm = 0)
  expect_true(all(is.na(res$T_IP.IC[res$trait == "empty"])))
  expect_false(all(is.na(res$T_IP.IC[res$trait == "trait"])))
})
