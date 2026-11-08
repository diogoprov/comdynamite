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
  res <- tstats(toy)
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
  res <- tstats(one_species, traits = c("trait", "other"))
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
  base <- tstats(ind)
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
    res <- tstats(sample_individuals(sim, local_adaptation = h, intraspecific_sd = 2))
    colMeans(res[c("T_IP.IC", "T_IC.IR", "T_PC.PR")], na.rm = TRUE)
  }
  none <- ratios(0)
  strong <- ratios(0.9)
  expect_gt(strong[["T_IP.IC"]], none[["T_IP.IC"]])
  expect_lt(strong[["T_PC.PR"]], none[["T_PC.PR"]])
  expect_lt(strong[["T_IC.IR"]], none[["T_IC.IR"]])
})
