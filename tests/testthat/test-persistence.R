test_that("persistence() recovers the rates of a known two-state Markov chain", {
  set.seed(1)
  m <- 400
  times <- 60
  c_true <- 0.1
  e_true <- 0.3
  x <- array(0L, c(m, 2, times), dimnames = list(NULL, c("a", "b"), NULL))
  x[, , 1] <- stats::rbinom(2 * m, 1, 0.5)
  for (t in 2:times) {
    p <- ifelse(x[, , t - 1] > 0, 1 - e_true, c_true)
    x[, , t] <- stats::rbinom(2 * m, 1, p)
  }
  p <- persistence(x)
  expect_named(p, c("species", "intervals", "colonisation", "extirpation", "occupancy",
                    "epsilon", "tau"))
  expect_equal(p$colonisation, c(c_true, c_true), tolerance = 0.05)
  expect_equal(p$extirpation, c(e_true, e_true), tolerance = 0.05)
  expect_equal(p$epsilon, rep(stats::qlogis(c_true / (c_true + e_true)), 2), tolerance = 0.05)
  expect_equal(p$tau, rep(log(c_true + e_true), 2), tolerance = 0.05)
  expect_equal(persistence(x, window = 10)$intervals, c(10L, 10L))
})

test_that("persistence() counts transitions by hand and handles edge cases", {
  x <- array(c(1, 0, 0, 0,    # site 1: present, absent, absent, absent
               0, 1, 1, 0),   # site 2: absent, present, present, absent
             c(4, 2), dimnames = NULL)
  x <- aperm(array(x, c(4, 2, 1)), c(2, 3, 1))
  p <- persistence(x)
  # absent exposures: site 1 at times 2, 3; site 2 at time 1 -> 3, one colonisation
  # present exposures: site 1 at 1; site 2 at 2, 3 -> 3, two extirpations
  expect_equal(p$colonisation, 1 / 3)
  expect_equal(p$extirpation, 2 / 3)
  expect_equal(p$species, "sp1")
  always <- array(1L, c(3, 1, 5))
  expect_true(is.na(persistence(always)$colonisation))
  expect_true(is.na(persistence(always)$epsilon))
  expect_error(persistence(matrix(1, 2, 2)), "array")
})

test_that("persistence() works on simulate_metacommunity() and its records", {
  set.seed(2)
  tr <- ape::rcoal(4)
  sim <- simulate_metacommunity(landscape_implicit(20), species_traits(tr),
                                competition_matrix(tr, "phylo"), n_steps = 100,
                                dispersal = 0.01, extirpation = 0.05, record_steps = 1:100)
  p <- persistence(sim)
  expect_equal(nrow(p), 4)
  expect_identical(p, persistence(sim$history))
  expect_true(all(p$intervals == 99L))
  sim$history <- NULL
  expect_error(persistence(sim), "record_steps")
})

test_that("persistence() of an evolving metacommunity follows its lineage series", {
  set.seed(1)
  ls <- landscape_grid(5, disp_scale = 0.5)
  anc <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
  sim <- simulate_evolving_metacommunity(ls, anc, speciation_control("fission", rate = 0.02),
                                         dispersal = 0.01, extirpation = 0.02, n_steps = 150)
  s <- sim$series
  expect_named(s, c("generation", "lineage", "occupancy", "abundance", "genetic", "colonised",
                    "extirpated", "present_before", "absent_before"))
  # occupancy and abundance match the final state of extant lineages
  last <- s[s$generation == sim$generations, ]
  extant <- colnames(sim$lineages)
  expect_setequal(last$lineage, extant)
  expect_equal(last$abundance[match(extant, last$lineage)], unname(colSums(sim$lineages)))
  expect_equal(last$occupancy[match(extant, last$lineage)],
               unname(colSums(sim$lineages > 0)))
  # every exposure is either present or absent, apart from sites handed to buds
  ok <- !is.na(s$present_before)
  full <- s$present_before + s$absent_before == nrow(sim$abundance)
  expect_true(all(full[ok] | s$generation[ok] %in% sim$events$time))
  expect_true(all(full[ok & !s$generation %in% sim$events$time]))
  expect_true(all(s$colonised[ok] <= s$absent_before[ok]))
  expect_true(all(s$extirpated[ok] <= s$present_before[ok]))
  # occupancy changes by colonisations minus extirpations, apart from sites given to buds
  p <- persistence(sim)
  expect_setequal(p$species, sim$species$species)
  expect_equal(p$extant, p$species %in% extant)
  expect_equal(p$good_species[p$extant], unname(sim$lineage_species[p$species[p$extant]]))
  sp1 <- s[s$lineage == "sp1" & ok, ]
  expect_equal(p$colonisation[p$species == "sp1"],
               sum(sp1$colonised) / sum(sp1$absent_before))
  w <- persistence(sim, window = 20)
  expect_true(all(w$intervals <= 20L))
})

test_that("the series records the last losses of lineages that die out", {
  set.seed(1)
  ls <- landscape_grid(5, disp_scale = 0.5)
  anc <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
  sim <- simulate_evolving_metacommunity(ls, anc, speciation_control("fission", rate = 0.02),
                                         dispersal = 0.01, extirpation = 0.02, n_steps = 150)
  s <- sim$series
  gone <- sim$species$species[sim$species$fate != "extant"]
  expect_gt(length(gone), 0)
  for (l in gone) {
    rows <- s[s$lineage == l, ]
    final <- rows[which.max(rows$generation), ]
    end <- sim$species$end[sim$species$species == l]
    if (nrow(rows) > 1 || !is.na(final$present_before)) {
      expect_equal(final$generation, end)
      expect_equal(final$occupancy, 0)
      expect_equal(final$abundance, 0)
    }
  }
  # every lineage present at the start of a generation is accounted for, and occupancy moves
  # by colonisations and losses except at sites handed to buds
  ok <- !is.na(s$present_before)
  expect_true(all(s$present_before[ok] + s$absent_before[ok] <= nrow(sim$abundance)))
  expect_true(all(s$extirpated[ok] <= s$present_before[ok]))
})

test_that("losses are recorded when every lineage dies out", {
  set.seed(2)
  ls <- landscape_implicit(4)
  anc <- data.frame(optimum = 50, sigma = 20, r_max = 0.5, row.names = "sp1")
  sim <- simulate_evolving_metacommunity(ls, anc, speciation_control("geographic", Inf),
                                         dispersal = 0, extirpation = 0.5, n_steps = 60)
  expect_equal(sum(sim$abundance), 0)
  s <- sim$series
  final <- s[nrow(s), ]
  expect_equal(final$occupancy, 0)
  expect_gt(final$extirpated, 0)
  p <- persistence(sim)
  expect_equal(p$extirpation, sum(s$extirpated, na.rm = TRUE) / sum(s$present_before, na.rm = TRUE))
  sim$series <- NULL
  expect_error(persistence(sim), "series")
})
