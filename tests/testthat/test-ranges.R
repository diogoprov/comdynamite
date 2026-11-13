range_sim <- function(mode = "fission", seed = 1, ...) {
  set.seed(seed)
  ls <- landscape_grid(5, disp_scale = 0.5)
  anc <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
  sc <- switch(mode, fission = speciation_control("fission", rate = 0.02, ...),
               point = speciation_control("point", rate = 5e-4, ...),
               geographic = speciation_control("geographic", 60, ...))
  list(sim = simulate_evolving_metacommunity(ls, anc, sc, dispersal = 0.01, extirpation = 0.02,
                                             n_steps = 200),
       areas = ifelse(ls$coords$x < 2.5, "W", "E"))
}

test_that("the log of ranges rebuilds the presence of every lineage", {
  r <- range_sim()
  sim <- r$sim
  g <- sim$ranges
  expect_named(g, c("generation", "lineage", "site", "present"))
  # replaying the changes gives the final presence of the extant lineages and nothing for the
  # extinct ones
  pres <- matrix(FALSE, nrow(sim$abundance), nrow(sim$species),
                 dimnames = list(NULL, sim$species$species))
  for (i in seq_len(nrow(g))) pres[g$site[i], g$lineage[i]] <- g$present[i]
  extant <- colnames(sim$lineages)
  expect_equal(unname(pres[, extant]), unname(sim$lineages > 0))
  expect_false(any(pres[, setdiff(colnames(pres), extant)]))
  # changes alternate in each site
  key <- paste(g$lineage, g$site)
  alternate <- tapply(g$present, key, function(p) all(p[-1] != p[-length(p)]))
  expect_true(all(alternate))
  # new lineages start at their birth
  first <- tapply(g$generation, g$lineage, min)
  births <- stats::setNames(sim$species$birth, sim$species$species)
  expect_equal(as.vector(first), unname(births[names(first)]))
})

test_that("range_history() maps ranges of areas along the full and reconstructed trees", {
  for (mode in c("fission", "point", "geographic")) {
    r <- range_sim(mode)
    sim <- r$sim
    for (tt in c("full", "reconstructed")) {
      h <- range_history(sim, r$areas, tree = tt)
      expect_s3_class(h, "simmap")
      expect_equal(unname(rowSums(h$mapped.edge)), h$edge.length)
      expect_equal(vapply(h$maps, sum, numeric(1)), h$edge.length)
      expect_true(all(colnames(h$mapped.edge) %in% c("W", "E", "W+E")))
    }
    # tips of the full tree end in the ranges of the extant lineages
    h <- range_history(sim, r$areas, tree = "full")
    extant <- colnames(sim$lineages)
    final <- vapply(extant, function(l) {
      paste(intersect(c("W", "E"), unique(r$areas[sim$lineages[, l] > 0])), collapse = "+")
    }, character(1))
    expect_equal(history_states(h)$tips[extant], final)
    hr <- range_history(sim, r$areas)
    expect_setequal(hr$tip.label, colnames(sim$abundance))
  }
  # the starting range is that of the first population
  r <- range_sim()
  h <- range_history(r$sim, factor(r$areas, c("E", "W")), tree = "full")
  expect_true(all(colnames(h$mapped.edge) %in% c("W", "E", "E+W")))
  expect_error(range_history(r$sim, r$areas[-1]), "one area per site")
  expect_error(range_history(r$sim, sub("W", "W+", r$areas)), "\\+")
  expect_error(range_history(as_metacommunity(r$sim), r$areas), "simulate_evolving")
  named <- stats::setNames(rev(r$areas), rev(seq_along(r$areas)))
  expect_equal(range_history(r$sim, named)$maps, range_history(r$sim, r$areas)$maps)
  old <- r$sim
  old$ranges <- NULL
  expect_error(range_history(old, r$areas), "ranges")
})

test_that("a lineage that changes range shows it on its branch", {
  r <- range_sim("geographic")
  h <- range_history(r$sim, r$areas, tree = "full")
  changes <- sum(lengths(h$maps) > 1)
  expect_gt(changes, 0)
  # one area per site: the timeline of the root follows its occupancy of each half
  tl <- range_timelines(r$sim$ranges, factor(r$areas, c("W", "E")))
  expect_equal(tl$sp1$time[1], 0)
  expect_true(all(diff(tl$sp1$time) > 0))
  expect_true(all(tl$sp1$state[-1] != tl$sp1$state[-nrow(tl$sp1)]))
})

test_that("ranges of several areas feed the functions of histories", {
  r <- range_sim()
  sim <- r$sim
  h <- range_history(sim, r$areas)
  rates <- history_rates(h, "state", state_rates = c(W = 1, E = 2, "W+E" = 3))
  expect_length(rates, nrow(h$edge))
  s0 <- history_rates(h, "sympatry", r = -0.1)
  s1 <- history_rates(h, "sympatry", r = -0.1, overlap = TRUE)
  expect_true(all(s1 <= s0 + 1e-12))
  expect_true(any(s1 < s0))
  mc <- as_metacommunity(sim)
  regions <- r$areas
  input <- suppressWarnings(as_herodotools(sim, regions))
  expect_equal(input$ancestral.area$area, history_states(h)$nodes)
  expect_identical(suppressWarnings(as_herodotools(sim, regions, history = h)), input)
  pool <- region_pool(h, regions)
  tips <- history_states(h)$tips
  expect_equal(unname(pool[regions == "W", "sp1"]), rep(grepl("W", tips[["sp1"]]),
                                                        sum(regions == "W")))
  ages <- arrival_ages(h, mc, regions)
  present <- mc$abundance > 0
  expect_true(all(!is.na(ages[present])))
  expect_true(all(is.na(ages[!present])))
  expect_true(all(ages[present] >= 0))
})

test_that("range_history() follows the representative lineages with protracted speciation", {
  set.seed(3)
  ls <- landscape_grid(5, disp_scale = 0.5)
  anc <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
  sim <- simulate_evolving_metacommunity(ls, anc,
                                         speciation_control("fission", rate = 0.05,
                                                            completion = 50),
                                         dispersal = 0.01, extirpation = 0.02, n_steps = 200)
  r <- list(areas = ifelse(ls$coords$x < 2.5, "W", "E"))
  expect_false(is.null(sim$tree))
  expect_gt(length(sim$lineage_species), length(sim$tree$tip.label))
  h <- range_history(sim, r$areas)
  expect_setequal(h$tip.label, sim$tree$tip.label)
  expect_equal(unname(rowSums(h$mapped.edge)), sim$tree$edge.length)
  rep <- range_representatives(sim)
  full <- range_history(sim, r$areas, tree = "full")
  expect_equal(unname(history_states(h)$tips[names(rep)]),
               unname(history_states(full)$tips[rep]))
})
