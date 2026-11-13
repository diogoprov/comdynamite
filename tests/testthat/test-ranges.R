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

test_that("nodes have the range of the mother before the split", {
  r <- range_sim("geographic")
  sim <- r$sim
  h <- range_history(sim, r$areas, tree = "full")
  st <- history_states(h)
  tl <- range_timelines(sim$ranges, factor(r$areas, c("W", "E")))
  depth <- round(ape::node.depth.edgelength(h) -
                   stats::median(ape::node.depth.edgelength(h)[seq_along(h$tip.label)] -
                                   sim$species$end[match(h$tip.label, sim$species$species)]), 6)
  # every budding: the node is the range of the mother just before it
  for (i in which(!is.na(sim$species$parent))) {
    mother <- sim$species$parent[i]
    g <- sim$species$birth[i]
    nodes <- which(abs(depth - g) < 1e-6)
    nodes <- nodes[nodes > ape::Ntip(h)]
    expect_true(range_before(tl[[mother]], g) %in% st$nodes[nodes - ape::Ntip(h)])
  }
  expect_equal(st$nodes[1], range_before(tl$sp1, min(sim$species$birth[-1])))
  # inner pieces of the reconstructed maps all have a duration
  hr <- range_history(sim, r$areas)
  inner <- unlist(lapply(hr$maps, function(m) m[-length(m)]))
  expect_true(all(inner > 0))
})

test_that("range_history() handles one lineage, starting trees and named factors", {
  set.seed(1)
  ls <- landscape_grid(4, disp_scale = 0.5)
  areas <- ifelse(ls$coords$x < 1.5, "W", "E")
  anc <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
  one <- simulate_evolving_metacommunity(ls, anc, speciation_control("geographic", Inf),
                                         dispersal = 0.01, n_steps = 30)
  expect_null(range_history(one, areas))
  h1 <- range_history(one, areas, tree = "full")
  expect_length(h1$maps, 1)
  expect_equal(sum(h1$maps[[1]]), h1$edge.length)
  # a starting tree of four species
  set.seed(2)
  tr <- ape::rcoal(4)
  tr$edge.length <- tr$edge.length * 10
  traits <- species_traits(tr)
  sim <- simulate_evolving_metacommunity(ls, traits, speciation_control("fission", rate = 0.02),
                                         dispersal = 0.01, extirpation = 0.02, n_steps = 150)
  full <- range_history(sim, areas, tree = "full")
  expect_equal(unname(rowSums(full$mapped.edge)), full$edge.length)
  extant <- colnames(sim$lineages)
  final <- vapply(extant, function(l) {
    paste(intersect(c("W", "E"), unique(areas[sim$lineages[, l] > 0])), collapse = "+")
  }, character(1))
  expect_equal(unname(history_states(full)$tips[extant]), unname(final))
  pieces <- unlist(lapply(full$maps, function(m) m[-length(m)]))
  expect_true(all(pieces > 1e-6))
  # names reorder factors too
  f <- factor(areas, c("E", "W"))
  named <- stats::setNames(rev(f), rev(seq_along(f)))
  expect_equal(range_history(sim, named, "full")$maps, range_history(sim, f, "full")$maps)
})

test_that("overlap in history_rates() counts lineages that share an area", {
  r <- range_sim("point")
  h <- range_history(r$sim, r$areas, tree = "full")
  brute <- function(history, r) {
    edge <- history$edge
    depth <- ape::node.depth.edgelength(history)
    pieces <- do.call(rbind, lapply(seq_len(nrow(edge)), function(e) {
      m <- history$maps[[e]]
      end <- depth[edge[e, 1]] + cumsum(unname(m))
      data.frame(edge = e, state = names(m), from = c(depth[edge[e, 1]], end[-length(end)]),
                 to = end, stringsAsFactors = FALSE)
    }))
    pieces <- pieces[pieces$to > pieces$from, ]
    breaks <- sort(unique(c(pieces$from, pieces$to)))
    integral <- numeric(nrow(edge))
    for (i in seq_len(length(breaks) - 1)) {
      alive <- which(pieces$from <= breaks[i] & pieces$to >= breaks[i + 1])
      parts <- strsplit(pieces$state[alive], "+", fixed = TRUE)
      n <- vapply(parts, function(p) sum(vapply(parts, function(q) any(p %in% q), TRUE)), 1)
      for (k in seq_along(alive)) {
        integral[pieces$edge[alive[k]]] <- integral[pieces$edge[alive[k]]] +
          exp(r * n[k]) * (breaks[i + 1] - breaks[i])
      }
    }
    integral
  }
  expect_equal(sympatry_integral(h, -0.1, overlap = TRUE), brute(h, -0.1))
})
