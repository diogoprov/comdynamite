ancestor <- function(optimum = 50, sigma = 20) {
  data.frame(optimum = optimum, sigma = sigma, r_max = 5, row.names = "sp1")
}
# two clusters of sites with no dispersal between them
two_islands <- function(per_island = 4, env = c(50, 50)) {
  xy <- cbind(lon = c(rep(0, per_island), rep(10, per_island)),
              lat = c(seq(0, 0.03, length.out = per_island), seq(0, 0.03, length.out = per_island)))
  landscape_points(xy, env = rep(env, each = per_island), disp_scale = 1, lonlat = TRUE)
}

test_that("without a split, the simulation is simulate_metacommunity()", {
  ls <- landscape_grid(5, disp_scale = 0.5)
  ev <- evolution_control(genetic_var = 10, ne_ratio = 0.5)
  for (imm in c("poisson", "multinomial")) {
    set.seed(1)
    a <- simulate_evolving_metacommunity(ls, ancestor(), speciation_control("geographic", Inf),
                                         evolution = ev, dispersal = 0.01, n_steps = 120,
                                         extirpation = 0.01, immigration = imm)
    set.seed(1)
    b <- simulate_metacommunity(ls, ancestor(), matrix(0.05, dimnames = list("sp1", "sp1")),
                                dispersal = 0.01, n_steps = 120, extirpation = 0.01,
                                seed_steps = numeric(0), evolution = ev, immigration = imm)
    expect_identical(unname(a$abundance), unname(b$abundance))
    expect_equal(unname(a$genetic), unname(b$evolution$genetic))
    expect_equal(nrow(a$species), 1)
    expect_null(a$tree)
  }
})

test_that("isolated groups of populations bud when their divergence reaches the threshold", {
  ls <- two_islands()
  set.seed(2)
  sim <- simulate_evolving_metacommunity(ls, ancestor(), speciation_control("geographic", 30),
                                         dispersal = 0.05, n_steps = 60,
                                         initial = matrix(50L, 8, 1,
                                                          dimnames = list(NULL, "sp1")))
  expect_equal(sim$events$time, 30)
  expect_equal(sim$events$parent, "sp1")
  expect_equal(sim$events$new, "sp2")
  # the mother keeps one island and its age; the bud takes the other
  expect_equal(sort(sim$species$species[sim$species$fate == "extant"]), c("sp1", "sp2"))
  expect_equal(sim$species$birth, c(0, 30))
  occupied <- lapply(c("sp1", "sp2"), function(s) which(sim$abundance[, s] > 0))
  expect_setequal(lapply(occupied, function(x) unique(x > 4)), list(FALSE, TRUE))
  expect_true(ape::is.ultrametric(sim$tree))
  expect_equal(unname(ape::branching.times(sim$tree)), 60 - 30)
  # within each island, populations exchange migrants and never split
  expect_equal(nrow(sim$events), 1)
})

test_that("several groups at once bud in order, with internal branches of length zero", {
  xy <- cbind(lon = rep(c(0, 10, 20), each = 3), lat = rep(c(0, 0.01, 0.02), 3))
  ls <- landscape_points(xy, env = rep(50, 9), disp_scale = 1, lonlat = TRUE)
  set.seed(10)
  sim <- simulate_evolving_metacommunity(ls, ancestor(), speciation_control("geographic", 20),
                                         dispersal = 0.05, n_steps = 40,
                                         initial = matrix(50L, 9, 1,
                                                          dimnames = list(NULL, "sp1")))
  expect_equal(nrow(sim$events), 1)
  expect_equal(sim$events$new, "sp2,sp3")
  expect_equal(sim$species$order, c(0L, 1L, 2L))
  expect_equal(ape::Ntip(sim$tree), 3)
  expect_true(ape::is.binary(sim$tree))
  expect_equal(sort(unname(ape::branching.times(sim$tree))), c(20, 20))
  expect_true(any(sim$tree$edge.length == 0))
  # a record at the step of the budding shows the new species
  set.seed(10)
  rec <- simulate_evolving_metacommunity(ls, ancestor(), speciation_control("geographic", 20),
                                         dispersal = 0.05, n_steps = 20, record_steps = 20,
                                         initial = matrix(50L, 9, 1,
                                                          dimnames = list(NULL, "sp1")))
  expect_equal(unname(rec$history[, , 1]), unname(rec$abundance))
})

test_that("pruning keeps the emigrants aligned with the alive species", {
  st <- new.env()
  st$k <- 2
  st$m <- 2
  st$alive <- 1:2
  st$n <- matrix(c(0L, 0L, 5L, 5L), 2, 2)
  st$at <- diag(2)
  st$em <- matrix(c(0, 0, 3, 4), 2, 2)
  st$dispersal <- st$extirpation <- c(0.1, 0.1)
  expect_equal(evolving_prune(st), 1L)
  expect_equal(st$em, matrix(c(3, 4), 2, 1))
})

test_that("gradual divergence lets strong gene flow hold a species together", {
  ls <- landscape_grid(4, disp_scale = 1)
  run <- function(h) {
    set.seed(3)
    simulate_evolving_metacommunity(ls, ancestor(),
                                    speciation_control("geographic", 20, divergence = "gradual",
                                                       homogenization = h),
                                    dispersal = 0.05, n_steps = 100)
  }
  expect_equal(nrow(run(10)$events), 0)
  expect_gt(nrow(run(0)$events), 0)
})

test_that("ecological speciation follows divergence of genetic means despite gene flow", {
  ls <- two_islands(env = c(10, 90))
  ev <- evolution_control(genetic_var = 50, ne_ratio = Inf, founder = FALSE)
  set.seed(4)
  sim <- simulate_evolving_metacommunity(ls, ancestor(sigma = 30),
                                         speciation_control("ecological", 30), evolution = ev,
                                         dispersal = 0.01, n_steps = 100,
                                         initial = matrix(50L, 8, 1,
                                                          dimnames = list(NULL, "sp1")))
  expect_equal(nrow(sim$events), 1)
  means <- sim$traits$optimum
  expect_gt(diff(range(means)), 30)
})

test_that("the trees and tables agree", {
  set.seed(5)
  ls <- landscape_grid(6, disp_scale = 0.5)
  ev <- evolution_control(genetic_var = 20, ne_ratio = 0.5)
  sim <- simulate_evolving_metacommunity(ls, ancestor(),
                                         speciation_control("geographic", 60,
                                                            migration_min = 0.002),
                                         evolution = ev, dispersal = 0.005, n_steps = 300,
                                         extirpation = 0.005, max_species = 2000)
  sp <- sim$species
  expect_equal(nrow(sp) - 1, sum(lengths(strsplit(sim$events$new, ","))))
  extant <- sp$species[sp$fate == "extant"]
  expect_setequal(colnames(sim$abundance), extant)
  expect_true(all(colSums(sim$abundance) > 0))
  expect_setequal(sim$tree$tip.label, extant)
  expect_setequal(sim$tree_full$tip.label, sp$species)
  expect_true(ape::is.ultrametric(sim$tree))
  expect_equal(max(ape::node.depth.edgelength(sim$tree_full)) + sim$tree_full$root.edge, 300)
  expect_equal(length(sim$richness), 300)
  expect_equal(dim(sim$alpha), rep(length(extant), 2))
  x <- as_metacommunity(sim)
  expect_s3_class(x, "mc_sim")
  expect_s3_class(as_mc_data(x), "mc_data")
  expect_output(print(sim), "mc_evolving")
})

test_that("several species at the start hang from their tree", {
  set.seed(6)
  tr <- ape::read.tree(text = "(a:100,b:100);")
  traits <- data.frame(optimum = c(45, 55), sigma = 20, r_max = 5, row.names = c("a", "b"))
  attr(traits, "tree") <- tr
  ls <- two_islands()
  sim <- simulate_evolving_metacommunity(ls, traits, speciation_control("geographic", 30),
                                         dispersal = 0.05, n_steps = 50, competition = "none",
                                         initial = matrix(50L, 8, 2,
                                                          dimnames = list(NULL, c("a", "b"))))
  expect_equal(max(ape::node.depth.edgelength(sim$tree_full)), 150)
  expect_true(all(c("a", "b") %in% sim$species$species))
  expect_error(simulate_evolving_metacommunity(ls, rbind(ancestor(), sp2 = ancestor()),
                                               speciation_control("geographic", 30)), "tree")
})

test_that("settings are checked", {
  ls <- landscape_grid(3)
  expect_error(speciation_control("geographic"), "threshold")
  expect_output(print(speciation_control("geographic", 10)), "cluster")
  expect_output(print(speciation_control("geographic", 10, divergence = "gradual")),
                "homogenization")
  expect_output(print(speciation_control("ecological", 10)), "ecological")
  expect_error(simulate_evolving_metacommunity(ls, ancestor(), list()), "speciation_control")
  expect_error(simulate_evolving_metacommunity(ls, ancestor(),
                                               speciation_control("ecological", 10)),
               "genetic variance")
  expect_error(simulate_evolving_metacommunity(ls, ancestor(),
                                               speciation_control("geographic", 10),
                                               evolution = list()), "evolution_control")
  set.seed(7)
  expect_error(simulate_evolving_metacommunity(landscape_grid(5, disp_scale = 0.2), ancestor(),
                                               speciation_control("geographic", 2),
                                               dispersal = 0.001, n_steps = 50,
                                               max_species = 3), "max_species")
})

test_that("edge cases of names, extinction, records and trees are handled", {
  ls <- two_islands()
  start <- function(nm) matrix(50L, 8, 1, dimnames = list(NULL, nm))
  # a species at the start named like a daughter would be
  anc2 <- ancestor()
  rownames(anc2) <- "sp2"
  set.seed(8)
  sim <- simulate_evolving_metacommunity(ls, anc2, speciation_control("geographic", 30),
                                         dispersal = 0.05, n_steps = 60, initial = start("sp2"))
  expect_equal(anyDuplicated(sim$species$species), 0)
  expect_equal(nrow(sim$events), 1)
  # heavy extirpation: every extant species has individuals, richness has one value per
  # generation, and the table agrees with it
  for (seed in 1:10) {
    set.seed(seed)
    s <- simulate_evolving_metacommunity(landscape_grid(4, disp_scale = 0.3), ancestor(),
                                         speciation_control("geographic", 5),
                                         dispersal = 0.001, extirpation = 0.15, n_steps = 30)
    expect_true(all(colSums(s$abundance) > 0))
    expect_length(s$richness, 30)
    expect_equal(s$richness[30], ncol(s$abundance))
    expect_false(anyNA(s$traits$optimum))
  }
  # records reach the output and the colonisation history
  set.seed(9)
  s <- simulate_evolving_metacommunity(ls, ancestor(), speciation_control("geographic", 30),
                                       dispersal = 0.05, n_steps = 60, record_steps = c(20, 60),
                                       initial = start("sp1"))
  expect_equal(dim(s$history), c(8, ncol(s$abundance), 2))
  expect_s3_class(colonisation_history(as_metacommunity(s)), "data.frame")
  expect_error(simulate_evolving_metacommunity(ls, ancestor(),
                                               speciation_control("geographic", 30),
                                               evolution = evolution_control(
                                                 0, steps_per_generation = 4),
                                               n_steps = 50), "multiple")
  traits <- data.frame(optimum = c(45, 55), sigma = 20, r_max = 5, row.names = c("a", "b"))
  attr(traits, "tree") <- ape::read.tree(text = "(a:100,b:50);")
  expect_error(simulate_evolving_metacommunity(ls, traits, speciation_control("geographic", 30),
                                               n_steps = 5), "ultrametric")
  attr(traits, "tree") <- ape::read.tree(text = "(a,b);")
  expect_error(simulate_evolving_metacommunity(ls, traits, speciation_control("geographic", 30),
                                               n_steps = 5), "branch lengths")
  attr(traits, "tree") <- ape::read.tree(text = "(a:1,c:1);")
  expect_error(simulate_evolving_metacommunity(ls, traits, speciation_control("geographic", 30),
                                               n_steps = 5), "tips")
})

test_that("migrants between sites of implicit landscapes are e / (m - 1)", {
  st <- new.env()
  st$em <- matrix(c(10, 10, 0, 0, 0), 5, 1)
  st$k <- 1
  st$alive <- 1
  flow <- evolving_flow(st, list(matrix(0, 5, 5)), NULL)[[1]]
  expect_equal(flow[1, 2], 10 / 4)
  expect_equal(flow[3, 1], 0)
  expect_equal(flow[1, 3], 10 / 4)
  expect_equal(unname(diag(flow)), rep(0, 5))
})

test_that("flow without emigrants is unchanged", {
  st <- new.env()
  st$em <- NULL
  expect_identical(evolving_flow(st, list(1), NULL), list(1))
})

test_that("divergence within a cluster is frozen and resets after reset_generations", {
  st <- new.env()
  st$k <- 1
  st$alive <- 1
  st$n <- matrix(c(100, 100), 2, 1)
  st$disp_matrix <- NULL
  d <- matrix(c(0, 25, 25, 0), 2, 2)
  tg <- matrix(0, 2, 2)
  sp <- speciation_control("geographic", 30, migration_min = 0.01, reset_generations = 3)
  linked <- matrix(c(0, 5, 5, 0), 2, 2)
  apart <- matrix(c(0, 0.1, 0.1, 0), 2, 2)
  step <- function(d, tg, f) evolving_divergence(st, list(d), list(tg), list(f), sp)
  one <- step(d, tg, linked)
  expect_equal(one$divergence[[1]][1, 2], 25)
  expect_equal(one$together[[1]][1, 2], 1)
  two <- step(one$divergence[[1]], one$together[[1]], linked)
  expect_equal(two$divergence[[1]][1, 2], 25)
  three <- step(two$divergence[[1]], two$together[[1]], linked)
  expect_equal(three$divergence[[1]][1, 2], 0)
  # a generation apart restarts the count and adds to divergence
  back <- step(two$divergence[[1]], two$together[[1]], apart)
  expect_equal(back$divergence[[1]][1, 2], 26)
  expect_equal(back$together[[1]][1, 2], 0)
})

test_that("ecological groups need a lasting difference with reduced gene flow", {
  st <- new.env()
  st$k <- 1
  st$alive <- 1
  st$n <- matrix(c(100, 100), 2, 1)
  st$z <- matrix(c(10, 40), 2, 1)
  st$disp_matrix <- NULL
  sp <- speciation_control("ecological", 20, min_generations = 2, migration_max = 0.01)
  d <- matrix(0, 2, 2)
  low <- matrix(c(0, 0.1, 0.1, 0), 2, 2)
  high <- matrix(c(0, 50, 50, 0), 2, 2)
  one <- evolving_divergence(st, list(d), list(d), list(low), sp)$divergence[[1]]
  expect_equal(one[1, 2], 1)
  expect_length(evolving_groups(st, 1, one, sp), 0)
  two <- evolving_divergence(st, list(one), list(d), list(low), sp)$divergence[[1]]
  expect_length(evolving_groups(st, 1, two, sp), 2)
  # gene flow above migration_max restarts the count
  reset <- evolving_divergence(st, list(one), list(d), list(high), sp)$divergence[[1]]
  expect_equal(reset[1, 2], 0)
})

test_that("drift alone does not split species in a uniform environment", {
  ls <- landscape_grid(4, disp_scale = 1)
  ev <- evolution_control(genetic_var = 50, ne_ratio = 0.05)
  run <- function(sp, n_steps) {
    set.seed(11)
    simulate_evolving_metacommunity(ls, ancestor(), sp, evolution = ev, dispersal = 0.05,
                                    n_steps = n_steps)
  }
  expect_equal(nrow(run(speciation_control("ecological", 5), 100)$events), 0)
  # without the guard, drift and founder sampling split it at once
  expect_gt(nrow(run(speciation_control("ecological", 5, min_generations = 1,
                                        migration_max = Inf), 10)$events), 0)
})

test_that("competition falls with the difference of phenotypes in each site", {
  st <- new.env()
  st$alive <- 1:2
  st$ev <- list(sigma = c(10, 10))
  st$at <- matrix(0.05, 2, 2)
  hook <- niche_crowding(st)
  n <- matrix(c(10, 0, 10, 5), 2, 2)
  same <- hook(n, matrix(50, 2, 2))
  expect_equal(same[n > 0], (n %*% st$at)[n > 0])
  apart <- hook(n, matrix(c(50, 50, 80, 80), 2, 2))
  expect_lt(apart[1, 1], same[1, 1])
  expect_equal(apart[1, 1], 0.05 * 10 + 0.05 * 10 * exp(-30^2 / 400))
})

test_that("point mutation starts new species from single individuals", {
  set.seed(1)
  ls <- landscape_grid(4, disp_scale = 0.5)
  run <- function(rate) {
    set.seed(12)
    simulate_evolving_metacommunity(ls, ancestor(), speciation_control("point", rate = rate),
                                    dispersal = 0.05, n_steps = 200, max_species = 2000)
  }
  # rate 0 draws nothing: the same as no speciation
  set.seed(12)
  none <- simulate_evolving_metacommunity(ls, ancestor(), speciation_control("geographic", Inf),
                                          dispersal = 0.05, n_steps = 200)
  expect_identical(run(0)$abundance, none$abundance)
  sim <- run(5e-4)
  expect_gt(nrow(sim$events), 0)
  expect_equal(nrow(sim$species) - 1, sum(lengths(strsplit(sim$events$new, ","))))
  expect_true(ape::is.ultrametric(sim$tree))
  expect_setequal(sim$tree$tip.label, colnames(sim$abundance))
  expect_true(all(colSums(sim$abundance) > 0))
  expect_equal(max(ape::node.depth.edgelength(sim$tree)) + sim$tree$root.edge, 200)
  # mothers emptied by mutants go extinct
  set.seed(1)
  hot <- simulate_evolving_metacommunity(ls, ancestor(), speciation_control("point", rate = 0.003),
                                         dispersal = 0.05, n_steps = 200, max_species = 5000)
  expect_true(all(colSums(hot$abundance) > 0))
  expect_equal(tail(hot$richness, 1), ncol(hot$abundance))
  # a mutant takes one individual from its mother
  st <- new.env()
  st$m <- 2
  st$k <- 1
  st$k_all <- 1
  st$alive <- 1
  st$species <- "sp1"
  st$n <- matrix(c(5L, 3L), 2, 1)
  for (nm in c("arrival", "propagule")) assign(nm, matrix(1L, 2, 1), st)
  for (nm in c("immigrants", "load")) assign(nm, matrix(0, 2, 1), st)
  st$present <- matrix(TRUE, 2, 1)
  st$z <- st$ph <- st$pool_mean <- matrix(50, 2, 1)
  st$history <- st$z_history <- NULL
  st$dispersal <- st$extirpation <- 0.1
  st$ev <- list(g = 1, b = 0, ve = 0, optimum = 50, sigma = 10, r_max = 5)
  st$fixed_z <- FALSE
  evolving_split(st, 1, list(2L, 2L), c("sp2", "sp3"), list(1L, 1L))
  expect_equal(st$n, matrix(c(5L, 1L, 0L, 1L, 0L, 1L), 2, 3))
  st$step <- 7
  evolving_split(st, 1, list(2L), "sp4", list(1L), mutants = TRUE)
  expect_equal(st$arrival[2, 4], 7L)
  expect_equal(st$propagule[2, 4], 1L)
  expect_equal(st$n[2, 1], 0L)
  expect_false(st$present[2, 1])
  expect_equal(st$z[2, 2:4], c(50, 50, 50))
  expect_true(all(is.na(st$z[1, 2:3])))
})

test_that("new species can shift their genetic means", {
  st <- new.env()
  st$m <- 2
  st$k <- 1
  st$k_all <- 1
  st$alive <- 1
  st$species <- "sp1"
  st$n <- matrix(c(5L, 3L), 2, 1)
  for (nm in c("arrival", "propagule")) assign(nm, matrix(1L, 2, 1), st)
  for (nm in c("immigrants", "load")) assign(nm, matrix(0, 2, 1), st)
  st$present <- matrix(TRUE, 2, 1)
  st$z <- st$ph <- st$pool_mean <- matrix(c(40, 60), 2, 1)
  st$history <- st$z_history <- NULL
  st$dispersal <- st$extirpation <- 0.1
  st$ev <- list(g = 1, b = 0, ve = 0, optimum = 50, sigma = 10, r_max = 5)
  st$fixed_z <- FALSE
  set.seed(13)
  evolving_split(st, 1, list(1:2), "sp2", jump_var = 100)
  shift <- st$z[1, 2] - 40
  expect_gt(abs(shift), 0)
  expect_equal(st$z[2, 2] - 60, shift)
  expect_equal(st$ev$optimum[2] - 50, shift)
  expect_equal(st$n[, 1], c(0L, 0L))
})

test_that("fission splits species into spatial or random groups", {
  # a chain of sites: the spatial group is a run of neighbours
  xy <- cbind(seq(0, 0.09, by = 0.01), 0)
  ls <- landscape_points(xy, env = rep(50, 10), disp_scale = 0.01)
  for (seed in 1:20) {
    set.seed(seed)
    g <- evolving_partition(1:10, "spatial", ls$dispersal)
    expect_equal(sort(c(g[[1]], g[[2]])), 1:10)
    expect_true(length(g[[1]]) > 0 && length(g[[2]]) > 0)
    expect_true(all(diff(g[[2]]) == 1))
    r <- evolving_partition(c(2, 5, 7), "random", NULL)
    expect_setequal(c(r[[1]], r[[2]]), c(2, 5, 7))
    expect_true(length(r[[1]]) > 0 && length(r[[2]]) > 0)
  }
  set.seed(14)
  sim <- simulate_evolving_metacommunity(landscape_grid(4, disp_scale = 0.5), ancestor(),
                                         speciation_control("fission", rate = 0.05),
                                         dispersal = 0.05, n_steps = 100, max_species = 2000)
  expect_gt(nrow(sim$events), 0)
  expect_true(ape::is.binary(sim$tree))
  expect_true(ape::is.ultrametric(sim$tree))
})

test_that("incipient lineages belong to the nearest good species", {
  set.seed(1)
  sp <- data.frame(species = paste0("sp", 1:5), parent = c(NA, "sp1", "sp2", "sp1", "sp4"),
                   birth = c(0, 10, 20, 30, 40), end = c(50, 50, 50, 45, 50),
                   fate = c("extant", "extant", "extant", "extinct", "extant"),
                   completion = c(0, 60, 25, 35, 70),
                   status = c("good", "incipient", "good", "good", "incipient"),
                   stringsAsFactors = FALSE)
  obs <- evolving_observed(sp)
  expect_equal(unname(obs$group[c("sp1", "sp2", "sp3", "sp5")]), c("sp1", "sp1", "sp3", "sp4"))
  # sp4 is extinct, so its incipient daughter sp5 represents it
  expect_equal(obs$rep, c(sp1 = "sp1", sp3 = "sp3", sp4 = "sp5"))
  ls <- landscape_grid(4, disp_scale = 0.5)
  run <- function(completion, dist = "fixed") {
    set.seed(15)
    simulate_evolving_metacommunity(ls, ancestor(),
                                    speciation_control("point", rate = 5e-4,
                                                       completion = completion,
                                                       completion_dist = dist),
                                    dispersal = 0.05, n_steps = 200, max_species = 2000)
  }
  inst <- run(0)
  never <- run(1000)
  expect_identical(nrow(never$species), nrow(inst$species))
  expect_equal(colnames(never$abundance), "sp1")
  expect_equal(unname(never$abundance[, 1]), unname(rowSums(inst$abundance)))
  expect_null(never$tree)
  expect_true(all(never$richness_good == 1))
  some <- run(30, "exponential")
  expect_output(print(some), "good species")
  expect_true(all(some$lineage_species %in% colnames(some$abundance)))
  expect_equal(sum(some$abundance), sum(some$lineages))
  if (!is.null(some$tree)) expect_setequal(some$tree$tip.label, colnames(some$abundance))
  # an incipient lineage keeps the species of its parent at its birth, even if the parent
  # completes later
  sp2 <- data.frame(species = c("sp1", "sp2", "sp3"), parent = c(NA, "sp1", "sp2"),
                    birth = c(0, 10, 20), end = 50, fate = "extant",
                    completion = c(0, 30, 80), status = c("good", "good", "incipient"),
                    stringsAsFactors = FALSE)
  expect_equal(unname(evolving_observed(sp2)$group["sp3"]), "sp1")
  expect_s3_class(as_mc_data(as_metacommunity(some)), "mc_data")
  # records sum the lineages of each species
  set.seed(15)
  rec <- simulate_evolving_metacommunity(ls, ancestor(),
                                         speciation_control("point", rate = 5e-4,
                                                            completion = 1000),
                                         dispersal = 0.05, n_steps = 200, record_steps = 200,
                                         max_species = 2000)
  expect_gt(ncol(rec$lineages), 1)
  expect_equal(unname(rec$history[, 1, 1]), unname(rowSums(rec$lineages)))
})

test_that("trees with many buds on one lineage are built", {
  q <- 3000
  sp <- data.frame(species = paste0("sp", seq_len(q + 1)), parent = c(NA, rep("sp1", q)),
                   birth = c(0, seq_len(q) / 10), end = q / 10 + 1,
                   fate = "extant", order = c(0L, rep(1L, q)), stringsAsFactors = FALSE)
  tr <- evolving_tree(sp, NULL)
  expect_equal(ape::Ntip(tr), q + 1)
  expect_true(ape::is.ultrametric(tr))
})

test_that("settings of the new modes are checked", {
  expect_error(speciation_control("point"), "rate")
  expect_error(speciation_control("point", rate = 2))
  expect_output(print(speciation_control("point", rate = 1e-3, jump_var = 4, completion = 10)),
                "completion")
  expect_output(print(speciation_control("fission", rate = 0.1)), "spatial")
})

test_that("split times of several buddings at once agree with one at a time", {
  st0 <- matrix(c(0, -5, -5, 0), 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
  many <- evolving_split_time(st0, c("a", "a", "b"), c("c", "d", "e"), 10)
  expect_equal(many["c", "a"], 10)
  expect_equal(many["c", "d"], 10)
  expect_equal(many["c", "b"], -5)
  expect_equal(many["c", "e"], -5)
  expect_equal(many["e", "b"], 10)
  expect_equal(many, t(many))
  expect_equal(unname(diag(many)), rep(0, 5))
})
