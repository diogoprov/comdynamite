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

test_that("isolated groups of populations split when their divergence reaches the threshold", {
  ls <- two_islands()
  set.seed(2)
  sim <- simulate_evolving_metacommunity(ls, ancestor(), speciation_control("geographic", 30),
                                         dispersal = 0.05, n_steps = 60,
                                         initial = matrix(50L, 8, 1,
                                                          dimnames = list(NULL, "sp1")))
  expect_equal(sim$events$time, 30)
  expect_equal(sim$events$parent, "sp1")
  expect_equal(sort(sim$species$species[sim$species$fate == "extant"]), c("sp2", "sp3"))
  occupied <- lapply(c("sp2", "sp3"), function(s) which(sim$abundance[, s] > 0))
  expect_setequal(lapply(occupied, function(x) unique(x > 4)), list(FALSE, TRUE))
  expect_true(ape::is.ultrametric(sim$tree))
  expect_equal(unname(ape::branching.times(sim$tree)), 60 - 30)
  # within each island, populations exchange migrants and never split
  expect_equal(nrow(sim$events), 1)
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
  sim <- simulate_evolving_metacommunity(ls, ancestor(), speciation_control("geographic", 60),
                                         evolution = ev, dispersal = 0.005, n_steps = 300,
                                         extirpation = 0.005, max_species = 2000)
  sp <- sim$species
  expect_equal(nrow(sp) - 1, sum(lengths(strsplit(sim$events$daughters, ","))))
  extant <- sp$species[sp$fate == "extant"]
  expect_setequal(colnames(sim$abundance), extant)
  expect_true(all(colSums(sim$abundance) > 0))
  expect_setequal(sim$tree$tip.label, extant)
  expect_setequal(sim$tree_full$tip.label, sp$species[sp$fate != "split"])
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
  expect_equal(flow[1, 2], 10 / 4 + 10 / 4)
  expect_equal(flow[1, 3], 10 / 4)
  expect_equal(unname(diag(flow)), rep(0, 5))
})

test_that("flow without emigrants is unchanged", {
  st <- new.env()
  st$em <- NULL
  expect_identical(evolving_flow(st, list(1), NULL), list(1))
})
