# a known history of areas: root 5 -> node 6 (edge 1) -> a (1), b (1); root -> c (2)
known_history <- function() {
  h <- ape::reorder.phylo(ape::read.tree(text = "((a:1,b:1):1,c:2);"), "cladewise")
  e <- paste(h$edge[, 1], h$edge[, 2])
  maps <- list("4 5" = c(A = 0.5, B = 0.5), "5 1" = c(B = 1), "5 2" = c(B = 0.4, A = 0.6),
               "4 3" = c(A = 2))
  h$maps <- unname(maps[e])
  class(h) <- c("simmap", "phylo")
  h
}

test_that("region_pool limits species to the sites of their region", {
  h <- known_history()
  pool <- region_pool(h, c("A", "A", "B"))
  expect_equal(dimnames(pool), list(c("site1", "site2", "site3"), c("a", "b", "c")))
  expect_equal(unname(pool[, "a"]), c(FALSE, FALSE, TRUE))
  expect_equal(unname(pool[, "c"]), c(TRUE, TRUE, FALSE))
  expect_equal(rownames(region_pool(h, c(x = "A", y = "B"))), c("x", "y"))
})

test_that("arrival_ages gives the time since each lineage entered its region", {
  h <- known_history()
  comm <- matrix(c(0, 1, 1, 1, 1, 0, 1, 0, 1), 3, dimnames = list(NULL, c("a", "b", "c")))
  ages <- arrival_ages(h, comm, c("A", "B", "A"))
  # a is in B since node 5 entered B at age 1.5; b entered A at age 0.6; c in A since the root
  expect_equal(ages[2, "a"], 1.5)
  expect_equal(ages[1, "b"], 0.6)
  expect_equal(ages[1, "c"], 2)
  expect_true(is.na(ages[1, "a"]))   # absent
  expect_true(is.na(ages[3, "a"]))   # present, but its tip is in B and site 3 in A
  expect_equal(ages[3, "c"], 2)
  # not ultrametric: the time from the tip to the present is added
  short <- h
  short$edge.length[short$edge[, 2] == 3] <- 1.5
  short$maps[[which(short$edge[, 2] == 3)]] <- c(A = 1.5)
  expect_equal(arrival_ages(short, comm, c("A", "B", "A"))[1, "c"], 2)
})

test_that("as_herodotools gives the inputs of Herodotools with the true ancestral areas", {
  h <- known_history()
  comm <- matrix(c(0, 1, 0, 1, 0, 0, 1, 0, 0), 3, dimnames = list(NULL, c("a", "b", "c")))
  expect_warning(input <- as_herodotools(comm, c("A", "B", "A"), history = h),
                 "sites without species dropped: site3")
  expect_named(input, c("W", "tree", "biogeo", "ancestral.area"))
  expect_equal(input$ancestral.area$area, c("A", "B"))
  expect_equal(rownames(input$ancestral.area), c("N4", "N5"))
  expect_equal(input$biogeo$region, c("A", "B"))
  expect_equal(class(input$tree), "phylo")
  expect_null(input$tree$maps)
  expect_equal(unname(input$W[2, ]), c(1, 0, 0))
  without <- as_herodotools(comm[1:2, ], c("A", "B"), tree = ape::as.phylo(h))
  expect_null(without$ancestral.area)
  expect_error(as_herodotools(comm, c("A", "B", "A")), "no tree")
  expect_error(as_herodotools(comm[, 1:2], c("A", "B", "A"), history = h), "same species")
  expect_error(as_herodotools(comm[1:2, ], c("A", "AB"), history = h), "contained")
  expect_error(region_pool(ape::rcoal(3), "A"), "simmap")
})

test_that("a simulation with a region pool goes through as_herodotools", {
  skip_if_not_installed("treeshapesim")
  set.seed(8)
  tree <- treeshapesim::sim_tree_shape(10)
  areas <- trait_mk(tree, q = 0.5, states = c("A", "B"))
  regions <- rep(c("A", "B"), each = 6)
  sim <- simulate_metacommunity(landscape_implicit(12), species_traits(tree),
                                competition_matrix(tree, "none"), n_steps = 100,
                                pool = region_pool(areas, regions))
  input <- suppressWarnings(as_herodotools(sim, regions, history = areas))
  occurring <- input$W == 1
  tips <- areas$states[colnames(input$W)]
  expect_true(all(outer(input$biogeo$region, tips, "==")[occurring]))
  ages <- arrival_ages(areas, sim, regions)
  expect_true(all(is.finite(ages[sim$abundance > 0])))
  expect_true(all(ages[is.finite(ages)] <= max(ape::node.depth.edgelength(tree))))
})
