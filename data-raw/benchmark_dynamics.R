# Time of simulate_metacommunity() and where it goes. Run from the root of the package.
suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(1)
run <- function(species, sites, steps, grid = FALSE) {
  tree <- treeshapesim::sim_tree_shape(species)
  land <- if (grid) landscape_grid(side = round(sqrt(sites))) else landscape_implicit(sites)
  seconds <- system.time(simulate_metacommunity(land, species_traits(tree),
                                                competition_matrix(tree, "phylo"), dispersal = 0.01,
                                                n_steps = steps))[["elapsed"]]
  data.frame(species, sites, steps, grid, seconds, ms_per_step = 1000 * seconds / steps)
}
print(rbind(run(36, 20, 300), run(50, 100, 2200), run(50, 100, 2200, TRUE), run(200, 100, 500),
            run(50, 400, 500, TRUE), run(200, 400, 300, TRUE)))
tree <- treeshapesim::sim_tree_shape(50)
profile <- tempfile()
Rprof(profile)
invisible(simulate_metacommunity(landscape_implicit(100), species_traits(tree),
                                 competition_matrix(tree, "phylo"), n_steps = 2200))
Rprof(NULL)
print(head(summaryRprof(profile)$by.self))
