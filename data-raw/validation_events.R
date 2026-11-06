# Type I error of biogeographic events in node_test() when sites of an area are neighbours.
# The community follows a hidden autocorrelated environment and the areas are placed at
# random, so no event has an effect. Compares free and torus permutations.
library(comdynamite)

side <- 15
n_rep <- as.integer(Sys.getenv("N_REP", "200"))
n_perm <- 199
torus <- torus_permutations(side)

block <- function(width) {
  corner <- sample(side, 2)
  x <- (corner[1] + seq_len(width) - 2) %% side
  y <- (corner[2] + seq_len(width) - 2) %% side
  as.vector(outer(x + 1, side * y, `+`))
}

one_rep <- function(i) {
  set.seed(i)
  tree <- ape::rcoal(12)
  hidden <- landscape_grid(side, env_bandwidth = 3)$env
  traits <- species_traits(tree, lambda = 1, optima_range = c(min(hidden), max(hidden)),
                           rmax_range = c(0, 30))
  comm <- simulate_static(traits, hidden, incidence = TRUE)$abundance
  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  nodes <- node_contrasts(ape::keep.tip(tree, colnames(comm)))
  # two blocks that do not overlap
  repeat {
    one <- block(5)
    two <- block(5)
    if (!length(intersect(one, two))) break
  }
  event <- rep(0, side^2)
  event[one] <- 1
  event[two] <- -1
  events <- cbind(event = event)
  noise <- cbind(noise = rnorm(side^2))
  free <- node_test(comm, nodes, noise, events, n_perm = n_perm)
  restricted <- node_test(comm, nodes, noise, events, n_perm = n_perm,
                          site_permutations = torus[sample(nrow(torus), n_perm), ])
  p <- function(x) x$biogeography$p_partial_site[, "event"]
  data.frame(rep = i, node = seq_along(p(free)), free = p(free), torus = p(restricted))
}

results <- do.call(rbind, parallel::mclapply(seq_len(n_rep), one_rep,
                                             mc.cores = as.integer(Sys.getenv("CORES", "8"))))
saveRDS(results, Sys.getenv("OUT", "data-raw/validation_events.rds"))
