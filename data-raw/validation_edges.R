# Type I error and power of biogeographic events in node_test() in landscapes with edges. The
# community follows a hidden autocorrelated environment. Areas are placed at random and, with
# effect > 0, move the environment that sites offer. Compares free permutations,
# patch_permutations() and, in the square, permutations of a torus, which ignore the edges.
library(comdynamite)

n_rep <- as.integer(Sys.getenv("N_REP", "200"))
n_perm <- 199

grid <- as.matrix(expand.grid(1:70, 1:30))
middle <- 15 + 9 * sin(grid[, 1] / 9)
half <- 3 + 1.5 * cos(grid[, 1] / 5)
landscapes <- list(square = as.matrix(expand.grid(1:20, 1:20)),
                   strip = as.matrix(expand.grid(1:50, 1:8)),
                   band = unname(grid[abs(grid[, 2] - middle) <= half, ]))
torus <- torus_permutations(20)

smooth_field <- function(xy, bandwidth) {
  d <- as.matrix(dist(xy))
  drop(exp(-d^2 / (2 * bandwidth^2)) %*% rnorm(nrow(xy)))
}

# two squares of 25 sites around sites drawn at random, which do not overlap
blocks <- function(xy) {
  repeat {
    centres <- xy[sample(nrow(xy), 2), ]
    one <- pmax(abs(xy[, 1] - centres[1, 1]), abs(xy[, 2] - centres[1, 2])) <= 2
    two <- pmax(abs(xy[, 1] - centres[2, 1]), abs(xy[, 2] - centres[2, 2])) <= 2
    if (!any(one & two) && sum(one) >= 15 && sum(two) >= 15) return(one * 1 - two * 1)
  }
}

# irregular areas of several patches from the tails of a smooth field
tails <- function(xy) {
  f <- smooth_field(xy, 2.5)
  limits <- quantile(f, c(0.1, 0.9))
  (f >= limits[2]) * 1 - (f <= limits[1]) * 1
}

# one area of several patches against all other sites
others <- function(xy) {
  f <- smooth_field(xy, 2.5)
  ifelse(f >= quantile(f, 0.9), 1, -1)
}

one_rep <- function(i, land, kind, effect) {
  set.seed(i)
  xy <- landscapes[[land]]
  tree <- treeshapesim::sim_tree_shape(12)
  hidden <- smooth_field(xy, 3)
  hidden <- (hidden - min(hidden)) / diff(range(hidden)) * 100
  event <- get(kind)(xy)
  traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100), rmax_range = c(0, 30))
  comm <- simulate_static(traits, hidden + effect * event, incidence = TRUE)$abundance
  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  if (ncol(comm) < 4) return(NULL)
  nodes <- node_contrasts(ape::keep.tip(tree, colnames(comm)))
  noise <- cbind(noise = rnorm(nrow(xy)))
  events <- cbind(event = event)
  p <- function(...) {
    node_test(comm, nodes, noise, events, n_perm = n_perm, ...)$biogeography$p_partial_site[, 1]
  }
  data.frame(rep = i, land = land, kind = kind, effect = effect, node = seq_len(ncol(nodes)),
             free = p(),
             patch = p(area_permutations = patch_permutations(xy, events, n_perm)),
             torus = if (land == "square") {
               p(site_permutations = torus[sample(nrow(torus), n_perm), ])
             } else {
               NA_real_
             })
}

design <- expand.grid(i = seq_len(n_rep), land = names(landscapes),
                      kind = c("blocks", "tails", "others"), effect = c(0, 20),
                      stringsAsFactors = FALSE)
results <- parallel::mclapply(seq_len(nrow(design)), function(k) {
  one_rep(design$i[k], design$land[k], design$kind[k], design$effect[k])
}, mc.cores = as.integer(Sys.getenv("CORES", "8")))
failed <- vapply(results, inherits, logical(1), "try-error")
if (any(failed)) stop(sum(failed), " replicates failed: ", results[failed][[1]])
saveRDS(do.call(rbind, results), Sys.getenv("OUT", "data-raw/validation_edges.rds"))
