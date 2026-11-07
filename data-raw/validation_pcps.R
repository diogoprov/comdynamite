# Checks pcps_test() against the package PCPS, then measures how often it rejects when the
# environment has no relation to the communities, with free permutations of sites and with
# torus_permutations(). Saves data-raw/validation_pcps.rds.
#
# The package PCPS is not a dependency of comdynamite: install it to run the first part.
suppressMessages(devtools::load_all(".", quiet = TRUE))
library(PCPS)
set.seed(20261107)

# 1. Same statistics as PCPS, on random communities ---------------------------------------
compare <- function(i) {
  tree <- treeshapesim::sim_tree_shape(8)
  env <- stats::runif(15)
  comm <- matrix(stats::rpois(15 * 8, 3) + rbinom(15 * 8, 1, 0.5), 15,
                 dimnames = list(paste0("s", 1:15), tree$tip.label))
  distances <- ape::cophenetic.phylo(tree)[colnames(comm), colnames(comm)]
  frame <- data.frame(env = env, row.names = rownames(comm))
  theirs <- c(
    matrix.p.sig(comm, distances, frame, FUN = FUN.ADONIS2.global, runs = 1, method.p = "bray",
                 sqrt.p = TRUE, formula = p.dist ~ env)$obs.statistic,
    pcps.sig(comm, distances, frame, method = "bray", squareroot = TRUE, FUN = FUN.GLM,
             choices = 1, runs = 1, formula = pcps.1 ~ env)$obs.statistic,
    pcps.sig(comm, distances, frame, method = "bray", squareroot = TRUE, FUN = FUN.GLM,
             choices = 2, runs = 1, formula = pcps.2 ~ env)$obs.statistic
  )
  max(abs(pcps_test(comm, tree, env, n_perm = 1)$tests$statistic - theirs))
}
statistics <- max(vapply(1:20, compare, numeric(1)))

# the same p-values, up to sampling error, on one community with a gradient
tree <- treeshapesim::sim_tree_shape(10)
traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100), rmax_range = c(0, 30))
env <- seq(0, 100, length.out = 25)
comm <- simulate_static(traits, env + stats::rnorm(25, 0, 25), incidence = FALSE)$abundance
occupied <- rowSums(comm) > 0
comm <- comm[occupied, colSums(comm) > 0]
env <- env[occupied]
distances <- ape::cophenetic.phylo(tree)[colnames(comm), colnames(comm)]
frame <- data.frame(env = env, row.names = rownames(comm))
ours <- pcps_test(comm, tree, env, n_perm = 999)$tests
pick <- function(x) c(p_site = x$p.site.shuffle, p_species = x$p.taxa.shuffle)
theirs <- rbind(
  whole = pick(matrix.p.sig(comm, distances, frame, FUN = FUN.ADONIS2.global, runs = 999,
                            method.p = "bray", sqrt.p = TRUE, formula = p.dist ~ env)),
  pcps1 = pick(pcps.sig(comm, distances, frame, method = "bray", squareroot = TRUE, FUN = FUN.GLM,
                        choices = 1, runs = 999, formula = pcps.1 ~ env)),
  pcps2 = pick(pcps.sig(comm, distances, frame, method = "bray", squareroot = TRUE, FUN = FUN.GLM,
                        choices = 2, runs = 999, formula = pcps.2 ~ env))
)
p_values <- data.frame(test = ours$test, ours = ours[c("p_site", "p_species")], theirs)

# 2. Rejections with autocorrelated fields ------------------------------------------------
side <- 10
xy <- as.matrix(expand.grid(x = seq_len(side), y = seq_len(side)))
field <- function(bandwidth, periodic) {
  dx <- abs(outer(xy[, 1], xy[, 1], "-"))
  dy <- abs(outer(xy[, 2], xy[, 2], "-"))
  if (periodic) {
    dx <- pmin(dx, side - dx)
    dy <- pmin(dy, side - dy)
  }
  z <- drop(exp(-(dx^2 + dy^2) / (2 * bandwidth^2)) %*% stats::rnorm(nrow(xy)))
  (z - mean(z)) / stats::sd(z)
}
replicate_once <- function(i, bandwidth, periodic, relation) {
  tree <- treeshapesim::sim_tree_shape(20)
  traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100), rmax_range = c(0, 30))
  hidden <- field(bandwidth, periodic)
  env <- relation * hidden + sqrt(1 - relation^2) * field(bandwidth, periodic)
  scaled <- (hidden - min(hidden)) / diff(range(hidden)) * 100
  comm <- simulate_static(traits, scaled, incidence = TRUE)$abundance
  occupied <- rowSums(comm) > 0
  moves <- torus_permutations(side, n_perm = 199)
  free <- pcps_test(comm[occupied, ], tree, env[occupied], n_perm = 199)$tests
  torus <- pcps_test(comm, tree, env, n_perm = 199, site_permutations = moves)$tests
  data.frame(test = free$test, free_site = free$p_site, free_max = free$p_max,
             torus_site = torus$p_site, torus_max = torus$p_max, rep = i)
}
scenarios <- rbind(
  expand.grid(bandwidth = c(1.5, 3), periodic = c(TRUE, FALSE), relation = 0, n_rep = 300),
  expand.grid(bandwidth = c(1.5, 3), periodic = c(TRUE, FALSE), relation = 0.7, n_rep = 100)
)
runs <- lapply(seq_len(nrow(scenarios)), function(s) {
  run <- scenarios[s, ]
  out <- parallel::mclapply(seq_len(run$n_rep), function(i) {
    tryCatch(replicate_once(i, run$bandwidth, run$periodic, run$relation),
             error = function(e) NULL)
  }, mc.cores = 6, mc.set.seed = TRUE)
  message(s, " of ", nrow(scenarios), ": ", sum(lengths(out) > 0), " replicates")
  cbind(run[c("bandwidth", "periodic", "relation")], do.call(rbind, out), row.names = NULL)
})
runs <- do.call(rbind, runs)
methods <- c("free_site", "free_max", "torus_site", "torus_max")
rates <- stats::aggregate(runs[methods], runs[c("relation", "periodic", "bandwidth", "test")],
                          function(p) round(100 * mean(p <= 0.05, na.rm = TRUE), 1))
rates$n_rep <- stats::aggregate(rep ~ relation + periodic + bandwidth + test, runs, length)$rep
rates <- rates[order(rates$relation, rates$test, !rates$periodic, rates$bandwidth), ]

saveRDS(list(statistics = statistics, p_values = p_values, rates = rates),
        "data-raw/validation_pcps.rds")
print(statistics)
print(p_values)
print(rates, row.names = FALSE)
