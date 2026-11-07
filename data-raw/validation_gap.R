# The default `gap` of patch_permutations() and rigid_permutations() with sparse sites, with
# communities that follow a hidden autocorrelated environment (kernel bandwidth of 3 sites)
# and areas without effect. Compares the former default of patch_permutations() (1.5 times
# the median distance to the nearest site, `patch_nn`) with the current one (1.5 times the
# median distance to the fourth nearest site, `patch`), and evaluates rigid_permutations()
# with its default, which also follows the fourth nearest site, and with `independent = 10`.
# The rule was chosen on areas in opposite corners with 25% and 40% of the cells of the
# square (tuning); the other scenarios were held out. 150 replicates and 99 permutations.
# About 25 min on 6 cores. Run from the root of the package, with the package installed.
source("data-raw/validation_helpers.R")
n_rep <- as.integer(Sys.getenv("N_REP", "150"))
cores <- as.integer(Sys.getenv("CORES", "6"))
out_file <- Sys.getenv("OUT", "data-raw/validation_gap.rds")
n_perm <- 99

landscapes <- list(
  sparse25 = function() sparse_sites(square, 0.25),
  sparse40 = function() sparse_sites(square, 0.4),
  sparse60 = function() sparse_sites(square, 0.6),
  sparse80 = function() sparse_sites(square, 0.8),
  strip40 = function() sparse_sites(strip, 0.4),
  uniform30 = function() cbind(runif(120, 0.5, 20.5), runif(120, 0.5, 20.5))
)
# a scenario is a landscape and a layout, as in "sparse25_corners"
scenario_data <- function(scenario) {
  layout <- sub("^.*_", "", scenario)
  land <- sub("_[a-z]+$", "", scenario)
  if (land == "clusters") {
    cl <- cluster_sites()
    xy <- cl$xy
  } else {
    xy <- landscapes[[land]]()
  }
  event <- if (layout == "far") {
    far_clusters(xy, cl$cluster)
  } else if (layout %in% names(fixed)) {
    fixed[[layout]](xy)
  } else {
    get(layout)(xy)
  }
  list(xy = xy, event = event)
}
# the former default of patch_permutations()
patch_nn <- function(xy) {
  d <- as.matrix(dist(xy))
  1.5 * median(apply(d + diag(Inf, nrow(d)), 1, min))
}
methods <- list(
  patch_nn = list(method = "patch", gap = NULL),
  patch = list(method = "patch"),
  rigid = list(method = "rigid"),
  rigid_i10 = list(method = "rigid", independent = 10)
)
# patch_nn takes its gap from the landscape of each replicate
method_perms <- function(spec, xy, events, n_perm) {
  if (spec$method == "patch" && "gap" %in% names(spec) && is.null(spec$gap)) {
    spec$gap <- patch_nn(xy)
  }
  tryCatch(suppressWarnings(
    if (spec$method == "patch") {
      do.call(patch_permutations, c(list(xy, events, n_perm), spec[setdiff(names(spec), "method")]))$event
    } else {
      do.call(rigid_permutations, c(list(xy, events, n_perm), spec[setdiff(names(spec), "method")]))$event
    }), error = function(e) NULL)
}
sc <- function(scenario) data.frame(scenario = scenario, effect = 0, bandwidth = 3)
design <- rbind(
  # tuning
  sc(c("sparse25_corners", "sparse40_corners")),
  # held out
  sc(c("sparse40_tails", "sparse60_adjacent", "sparse80_blocks", "clusters_tails"))
)
cat(nrow(design), "scenarios\n")
rows <- run_design(design, scenario_data, methods, n_rep, n_perm, cores)
rates <- rate_table(rows, names(methods))
print(rates, digits = 3)

# the gap of each rule on one landscape of each kind
set.seed(1)
gaps <- t(sapply(c(names(landscapes), "clusters"), function(land) {
  xy <- if (land == "clusters") cluster_sites()$xy else landscapes[[land]]()
  d <- as.matrix(dist(xy))
  c(sites = nrow(xy), patch_nn = patch_nn(xy),
    patch = 1.5 * comdynamite:::neighbour_distance(d, 4),
    rigid = comdynamite:::neighbour_distance(d, 4))
}))
print(round(gaps, 2))
saveRDS(list(rates = rates, gaps = gaps), out_file)
