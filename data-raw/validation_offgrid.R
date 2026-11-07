# Type I error and power of rigid_permutations() and patch_permutations() with coordinates off
# a regular lattice, with communities that follow a hidden autocorrelated environment (kernel
# bandwidth of 3 sites) and areas without effect or, with EFFECT 20, moving the environment of
# their sites. Landscapes: the square of 20 by 20 with jitter of up to 0.3 on each coordinate
# (used while building the method, with the adjacent and blocks layouts), 400 sites placed
# at random in the square, the irregular band of 436 sites with jitter, and 300 sites placed
# at random in a crescent (held out). 150 replicates and 99 permutations. About 30 min on
# 6 cores. Run from the root of the package, with the package installed.
source("data-raw/validation_helpers.R")
n_rep <- as.integer(Sys.getenv("N_REP", "150"))
cores <- as.integer(Sys.getenv("CORES", "6"))
out_file <- Sys.getenv("OUT", "data-raw/validation_offgrid.rds")
n_perm <- 99

landscapes <- list(
  jitter = function() square + matrix(runif(800, -0.3, 0.3), ncol = 2),
  uniform = function() cbind(runif(400, 0.5, 20.5), runif(400, 0.5, 20.5)),
  band_jitter = function() band + matrix(runif(2 * nrow(band), -0.3, 0.3), ncol = 2),
  crescent = function() {
    xy <- cbind(runif(3000, -12, 12), runif(3000, -12, 12))
    r <- sqrt(rowSums(xy^2))
    inner <- sqrt((xy[, 1] - 4)^2 + xy[, 2]^2)
    xy[r <= 12 & inner >= 7, ][1:300, ]
  }
)
# a scenario is a landscape and a layout, as in "band_jitter_blocks"
scenario_data <- function(scenario) {
  layout <- sub("^.*_", "", scenario)
  land <- sub("_[a-z]+$", "", scenario)
  xy <- landscapes[[land]]()
  event <- if (layout %in% names(fixed)) fixed[[layout]](xy) else get(layout)(xy)
  list(xy = xy, event = event)
}
methods <- list(
  patch = list(method = "patch"),
  rigid = list(method = "rigid")
)
sc <- function(scenario, effect = 0) data.frame(scenario = scenario, effect = effect, bandwidth = 3)
design <- rbind(
  sc(c("jitter_adjacent", "jitter_blocks")),
  sc(c("uniform_adjacent", "uniform_blocks")),
  sc(c("band_jitter_blocks", "crescent_blocks")),
  sc("jitter_adjacent", effect = 20)
)
cat(nrow(design), "scenarios\n")
rows <- run_design(design, scenario_data, methods, n_rep, n_perm, cores)
rates <- rate_table(rows, names(methods))
print(rates, digits = 3)

# placements and time of 999 permutations on about 750 sites off a lattice: 750 sites at
# random in a rectangle of 60 by 40, two round areas of radius 6, and the same with jitter
# on a lattice of 30 by 25
set.seed(1)
timing <- lapply(list(uniform750 = cbind(runif(750, 0, 60), runif(750, 0, 40)),
                      jitter750 = as.matrix(expand.grid(1:30, 1:25)) + matrix(runif(1500, -0.3, 0.3), ncol = 2)),
                 function(xy) {
  centres <- rbind(c(15, 15), c(35, 25)) * c(max(xy[, 1]) / 60, max(xy[, 2]) / 40)
  d1 <- sqrt((xy[, 1] - centres[1, 1])^2 + (xy[, 2] - centres[1, 2])^2)
  d2 <- sqrt((xy[, 1] - centres[2, 1])^2 + (xy[, 2] - centres[2, 2])^2)
  event <- cbind(event = ifelse(d1 < 6, 1, ifelse(d2 < 6, -1, 0)))
  t0 <- Sys.time()
  drawn <- suppressWarnings(rigid_permutations(xy, event, 999))
  data.frame(sites = nrow(xy), event_sites = sum(event != 0), placements = nrow(drawn$event),
             seconds = round(as.numeric(Sys.time() - t0, units = "secs")))
})
timing <- do.call(rbind, timing)
print(timing)
saveRDS(list(rates = rates, timing = timing), out_file)
