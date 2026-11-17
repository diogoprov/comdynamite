# Identifiability pilot of issue #97: can monthly tadpole and calling-male counts in a
# design of 13 ponds and 15 species separate four community models (SS, NT, PE, PM), and how
# much do a second and a third field year add? Each run is summarised for 1, 2 and 3 field
# years, so the three designs share the same simulations. Output: data-raw/seasonal_pilot/pilot.rds.

suppressMessages(devtools::load_all(quiet = TRUE))
source("data-raw/seasonal_pilot/seasonal_model.R")

n_per_model <- as.integer(Sys.getenv("N_PER_MODEL", "3000"))
cores <- as.integer(Sys.getenv("CORES", max(1, parallel::detectCores() - 2)))
models <- c("SS", "NT", "PE", "PM")
design <- seasonal_design(n_ponds = 13, n_species = 15, seed = 1)   # one fixed design

RNGkind("L'Ecuyer-CMRG")
set.seed(97)
jobs <- rep(models, each = n_per_model)
runs <- parallel::mclapply(seq_along(jobs), function(j) {
  p <- seasonal_prior(jobs[j])
  sp <- seasonal_species(design$n_species)
  r <- seasonal_run(design, p, sp)
  s <- lapply(1:3, function(y) suppressWarnings(seasonal_summarise(r, design, y)))
  list(model = jobs[j], par = unlist(p[-1]), stats = s)
}, mc.cores = cores)

pilot <- list(design = design, runs = runs, n_per_model = n_per_model,
              date = Sys.Date())
saveRDS(pilot, "data-raw/seasonal_pilot/pilot.rds", compress = "xz")
