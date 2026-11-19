# Block B, Hajian-Forooshani & Chase (2026): more replicates of the published grid and a sweep
# of the pulse, with one_run() of data-raw/replication_hajian2026.R.
# * replicates: extent (1 and 5 to 100 patches) x dispersal (1e-5, 1e-4, 1e-3) x niche (0.2,
#   0.25, 0.3, 0.4) x 500 replicates, and the undisturbed control; new seeds, in separate
#   ranges for the grid and the control;
# * sweep: pulse duration (10, 50, 101, 200, 400 steps from step 500) x value (0.1 to 0.5) x
#   dispersal (1e-5 to 1e-1) x extent (1, 20, 40, 60, 80, 100) x 100 replicates, niche 0.25,
#   with alpha and gamma 100 steps after the end of the pulse and at step 1000.
# Each run keeps the changes at steps 601, 700, 800 and 1000 (replicates) or at the end of the
# pulse plus 100 and at 1000 (sweep), relative to step 499. Run from the root of the package;
# see run_chunks.R for CORES, CHUNK, OUT and MAX_ROWS. N_REP and N_SWEEP change the replicates.
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("data-raw/replication_hajian2026.R")
source("data-raw/blockB/run_chunks.R")
n_rep <- as.integer(Sys.getenv("N_REP", "500"))
n_sweep <- as.integer(Sys.getenv("N_SWEEP", "100"))

grid <- expand.grid(rep = seq_len(n_rep), n_disturbed = c(0, 1, seq(5, 100, 5)),
                    dispersal = c(1e-5, 1e-4, 1e-3), niche = c(0.2, 0.25, 0.3, 0.4))
grid$seed <- ifelse(grid$n_disturbed == 0, 3000000, 1000000) + seq_len(nrow(grid))
run_chunks(grid, function(d) {
  out <- one_run(d$dispersal, d$n_disturbed, d$niche, d$seed,
                 record = c(499, 601, 700, 800, 1000))
  cbind(d, out, row.names = NULL)
}, "hajian_replicates")

sweep <- expand.grid(duration = c(10, 50, 101, 200, 400), value = seq(0.1, 0.5, 0.1),
                     dispersal = 10^(-5:-1), n_disturbed = c(1, 20, 40, 60, 80, 100),
                     rep = seq_len(n_sweep))
sweep$seed <- 5000000 + seq_len(nrow(sweep))
run_chunks(sweep, function(d) {
  after <- 500 + d$duration - 1 + 100
  out <- one_run(d$dispersal, d$n_disturbed, 0.25, d$seed, duration = d$duration,
                 value = d$value, record = unique(c(499, after, 1000)))
  out$when <- ifelse(out$step == 1000 & after != 1000, "end", "after")
  if (after == 1000) out <- rbind(out, transform(out, when = "end"))
  cbind(d[rep(1, nrow(out)), ], out, row.names = NULL)
}, "hajian_sweep")
message("Hajian block B done")
