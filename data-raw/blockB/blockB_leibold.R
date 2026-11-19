# Block B, Leibold et al. (2019): more replicates of the published design and sweeps of its
# central parameters, with the functions of data-raw/replication_leibold2019.R.
# * replicates: the 20 cells of the replication (6 G x 3 Nm and 2 controls) with 500
#   replicates; seeds 1000 + rep as in the replication, so its 10 replicates are a subset;
# * sweep_pext: G (6) x Nm (1, 3, 10, 30, 100) x patch extinction (0, 0.001, 0.005, 0.01,
#   0.05) x 100 replicates, whether local dominance and regional neutrality depend on the
#   values chosen;
# * sweep_vm: mutation from no standing variation, Vm (0.01, 0.05, 0.1, 0.5, 1) x Nm (5
#   values) x 100 replicates, the open question of Vm above 0.5.
# Summaries at generations 1000, 2500, 5000, 7500 and 10000 (summarise_leibold()), and the
# classes of patches at generation 5000. Run from the root of the package; see run_chunks.R
# for CORES, CHUNK, OUT and MAX_ROWS. N_REP and N_SWEEP change the replicates.
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("data-raw/replication_leibold2019.R")
source("data-raw/blockB/run_chunks.R")
n_rep <- as.integer(Sys.getenv("N_REP", "500"))
n_sweep <- as.integer(Sys.getenv("N_SWEEP", "100"))
keep_steps <- c(1000, 2500, 5000, 7500, 10000)

one <- function(d) {
  run <- run_leibold(d$G, d$landscape_type, seed = d$seed, nm = d$nm, p_ext = d$p_ext,
                     mutation_var = if (d$vm > 0) d$vm)
  idx <- match(keep_steps, run$steps)
  run$steps <- run$steps[idx]
  run$abundance <- run$abundance[, , idx, drop = FALSE]
  run$genetic <- run$genetic[, , idx, drop = FALSE]
  cls <- classify_leibold(run)
  data.frame(d[rep(1, length(idx)), ], summarise_leibold(run), preadapted = cls[["preadapted"]],
             switched = cls[["switched"]], shared = cls[["shared"]], row.names = NULL)
}

nm_grid <- c(1, 3, 10, 30, 100)
g_grid <- c(0, 20, 50, 100, 200, 400)
replicates <- rbind(
  expand.grid(G = g_grid, nm = c(1, 10, 100), landscape_type = "two", rep = seq_len(n_rep),
              stringsAsFactors = FALSE),
  expand.grid(G = c(0, 200), nm = 10, landscape_type = "one", rep = seq_len(n_rep),
              stringsAsFactors = FALSE))
replicates$p_ext <- 0.01
replicates$vm <- 0
replicates$seed <- 1000 + replicates$rep

sweep_pext <- expand.grid(G = g_grid, nm = nm_grid, p_ext = c(0, 0.001, 0.005, 0.01, 0.05),
                          landscape_type = "two", rep = seq_len(n_sweep),
                          stringsAsFactors = FALSE)
sweep_pext$vm <- 0
sweep_pext$seed <- 100000 + sweep_pext$rep

sweep_vm <- expand.grid(vm = c(0.01, 0.05, 0.1, 0.5, 1), nm = nm_grid, landscape_type = "two",
                        rep = seq_len(n_sweep), stringsAsFactors = FALSE)
sweep_vm$G <- 1e-6
sweep_vm$p_ext <- 0.01
sweep_vm$seed <- 200000 + sweep_vm$rep

run_chunks(replicates, one, "leibold_replicates")
run_chunks(sweep_pext, one, "leibold_sweep_pext")
run_chunks(sweep_vm, one, "leibold_sweep_vm")
message("Leibold block B done")
