# Block B, Urban & De Meester (2009): more replicates of the published experiments and a sweep
# of the mutational variance and the environmental contrast, with the functions of
# data-raw/replication_urban2009.R.
# * fig2: 7 genetic regimes x delays -20 to 20 x 500 replicates;
# * fig4: dispersal of each species (11 x 11) x 500 replicates;
# * fig3: dispersal 0.01, 0.05, 0.10 x 500 replicates, 2000 generations;
# seeds as in the replication, so its replicates are a subset;
# * sweep: Vm (0.01, 0.025, 0.05, 0.1, 0.25) x environment of the new site (50, 55, 60, 65,
#   70, sources at 40) x 4 regimes (mutation and standing variation, one or both species) x
#   delays (-20 to 20 by 5) x 200 replicates: whether the advantage of the mutation-evolver
#   and the priority of standing variation depend on these values.
# Run from the root of the package; see run_chunks.R for CORES, CHUNK, OUT and MAX_ROWS.
# N_REP and N_SWEEP change the replicates.
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("data-raw/replication_urban2009.R")
source("data-raw/blockB/run_chunks.R")
n_rep <- as.integer(Sys.getenv("N_REP", "500"))
n_sweep <- as.integer(Sys.getenv("N_SWEEP", "200"))
m_grid <- signif(exp(seq(log(0.005), log(0.5), length.out = 11)), 3)

d2 <- expand.grid(genetic = names(genetics), delay = -20:20, rep = seq_len(n_rep),
                  stringsAsFactors = FALSE)
run_chunks(d2, function(d) {
  sim <- run_urban(d$genetic, d$delay, seed = 10000 * d$rep + d$delay + 50)
  data.frame(d, n1 = sim$abundance[3, 1], n2 = sim$abundance[3, 2])
}, "urban_fig2")

d4 <- expand.grid(m1 = m_grid, m2 = m_grid, rep = seq_len(n_rep))
run_chunks(d4, function(d) {
  sim <- run_urban("mutation_both", 0, c(d$m1, d$m2), seed = 20000 * d$rep +
                     100 * match(d$m1, m_grid) + match(d$m2, m_grid))
  data.frame(d, n1 = sim$abundance[3, 1], n2 = sim$abundance[3, 2],
             src1 = sim$abundance[1, 1], src2 = sim$abundance[2, 2])
}, "urban_fig4")

d3 <- expand.grid(m = c(0.01, 0.05, 0.10), rep = seq_len(n_rep))
run_chunks(d3, function(d) {
  sim <- run_urban("mutation_both", 0, rep(d$m, 2), seed = 30000 * d$rep + 100 * d$m,
                   n_steps = 2000, record_steps = c(500, 1000, 1500, 2000))
  h <- sim$history
  data.frame(d[rep(1, 4), ], step = c(500, 1000, 1500, 2000), new1 = h[3, 1, ],
             new2 = h[3, 2, ], src1_sp1 = h[1, 1, ], src1_sp2 = h[1, 2, ],
             src2_sp1 = h[2, 1, ], src2_sp2 = h[2, 2, ], row.names = NULL)
}, "urban_fig3")

sweep <- expand.grid(vm = c(0.01, 0.025, 0.05, 0.1, 0.25), env_new = c(50, 55, 60, 65, 70),
                     genetic = c("mutation_sp1", "mutation_both", "standing_sp1",
                                 "standing_both"),
                     delay = seq(-20, 20, 5), rep = seq_len(n_sweep), stringsAsFactors = FALSE)
sweep$seed <- 50000000 + seq_len(nrow(sweep))
run_chunks(sweep, function(d) {
  sim <- run_urban(d$genetic, d$delay, seed = d$seed, env_sites = c(40, 40, d$env_new),
                   vm_value = d$vm)
  data.frame(d, n1 = sim$abundance[3, 1], n2 = sim$abundance[3, 2])
}, "urban_sweep")
message("Urban block B done")
