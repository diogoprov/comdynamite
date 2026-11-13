# Replication of the individual-based model of Urban & De Meester (2009, Proc R Soc B 276:
# 4129-4138, doi:10.1098/rspb.2009.1382) with the quantitative-genetic engine of comdynamite.
#
# Their model (Section 2 and ESM Table S1): three patches of 1000 microsites with lottery
# competition; two source patches with environment 0.4, one filled with each of two initially
# equivalent species adapted to it, and a vacant patch with environment 0.6; survival
# (1 - b) exp(-((z - e) / omega)^2) with baseline mortality b = 0.1 and omega = 0.5; birth
# rate r = 0.5; dispersal m = 0.05 to each of the other two patches in equal halves; a haploid
# trait from 20 additive loci (0 or 1, phenotype the mean) plus a non-genetic deviation of SD
# 1/40; microsite environments with SD 1/40; mutation 1e-4 per locus (1e-6 in ESM Fig. S3);
# evolution by mutation alone or from standing variation (alleles drawn with probability equal
# to the mean phenotype); 5000 time steps; colonisation of the new patch delayed for one species
# by -200 to 200 steps in steps of 5; 100 runs per setting.
#
# Translation used here (environment and traits multiplied by 100):
# * sites with environment 40, 40 and 60, dispersal split equally between the other two sites
#   (landscape_matrix()), species 1 filling site 1 and species 2 site 2 at the start;
# * a Gaussian niche with sigma = omega / sqrt(2) x 100 = 35.4, since their survival is
#   exp(-((z - e) / omega)^2); the non-genetic and microsite deviations (SD 2.5 each) enter as
#   V_E = 12.5;
# * Beverton-Holt growth with r_max = 5 and alpha = 0.004 for all pairs (K = 1000, equal
#   competition within and between species);
# * one step = one generation, taken as 10 of their time steps (1 / b, the lifespan without
#   selection), so 500 generations, delays of -20 to 20 generations, dispersal 0.05 per
#   generation (offspring disperse once);
# * the new site is opened to each species at its introduction time with `pool_from`, and is
#   then colonised by dispersal from the sources; the first species can enter at generation 21,
#   after 20 generations in the sources;
# * genetic variance per population with mutation and drift (`mutation_var`): mutation adds
#   V_m = L mu a^2 = 20 x 1e-4 x 5^2 = 0.05 per generation (a = 5, the effect of one allele);
#   `ne_ratio = 0.5` makes the loss 1 / N per generation, as for a haploid population of N, so
#   the equilibrium in a full site is 2 x 0.5 x 1000 x 0.05 = 50; evolution by mutation starts
#   without variance (1e-6), and standing variation starts at L p (1 - p) a^2 = 120 for
#   p = 0.4; a species that does not evolve has G = 0.
# * Two more scenarios keep a constant G = 120 (no mutation), as in the rest of the package.
#   With `mutation_var`, selection does not erode the variance, so standing variation persists
#   and carries a variance load (lower mean fitness), unlike in their model, where selection
#   on 20 loci erodes it.
#
# Not translated: loci, alleles bounded between 0 and 1, lottery recruitment, overlapping
# generations and survival selection, clonal inheritance (as such), individual microsites.
#
# Experiments, as in their figures:
# * Fig. 2: final abundance of each species in the new site against the delay, for no
#   evolution, species 1 alone evolving, and both evolving, by mutation or standing variation;
# * Fig. 4: both species evolving by mutation, the new site open to both from generation 21,
#   dispersal of each species on a grid of 11 values from 0.005 to 0.5, response the relative
#   proportional dominance (N1 - N2) / 1000 in the new site;
# * Fig. 3: both evolving by mutation, dispersal 0.01, 0.05 or 0.10, 2000 generations: does the
#   winner of the new site lose its own source site to the other species?

devtools::load_all(quiet = TRUE)

species <- c("sp1", "sp2")
env <- c(40, 40, 60)
landscape <- landscape_matrix(matrix(1, 3, 3), env = env)
traits <- data.frame(optimum = c(40, 40), sigma = 0.5 / sqrt(2) * 100, r_max = 5,
                     row.names = species)
alpha <- matrix(0.004, 2, 2, dimnames = list(species, species))
initial <- cbind(sp1 = c(1000L, 0L, 0L), sp2 = c(0L, 1000L, 0L))
vm <- 20 * 1e-4 * 5^2
genetics <- list(
  none = list(g = c(0, 0), vm = c(0, 0)),
  mutation_sp1 = list(g = c(1e-6, 0), vm = c(vm, 0)),
  mutation_both = list(g = c(1e-6, 1e-6), vm = c(vm, vm)),
  standing_sp1 = list(g = c(120, 0), vm = c(vm, 0)),
  standing_both = list(g = c(120, 120), vm = c(vm, vm)),
  constant_sp1 = list(g = c(120, 0), vm = NULL),
  constant_both = list(g = c(120, 120), vm = NULL)
)

# `recruitment`, `mortality` and `variance_selection` are for data-raw/mechanisms_urban2009.R;
# with overlapping generations (mortality < 1), delays, n_steps, start and record_steps stay in
# generations and are converted to steps (1 / mortality per generation)
run_urban <- function(genetic, delay = 0, dispersal = c(0.05, 0.05), seed, n_steps = 500,
                      start = 21, record_steps = NULL, recruitment = "beverton_holt",
                      mortality = 1, variance_selection = FALSE) {
  set.seed(seed)
  gs <- genetics[[genetic]]
  per <- 1 / mortality
  open <- round((start - 1) * per) + 1
  from <- cbind(sp1 = c(1, 1, open + round(max(0, -delay) * per)),
                sp2 = c(1, 1, open + round(max(0, delay) * per)))
  evo <- if (all(gs$g == 0)) NULL else {
    evolution_control(genetic_var = gs$g, env_var = 12.5, ne_ratio = 0.5,
                      mutation_var = gs$vm,
                      variance_selection = variance_selection && !is.null(gs$vm))
  }
  sim <- simulate_metacommunity(landscape, traits, alpha, dispersal = dispersal,
                                n_steps = round(n_steps * per), seed_steps = integer(0),
                                initial = initial, pool_from = from,
                                record_steps = if (!is.null(record_steps)) round(record_steps * per),
                                evolution = evo, recruitment = recruitment, microsites = 1000,
                                mortality = mortality)
  sim
}

if (sys.nframe() == 0) {
  n_rep <- 100
  cores <- max(1, parallel::detectCores() - 2)
  run_all <- function(design, f) {
    out <- parallel::mclapply(seq_len(nrow(design)), function(i) f(design[i, ]),
                              mc.cores = cores, mc.preschedule = TRUE)
    bad <- vapply(out, inherits, logical(1), "try-error")
    if (any(bad)) stop("runs failed: ", out[[which(bad)[1]]])
    do.call(rbind, out)
  }
  t0 <- Sys.time()
  # Fig. 2: delays (positive: species 1 enters first)
  d2 <- expand.grid(genetic = names(genetics), delay = -20:20, rep = seq_len(n_rep),
                    stringsAsFactors = FALSE)
  fig2 <- run_all(d2, function(d) {
    sim <- run_urban(d$genetic, d$delay, seed = 10000 * d$rep + d$delay + 50)
    data.frame(d, n1 = sim$abundance[3, 1], n2 = sim$abundance[3, 2])
  })
  message("Fig. 2: ", format(Sys.time() - t0))
  # Fig. 4: dispersal of each species
  m_grid <- signif(exp(seq(log(0.005), log(0.5), length.out = 11)), 3)
  d4 <- expand.grid(m1 = m_grid, m2 = m_grid, rep = seq_len(n_rep))
  fig4 <- run_all(d4, function(d) {
    sim <- run_urban("mutation_both", 0, c(d$m1, d$m2), seed = 20000 * d$rep +
                       100 * match(d$m1, m_grid) + match(d$m2, m_grid))
    data.frame(d, n1 = sim$abundance[3, 1], n2 = sim$abundance[3, 2],
               src1 = sim$abundance[1, 1], src2 = sim$abundance[2, 2])
  })
  message("Fig. 4: ", format(Sys.time() - t0))
  # Fig. 3: gene flow back into the sources
  d3 <- expand.grid(m = c(0.01, 0.05, 0.10), rep = seq_len(40))
  fig3 <- run_all(d3, function(d) {
    sim <- run_urban("mutation_both", 0, rep(d$m, 2), seed = 30000 * d$rep + 100 * d$m,
                     n_steps = 2000, record_steps = seq(50, 2000, 50))
    h <- sim$history
    z <- sim$evolution$history
    do.call(rbind, lapply(seq_len(dim(h)[3]), function(s) {
      data.frame(d, step = as.integer(sub("step", "", dimnames(h)[[3]][s])),
                 new1 = h[3, 1, s], new2 = h[3, 2, s], src1_sp1 = h[1, 1, s],
                 src1_sp2 = h[1, 2, s], src2_sp1 = h[2, 1, s], src2_sp2 = h[2, 2, s],
                 z_src1_sp1 = z[1, 1, s], z_src2_sp2 = z[2, 2, s])
    }))
  })
  message("Fig. 3: ", format(Sys.time() - t0))
  saveRDS(list(fig2 = fig2, fig4 = fig4, fig3 = fig3, m_grid = m_grid),
          "data-raw/replication_urban2009.rds")
}
