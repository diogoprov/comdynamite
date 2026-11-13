# Three mechanisms missing from the engine, tested on the departures of the replication of
# Urban & De Meester (2009) (data-raw/replication_urban2009.R, article "Replicating Urban &
# De Meester (2009)"):
# * lottery recruitment into the microsites vacated by mortality (`recruitment = "lottery"`,
#   1000 microsites per site), in place of Beverton-Holt growth;
# * overlapping generations: a lottery with baseline mortality 0.1 per step (their b), so one
#   generation is 10 steps, against a lottery with discrete generations (mortality 1);
# * a genetic variance built from new mutations and eroded by selection: `mutation_var` with
#   `variance_selection = TRUE` (the Bulmer reduction G^2 / (sigma^2 + G + VE) per generation,
#   permanent as under their clonal inheritance), against mutation and drift alone.
# The mechanisms are crossed (3 recruitment settings x 2 variance settings) on the three
# experiments, reduced: the colonisation delay (Fig. 2, delays of -10 to 10 generations, 100
# runs), a 6 x 6 sub-grid of the dispersal grid (Fig. 4, 50 runs) and the gene flow back into
# the sources (Fig. 3, 20 runs of 2000 generations). Delays, run lengths and records are in
# generations; with overlapping generations they are converted to steps by run_urban().

source("data-raw/replication_urban2009.R")

mechanisms <- data.frame(
  recruitment = c("beverton_holt", "lottery", "lottery"),
  mortality = c(1, 1, 0.1),
  label = c("Beverton-Holt", "lottery, discrete", "lottery, overlapping")
)
mech_grid <- function(genetic, vsel = c(FALSE, TRUE)) {
  out <- merge(mechanisms, expand.grid(genetic = genetic, variance_selection = vsel,
                                       stringsAsFactors = FALSE))
  # without genetic variance there is nothing to erode
  out[!(out$genetic == "none" & out$variance_selection), ]
}

if (sys.nframe() == 0) {
  cores <- max(1, parallel::detectCores() - 2)
  run_all <- function(design, f) {
    out <- parallel::mclapply(seq_len(nrow(design)), function(i) f(design[i, ]),
                              mc.cores = cores, mc.preschedule = FALSE)
    bad <- vapply(out, inherits, logical(1), "try-error")
    if (any(bad)) stop("runs failed: ", out[[which(bad)[1]]])
    do.call(rbind, out)
  }
  t0 <- Sys.time()
  # slow runs (overlapping generations) first, so that the cores stay busy
  by_cost <- function(d) d[order(d$mortality), ]

  d2 <- merge(mech_grid(c("none", "mutation_sp1", "mutation_both", "standing_sp1",
                          "standing_both")),
              expand.grid(delay = c(-10, -5, -2, -1, 0, 1, 2, 5, 10), rep = seq_len(100)))
  fig2 <- run_all(by_cost(d2), function(d) {
    sim <- run_urban(d$genetic, d$delay, seed = 10000 * d$rep + d$delay + 50,
                     recruitment = d$recruitment, mortality = d$mortality,
                     variance_selection = d$variance_selection)
    data.frame(d, n1 = sim$abundance[3, 1], n2 = sim$abundance[3, 2])
  })
  message("Fig. 2: ", format(Sys.time() - t0))

  m_sub <- c(0.005, 0.0126, 0.0315, 0.0792, 0.199, 0.5)
  d4 <- merge(mech_grid(c("none", "mutation_both")),
              expand.grid(m1 = m_sub, m2 = m_sub, rep = seq_len(50)))
  fig4 <- run_all(by_cost(d4), function(d) {
    sim <- run_urban(d$genetic, 0, c(d$m1, d$m2),
                     seed = 20000 * d$rep + 100 * match(d$m1, m_sub) + match(d$m2, m_sub),
                     recruitment = d$recruitment, mortality = d$mortality,
                     variance_selection = d$variance_selection)
    data.frame(d, n1 = sim$abundance[3, 1], n2 = sim$abundance[3, 2],
               src1 = sim$abundance[1, 1], src2 = sim$abundance[2, 2])
  })
  message("Fig. 4: ", format(Sys.time() - t0))

  d3 <- merge(mech_grid("mutation_both"), expand.grid(m = c(0.01, 0.05, 0.10),
                                                      rep = seq_len(20)))
  fig3 <- run_all(by_cost(d3), function(d) {
    sim <- run_urban("mutation_both", 0, rep(d$m, 2), seed = 30000 * d$rep + 100 * d$m,
                     n_steps = 2000, record_steps = c(500, 1000, 2000),
                     recruitment = d$recruitment, mortality = d$mortality,
                     variance_selection = d$variance_selection)
    h <- sim$history
    z <- sim$evolution$history
    do.call(rbind, lapply(seq_len(dim(h)[3]), function(s) {
      data.frame(d, step = c(500, 1000, 2000)[s],
                 new1 = h[3, 1, s], new2 = h[3, 2, s], src1_sp1 = h[1, 1, s],
                 src1_sp2 = h[1, 2, s], src2_sp1 = h[2, 1, s], src2_sp2 = h[2, 2, s],
                 z_src1_sp1 = z[1, 1, s], z_src2_sp2 = z[2, 2, s])
    }))
  })
  message("Fig. 3: ", format(Sys.time() - t0))
  saveRDS(list(fig2 = fig2, fig4 = fig4, fig3 = fig3, m_sub = m_sub),
          "data-raw/mechanisms_urban2009.rds")
}
