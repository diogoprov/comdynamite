# Qualitative check of geographic speciation in simulate_evolving_metacommunity() against the
# prediction of Desjardins-Proulx & Gravel (2012, Am. Nat. 179: 137-144) that peripheral and
# poorly connected communities accumulate more species formed by speciation. Design (ours): 60
# sites at random positions in a square, a uniform environment (no selection gradient), one
# ancestor, geographic speciation with each rule of divergence. Connectivity of a site: the sum
# over the other sites of the dispersal kernel before normalisation, exp(-d / disp_scale).
# Endemism of a site: the share of the species present in it that occur in no other site. The
# check is the Spearman correlation between connectivity and endemism across sites, in each
# replicate. Run from the root of the package; REPS changes the size.
# Result (20 replicates, 2026-11-12, after the fixes of the audit): no support. Spearman rho
# between connectivity and endemism averaged -0.030 with the cluster rule (negative in 55% of
# replicates) and -0.001 with the gradual rule (35%), with about 365 extant species and 1100
# splits per run (before the fixes: -0.076 and -0.068). With threshold 300 and dispersal 0.05
# the cluster rule gave almost no speciation and the gradual rule a mean rho of +0.04. The
# prediction is not reproduced; endemism may be a poor proxy for the species formed in situ,
# which the simulation does not yet record by site, and the design differs from theirs.
suppressMessages(devtools::load_all(".", quiet = TRUE))
reps <- as.integer(Sys.getenv("REPS", "20"))
ancestor <- data.frame(optimum = 50, sigma = 30, r_max = 5, row.names = "sp1")
grid <- expand.grid(divergence = c("cluster", "gradual"), rep = seq_len(reps),
                    stringsAsFactors = FALSE)
run <- function(i) {
  g <- grid[i, ]
  set.seed(3000 + i)
  xy <- cbind(stats::runif(60), stats::runif(60))
  ls <- landscape_points(xy, env = rep(50, 60), disp_scale = 0.1)
  sim <- simulate_evolving_metacommunity(
    ls, ancestor, speciation_control("geographic", 150, divergence = g$divergence,
                                     homogenization = 0.5),
    dispersal = 0.02, n_steps = 800, max_species = 5000)
  d <- as.matrix(stats::dist(xy))
  kernel <- exp(-d / 0.1)
  diag(kernel) <- 0
  connectivity <- rowSums(kernel)
  occ <- sim$abundance > 0
  endemic <- occ[, colSums(occ) == 1, drop = FALSE]
  share <- ifelse(rowSums(occ) > 0, rowSums(endemic) / rowSums(occ), NA_real_)
  data.frame(g, species = ncol(occ), splits = nrow(sim$events),
             extinct = sum(sim$species$fate == "extinct"),
             rho = suppressWarnings(stats::cor(connectivity, share, method = "spearman",
                                               use = "complete.obs")))
}
started <- Sys.time()
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), run,
                                         mc.cores = as.integer(Sys.getenv("CORES", "6"))))
summary_table <- stats::aggregate(cbind(species, splits, extinct, rho) ~ divergence, res,
                                  function(x) round(mean(x, na.rm = TRUE), 3), na.action = na.pass)
summary_table$rho_negative <- tapply(res$rho < 0, res$divergence, mean, na.rm = TRUE)
saveRDS(list(results = res, summary = summary_table, reps = reps,
             minutes = as.numeric(difftime(Sys.time(), started, units = "mins"))),
        "data-raw/evolving_peripheral.rds")
print(summary_table)
