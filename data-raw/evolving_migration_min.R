# Calibration of migration_min, the migration rate that links two populations in the cluster
# rule of speciation_control("geographic"). Design: that of data-raw/evolving_peripheral.R (60
# sites at random positions, uniform environment, one ancestor, dispersal 0.02, disp_scale 0.1,
# threshold 150, 800 generations), with migration_min from 0.001 to 0.05, and the gradual rule
# as reference. Reports extant species, buddings, extinctions and the Spearman rho between
# connectivity and endemism. Run from the root of the package; REPS changes the size.
# Result (5 replicates, 2026-11-12): the strongest link of each population had a median rate of
# 0.006 (10%-90%: 0.0035-0.013) with about 80 individuals per population. With migration_min
# 0.001 and 0.002 all populations stayed linked and nothing split; at 0.005, 75 extant species
# (88 buddings, mean range 26 sites); at 0.01, 315; at 0.02 and 0.05, 343 and 363, close to
# the gradual rule (365), as almost no link holds. The default was set to 0.005, inside the
# transition; rho between connectivity and endemism ranged from -0.12 to +0.09.
suppressMessages(devtools::load_all(".", quiet = TRUE))
reps <- as.integer(Sys.getenv("REPS", "5"))
ancestor <- data.frame(optimum = 50, sigma = 30, r_max = 5, row.names = "sp1")
grid <- rbind(expand.grid(divergence = "cluster", migration_min = c(0.001, 0.002, 0.005, 0.01,
                                                                    0.02, 0.05),
                          rep = seq_len(reps), stringsAsFactors = FALSE),
              expand.grid(divergence = "gradual", migration_min = NA, rep = seq_len(reps),
                          stringsAsFactors = FALSE))
run <- function(i) {
  g <- grid[i, ]
  set.seed(3000 + g$rep)
  xy <- cbind(stats::runif(60), stats::runif(60))
  ls <- landscape_points(xy, env = rep(50, 60), disp_scale = 0.1)
  sp <- if (g$divergence == "cluster") {
    speciation_control("geographic", 150, migration_min = g$migration_min)
  } else {
    speciation_control("geographic", 150, divergence = "gradual", homogenization = 0.5)
  }
  sim <- simulate_evolving_metacommunity(ls, ancestor, sp, dispersal = 0.02, n_steps = 800,
                                         max_species = 5000)
  kernel <- exp(-as.matrix(stats::dist(xy)) / 0.1)
  diag(kernel) <- 0
  occ <- sim$abundance > 0
  endemic <- occ[, colSums(occ) == 1, drop = FALSE]
  share <- ifelse(rowSums(occ) > 0, rowSums(endemic) / rowSums(occ), NA_real_)
  data.frame(g, species = ncol(occ), splits = nrow(sim$events),
             extinct = sum(sim$species$fate == "extinct"),
             range = mean(colSums(occ)),
             rho = suppressWarnings(stats::cor(rowSums(kernel), share, method = "spearman",
                                               use = "complete.obs")))
}
started <- Sys.time()
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), run,
                                         mc.cores = as.integer(Sys.getenv("CORES", "6"))))
res$migration_min[is.na(res$migration_min)] <- -1
summary_table <- stats::aggregate(cbind(species, splits, extinct, range, rho) ~ divergence +
                                    migration_min, res,
                                  function(x) round(mean(x, na.rm = TRUE), 3),
                                  na.action = na.pass)
saveRDS(list(results = res, summary = summary_table, reps = reps,
             minutes = as.numeric(difftime(Sys.time(), started, units = "mins"))),
        "data-raw/evolving_migration_min.rds")
print(summary_table)
