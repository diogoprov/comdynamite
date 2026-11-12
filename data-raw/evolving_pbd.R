# Validation 2 of speciation in simulate_evolving_metacommunity(): protracted speciation against
# the protracted birth-death model (PBD; Etienne & Rosindell 2012). Etienne & Rosindell show
# that protracted speciation makes reconstructed trees look like diversification slowed down
# toward the present. Design (ours): random fission at rate 0.01 per species per generation
# (the initiation rate of the PBD, for species with two or more populations), completion after
# an exponential time with mean tau (completion rate 1 / tau, the lambda of the PBD), 6 x 6 grid,
# dispersal 0.05, patch extirpation 0.01, all species equivalent, 600 generations. Checks: the
# gamma statistic (Pybus & Harvey 2000) of the reconstructed tree should fall as tau grows, and
# the PBD fitted by maximum likelihood (package PBD, GPL-2, used only here; equal extinction of
# good and incipient species) should give completion rates that fall as tau grows, near 1 / tau.
# Run from the root of the package; REPS changes the size.
# Result (10 replicates, 2026-11-12; medians): gamma fell from -0.98 (tau = 0) to -3.28 (tau =
# 25) and -4.58 (tau = 100), the apparent slowdown expected under protracted speciation; with
# tau = 0 it is already slightly negative, perhaps from competition limiting richness. The PBD
# fits converged in 3 of 10 trees with tau = 0 (the completion rate runs to infinity) and in all
# with tau > 0. With tau = 25 the estimated completion rate (0.045) was close to 1 / tau (0.04)
# and the initiation rate (0.0077) close to the realised rate of fission per lineage (0.0099).
# With tau = 100 the completion rate was underestimated (0.0018 against 0.01) and initiation and
# extinction inflated (0.050 and 0.040), with two fits far off (b near 28 and 82): with long
# completion times few completed nodes are left to inform the fit. The ordering of the
# completion rates follows tau.
# Etienne, R. S., & Rosindell, J. (2012). Systematic Biology, 61, 204-213.
# https://doi.org/10.1093/sysbio/syr091
# Pybus, O. G., & Harvey, P. H. (2000). Proceedings of the Royal Society B, 267, 2267-2272.
# https://doi.org/10.1098/rspb.2000.1278
suppressMessages(devtools::load_all(".", quiet = TRUE))
reps <- as.integer(Sys.getenv("REPS", "10"))
ancestor <- data.frame(optimum = 50, sigma = 1e4, r_max = 5, row.names = "sp1")
grid <- expand.grid(tau = c(0, 25, 100), rep = seq_len(reps))
run <- function(i) {
  g <- grid[i, ]
  set.seed(7000 + i)
  ls <- landscape_grid(6, disp_scale = 0.5)
  sim <- simulate_evolving_metacommunity(
    ls, ancestor,
    speciation_control("fission", rate = 0.01, partition = "random", completion = g$tau,
                       completion_dist = "exponential"),
    dispersal = 0.05, n_steps = 600, extirpation = 0.01, max_species = 1e5)
  out <- data.frame(g, tips = if (is.null(sim$tree)) 0 else ape::Ntip(sim$tree),
                    gamma = NA_real_, b = NA_real_, mu = NA_real_, lambda = NA_real_,
                    lineages = nrow(sim$species),
                    initiations = nrow(sim$events) / sum(pmin(sim$species$end, 600) -
                                                           sim$species$birth))
  if (out$tips < 5) return(out)
  out$gamma <- ape::gammaStat(sim$tree)
  brts <- sort(ape::branching.times(sim$tree), decreasing = TRUE)
  fit <- tryCatch(suppressMessages(utils::capture.output(
    f <- PBD::pbd_ML(brts, initparsopt = c(0.01, 0.005, 0.05), idparsopt = 1:3, exteq = 1,
                     soc = 2, cond = 1, btorph = 1, verbose = FALSE))) , error = function(e) NULL)
  # pbd_ML() returns -1 for every parameter when the fit fails
  if (!is.null(fit) && exists("f") && f$b >= 0) {
    out$b <- f$b
    out$mu <- f$mu_1
    out$lambda <- f$lambda_1
  }
  out
}
started <- Sys.time()
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), run,
                                         mc.cores = as.integer(Sys.getenv("CORES", "6"))))
summary_table <- stats::aggregate(cbind(tips, gamma, b, mu, lambda, initiations) ~ tau, res,
                                  function(x) round(stats::median(x, na.rm = TRUE), 4),
                                  na.action = na.pass)
summary_table$fits <- tapply(!is.na(res$lambda), res$tau, sum)
saveRDS(list(results = res, summary = summary_table, reps = reps,
             minutes = as.numeric(difftime(Sys.time(), started, units = "mins"))),
        "data-raw/evolving_pbd.rds")
print(summary_table)
