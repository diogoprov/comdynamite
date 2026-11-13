# Three mechanisms missing from the engine, tested on the local dominance of Leibold et al.
# (2019) that the engine did not reproduce (data-raw/replication_leibold2019.R and
# data-raw/mechanism_leibold2019.R, article "Replicating Leibold et al. (2019)"):
# * lottery recruitment into the microsites vacated by mortality (`recruitment = "lottery"`,
#   400 microsites per patch, their carrying capacity), in place of Beverton-Holt growth;
# * overlapping generations: a lottery with mortality 0.1 per step, since "one generation
#   equals 10 iterations of mortality, reproduction, and dispersal" in their Methods, against a
#   lottery with discrete generations; patch extinction stays 0.01 per generation;
# * a genetic variance built from new mutations: `mutation_var` starting from no variance, as
#   their loci start without standing variation, with and without the erosion of the variance
#   by selection (`variance_selection`, the Bulmer reduction taken as permanent, which holds
#   for asexual reproduction; their default is sexual, where recombination restores most of
#   it, so this arm is an upper bound on erosion).
# Each recruitment setting is crossed with four genetic settings at Nm = 10: a constant
# G = 200; mutation and drift from G = 200 with Vm = 0.5 (equilibrium 2 Ne Vm = 200 in a
# full patch); mutation and drift from no variance; and the same with erosion by selection;
# plus runs without evolution. 10 runs of 10,000 generations each.

source("data-raw/replication_leibold2019.R")

recruitments <- data.frame(recruitment = c("beverton_holt", "lottery", "lottery"),
                           mortality = c(1, 1, 0.1))
genetics_l <- data.frame(
  genetic = c("none", "constant", "mutation_from_G", "mutation_from_zero",
              "mutation_selection"),
  G = c(0, 200, 200, 0, 0),
  vm = c(NA, NA, 0.5, 0.5, 0.5),
  variance_selection = c(FALSE, FALSE, FALSE, FALSE, TRUE)
)

if (sys.nframe() == 0) {
  cores <- max(1, parallel::detectCores() - 2)
  design <- merge(merge(recruitments, genetics_l), data.frame(rep = 1:10))
  design <- design[order(design$mortality), ]
  t0 <- Sys.time()
  out <- parallel::mclapply(seq_len(nrow(design)), function(i) {
    d <- design[i, ]
    run <- run_leibold(G = d$G, seed = d$rep, n_steps = 10000, every = 100, nm = 10,
                       mutation_var = if (is.na(d$vm)) NULL else d$vm,
                       recruitment = d$recruitment, mortality = d$mortality,
                       variance_selection = d$variance_selection)
    # fitness differences use the variance of the populations where it is known
    g_now <- if (is.null(run$variance)) d$G else mean(run$variance, na.rm = TRUE)
    summ <- summarise_leibold(run, var_total = 400 + g_now + 4)
    list(design = d, summary = summ, classes = classify_leibold(run, 5000),
         variance = run$variance, final = run$final, env = run$env)
  }, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(out, inherits, logical(1), "try-error")
  if (any(bad)) stop("runs failed: ", out[[which(bad)[1]]])
  message("Leibold mechanisms: ", format(Sys.time() - t0))
  saveRDS(out, "data-raw/mechanisms_leibold2019_recruitment.rds")
}
