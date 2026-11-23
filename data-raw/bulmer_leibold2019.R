# Block H, step (i): does a genetic variance with mutation, drift and the transient Bulmer
# reduction of sexual reproduction give the local dominance of Leibold et al. (2019) that the
# engine did not reproduce (article "Replicating Leibold et al. (2019)", section "Three missing
# mechanisms")? Their default reproduction is sexual, where recombination restores most of the
# reduction of the variance by stabilising selection; the earlier runs took that reduction as
# permanent (`variance_selection = TRUE`, clonal), an upper bound on erosion.
# `variance_selection = "sexual"` keeps the reduction as linkage disequilibrium, which halves
# each generation (Bulmer 1971). Same design as data-raw/mechanisms_leibold2019_recruitment.R:
# three recruitment settings (Beverton-Holt, lottery with discrete generations, lottery with
# mortality 0.1 per step) crossed with mutation and drift from G = 200 and from no variance,
# Vm = 0.5, at Nm = 10; 10 runs of 10,000 generations each. The other genetic settings are in
# data-raw/mechanisms_leibold2019_recruitment.rds, whose runs are unchanged by this option.

source("data-raw/mechanisms_leibold2019_recruitment.R")

genetics_s <- data.frame(
  genetic = c("sexual_from_G", "sexual_from_zero"),
  G = c(200, 0),
  vm = 0.5,
  variance_selection = "sexual"
)

if (sys.nframe() == 0) {
  cores <- max(1, parallel::detectCores() - 2)
  design <- merge(merge(recruitments, genetics_s), data.frame(rep = 1:10))
  design <- design[order(design$mortality), ]
  t0 <- Sys.time()
  out <- parallel::mclapply(seq_len(nrow(design)), function(i) {
    d <- design[i, ]
    run <- run_leibold(G = d$G, seed = d$rep, n_steps = 10000, every = 100, nm = 10,
                       mutation_var = d$vm, recruitment = d$recruitment,
                       mortality = d$mortality, variance_selection = d$variance_selection)
    g_now <- mean(run$variance, na.rm = TRUE)
    summ <- summarise_leibold(run, var_total = 400 + g_now + 4)
    list(design = d, summary = summ, classes = classify_leibold(run, 5000),
         variance = run$variance, disequilibrium = run$disequilibrium, final = run$final,
         env = run$env)
  }, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(out, inherits, logical(1), "try-error")
  if (any(bad)) stop("runs failed: ", out[[which(bad)[1]]])
  message("Leibold, sexual Bulmer reduction: ", format(Sys.time() - t0))
  saveRDS(out, "data-raw/bulmer_leibold2019.rds")
}
