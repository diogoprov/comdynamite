# Block H, step (ii): does a genetic architecture of a few loci of large effect give the local
# dominance of Leibold et al. (2019) that the engine did not reproduce with the infinitesimal
# model (article "Replicating Leibold et al. (2019)", sections "Three missing mechanisms" and
# "Erosion under sexual reproduction")? Their model has 10 diploid loci with two alleles and
# additive effects, the genotype is the mean of the 20 alleles (so one allele is worth 0.05 of
# the range from 0 to 1), mutation changes an allele with probability 1e-4, and there is no
# standing variation at the start: the species of patch type 0.4 has 8 of its 20 alleles at 1
# and that of type 0.6 has 12. With the trait multiplied by 100, as in
# data-raw/replication_leibold2019.R, one allele is worth 5 units: `loci = list(n = 10,
# effect = 5, mutation = 1e-4)` in `evolution_control()`, where the genetic variance of each
# population follows its allele frequencies. Same design as
# data-raw/mechanisms_leibold2019_recruitment.R: three recruitment settings (Beverton-Holt,
# lottery with discrete generations, lottery with mortality 0.1 per step) at Nm = 10, crossed
# with mutation 1e-4 (theirs), 1e-3 and 1e-2 for the 10 loci, plus a single locus of effect 10
# from 40 to 60 with mutation 1e-3; 10 runs of 10,000 generations each.
#
# Rscript data-raw/loci_leibold2019.R [generations] [runs] [recruitment settings, e.g. 1,2]

source("data-raw/mechanisms_leibold2019_recruitment.R")

genetics_loci <- data.frame(
  genetic = c("loci10_mu1e-4", "loci10_mu1e-3", "loci10_mu1e-2", "locus1_mu1e-3"),
  n = c(10, 10, 10, 1),
  effect = c(5, 5, 5, 10),
  origin = c(0, 0, 0, 40),
  mutation = c(1e-4, 1e-3, 1e-2, 1e-3)
)

recruitment_label <- function(d) {
  c("Beverton-Holt", "lottery, discrete", "lottery, overlapping")[
    match(paste(d$recruitment, d$mortality), c("beverton_holt 1", "lottery 1", "lottery 0.1"))]
}

# one row per recruitment and genetic setting: runs that kept both species at generation
# `at`, mean genetic variance at the end, mean absolute fitness difference within patch types
# at `at` and at generations 4000-5000, share of the variation in composition explained by
# patch type at generations 4000-5000, and classes of patches at `at` (the last four in the
# runs that kept both species), as in data-raw/figures_mechanisms.R
loci_table <- function(out, at = 5000) {
  rows <- lapply(out, function(r) {
    s <- r$summary
    both <- isTRUE(s$both[s$step == at])
    late <- s[s$step > at - 1000 & s$step <= at & s$both, ]
    fit <- function(x) mean((abs(x$fitdiff_e40) + abs(x$fitdiff_e60)) / 2, na.rm = TRUE)
    cls <- if (both) as.numeric(r$classes) else rep(NA_real_, 3)
    data.frame(recruitment = recruitment_label(r$design), genetic = r$design$genetic,
               both = both, G_end = mean(r$variance, na.rm = TRUE),
               fitdiff_at = if (both) fit(s[s$step == at, ]) else NA_real_,
               fitdiff_late = if (both) fit(late) else NA_real_,
               env_comp = if (both) mean(late$env_comp, na.rm = TRUE) else NA_real_,
               preadapted = cls[1], switched = cls[2], shared = cls[3])
  })
  d <- do.call(rbind, rows)
  d$recruitment <- factor(d$recruitment, unique(d$recruitment))
  d$genetic <- factor(d$genetic, unique(genetics_loci$genetic))
  n <- aggregate(both ~ recruitment + genetic, d, sum)
  means <- aggregate(cbind(G_end, fitdiff_at, fitdiff_late, env_comp, preadapted, switched,
                           shared) ~ recruitment + genetic, d,
                     function(v) mean(v, na.rm = TRUE), na.action = na.pass)
  merge(n, means)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  n_steps <- if (length(args)) as.integer(args[1]) else 10000
  n_rep <- if (length(args) > 1) as.integer(args[2]) else 10
  parts <- if (length(args) > 2) as.integer(strsplit(args[3], ",")[[1]]) else 1:3
  cores <- max(1, parallel::detectCores() - 2)
  design <- merge(merge(recruitments[parts, ], genetics_loci), data.frame(rep = seq_len(n_rep)))
  design <- design[order(design$mortality), ]
  at <- min(5000, n_steps)
  t0 <- Sys.time()
  out <- parallel::mclapply(seq_len(nrow(design)), function(i) {
    d <- design[i, ]
    t1 <- Sys.time()
    run <- run_leibold(G = 0, seed = d$rep, n_steps = n_steps, every = 100, nm = 10,
                       recruitment = d$recruitment, mortality = d$mortality,
                       loci = list(n = d$n, effect = d$effect, mutation = d$mutation,
                                   origin = d$origin))
    # fitness differences use the variance of the populations at the end
    g_now <- mean(run$variance, na.rm = TRUE)
    summ <- summarise_leibold(run, var_total = 400 + (if (is.finite(g_now)) g_now else 0) + 4)
    list(design = d, summary = summ, classes = classify_leibold(run, at),
         variance = run$variance, allele_frequency = run$allele_frequency, final = run$final,
         env = run$env, preadapted = run$preadapted,
         seconds = as.numeric(Sys.time() - t1, units = "secs"))
  }, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(out, inherits, logical(1), "try-error")
  if (any(bad)) stop("runs failed: ", out[[which(bad)[1]]])
  message("Leibold, loci of large effect: ", format(Sys.time() - t0))
  full <- n_steps == 10000 && n_rep == 10
  out_file <- if (full) "data-raw/loci_leibold2019.rds" else {
    file.path(tempdir(), "loci_leibold2019_pilot.rds")
  }
  saveRDS(out, out_file)
  message("saved ", out_file)
  print(loci_table(out, at), digits = 3, row.names = FALSE)
}
