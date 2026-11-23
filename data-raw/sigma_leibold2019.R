# Stronger selection per generation in the replication of Leibold et al. (2019).
# In their individual-based model (SI Appendix, IBM description and Table S2), fecundity does
# not depend on the phenotype (offspring per patch b sum(1 - N/K) P_D): selection acts only on
# survival, S = (1 - d) f per time step with d = 0.1 and f = exp(-(z - e)^2 / (2 omega^2)),
# omega = 0.2 (20 in the x100 units of the article), and a generation time of 1 / d = 10 steps
# for fully adapted individuals.
# Lifetime reproduction is then proportional to the expected lifespan, 1 / (1 - (1 - d) w),
# so relative fitness per generation is W(x) = d / (1 - (1 - d) w(x)), which is not Gaussian.
# Two Gaussian niche widths per generation bracket it: matching its curvature at the optimum,
# sigma = omega sqrt(d / (1 - d)) = 6.67, and matching its value at the other patch type
# (x = 20, W = 0.22), sigma = 20 / sqrt(2 log(1 / 0.22)) = 11.5. The article used sigma = 20.
# Design: the main replication (constant G of 0 to 400, Nm of 1, 10 and 100, Beverton-Holt)
# and their genetics (10 loci of effect 5, mutation 1e-4 per allele, no standing variation)
# with Beverton-Holt and with a lottery with mortality 0.1 (overlapping generations), for
# each sigma; 10 runs of 10,000 generations each. Seeds as in the other scripts.

source("data-raw/mechanisms_leibold2019_recruitment.R")

w_gen <- function(x, omega = 20, d = 0.1) d / (1 - (1 - d) * exp(-x^2 / (2 * omega^2)))
sigmas <- c(curvature = 20 * sqrt(0.1 / 0.9),
            other_patch = 20 / sqrt(2 * log(1 / w_gen(20))))

if (sys.nframe() == 0) {
  cores <- max(1, parallel::detectCores() - 2)
  constant <- expand.grid(G = c(0, 20, 50, 100, 200, 400), nm = c(1, 10, 100),
                          recruitment = "beverton_holt", mortality = 1, loci = FALSE,
                          rep = 1:10, stringsAsFactors = FALSE)
  loci_runs <- merge(data.frame(G = 0, nm = c(1, 10, 100), loci = TRUE),
                     merge(data.frame(recruitment = c("beverton_holt", "lottery"),
                                      mortality = c(1, 0.1)),
                           data.frame(rep = 1:10)))
  design <- merge(rbind(constant, loci_runs[names(constant)]),
                  data.frame(sigma_name = names(sigmas), sigma = unname(sigmas)))
  design <- design[order(design$mortality), ]
  t0 <- Sys.time()
  out <- parallel::mclapply(seq_len(nrow(design)), function(i) {
    d <- design[i, ]
    loci <- if (d$loci) list(n = 10, effect = 5, mutation = 1e-4)
    seed <- if (d$loci) d$rep else 1000 + d$rep
    run <- run_leibold(G = d$G, seed = seed, n_steps = 10000, every = 100, nm = d$nm,
                       recruitment = d$recruitment, mortality = d$mortality, loci = loci,
                       sigma = d$sigma)
    g_now <- if (is.null(run$variance)) d$G else mean(run$variance, na.rm = TRUE)
    if (!is.finite(g_now)) g_now <- 0
    summ <- summarise_leibold(run, var_total = d$sigma^2 + g_now + 4)
    list(design = d, summary = summ, classes = classify_leibold(run, 5000),
         variance = run$variance, final = run$final, env = run$env)
  }, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(out, inherits, logical(1), "try-error")
  if (any(bad)) stop("runs failed: ", out[[which(bad)[1]]])
  message("Leibold, stronger selection: ", format(Sys.time() - t0))
  saveRDS(out, "data-raw/sigma_leibold2019.rds")
}
