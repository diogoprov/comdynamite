# Pilot for ABC: can local adaptation (h) and character displacement (c) be told apart from
# the traits of individuals? Three sets of summary statistics are compared:
#   ratios:  the three T-statistics of Violle et al. (2012), means over sites;
#   null:    ratios plus the mean standardised effect size of T_IP.IC under the null model that
#            exchanges population means within species (tstats(null_ip = "species"));
#   shifts:  null plus the estimates of h, c and log(w) from population_shifts().
# Criterion set before running, as in data-raw/pilot_identifiability.R: a parameter is
# recoverable with R2 >= 0.5 and coverage of the 90% interval between 85% and 95%, in
# pseudo-observations left out of the reference table. The criterion is a choice of the
# authors, not from the literature.
# Design: 200 metacommunities simulated once (30 species, 20 sites, no competition in the
# dynamics); each draw takes one at random and samples 10 individuals per population with
# intraspecific_sd = 2, h ~ U(0, 1) and c ~ U(0, 10). The kernel width is the default of
# sample_individuals(), the mean niche width of the species. Pseudo-observations come from 40
# metacommunities never used in the reference table (the metacommunity is unknown), and
# coverage above 95% (intervals too wide) is reported apart from coverage below 85%.
# History: a first run drew reference and test sets from the same 200 metacommunities with
# N_TEST = 200 (data-raw/pilot_tstats_identifiability_shared.rds). An independent audit pointed
# out that nearest neighbours then come largely from the metacommunity of the test draw, and
# that 200 pseudo-observations leave the 85-95% window only about 2.4 standard errors wide, so
# the design was changed to the one above before the final run. The first run gave R2 of 0.93
# (c) and 0.99 (h) with the shifts statistics, recoverable at tol = 0.05.
# Run from the root of the package, with the package abc installed. N_REF and N_TEST change the
# size, N_IND the individuals per population, CORES the number of workers and TOLS the
# tolerances. The pilot as planned is the default.
suppressMessages(devtools::load_all(".", quiet = TRUE))
n_ref <- as.integer(Sys.getenv("N_REF", "5000"))
n_test <- as.integer(Sys.getenv("N_TEST", "500"))
n_ind <- as.integer(Sys.getenv("N_IND", "10"))
cores <- as.integer(Sys.getenv("CORES", "6"))
output <- Sys.getenv("OUTPUT", "data-raw/pilot_tstats_identifiability.rds")
tolerances <- as.numeric(strsplit(Sys.getenv("TOLS", "0.01,0.02,0.05"), ",")[[1]])
n_meta <- 200
n_meta_test <- 40

started <- Sys.time()
metacommunities <- parallel::mclapply(seq_len(n_meta), function(m) {
  set.seed(m)
  tree <- treeshapesim::sim_tree_shape(30)
  simulate_metacommunity(landscape_implicit(20), species_traits(tree),
                         competition_matrix(tree, "none"), dispersal = 0.01, n_steps = 300,
                         seed_steps = 1:10)
}, mc.cores = cores)

set.seed(20261108)
n <- n_ref + n_test
parameters <- data.frame(local_adaptation = runif(n, 0, 1), displacement = runif(n, 0, 10),
                         metacommunity = c(sample.int(n_meta - n_meta_test, n_ref, TRUE),
                                           n_meta - n_meta_test +
                                             sample.int(n_meta_test, n_test, TRUE)))
limits <- rbind(local_adaptation = c(0, 1), displacement = c(0, 10))
names_parameters <- rownames(limits)

ratios <- c("T_IP.IC", "T_IC.IR", "T_PC.PR")
one_draw <- function(i) {
  set.seed(1000 + i)
  p <- parameters[i, ]
  ind <- sample_individuals(metacommunities[[p$metacommunity]], n = n_ind,
                            intraspecific_sd = 2, local_adaptation = p$local_adaptation,
                            displacement = p$displacement)
  t_site <- tstats(ind, n_perm = 49, null_ip = "species")
  shifts <- population_shifts(ind)
  c(colMeans(t_site[ratios], na.rm = TRUE),
    ses_species = mean(t_site$ses_T_IP.IC, na.rm = TRUE),
    h_hat = shifts$local_adaptation, c_hat = shifts$displacement,
    log_w_hat = log(shifts$width))
}
statistics <- do.call(rbind, parallel::mclapply(seq_len(n), one_draw, mc.cores = cores))
table <- cbind(parameters, statistics)
minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))

views <- list(ratios = ratios, null = c(ratios, "ses_species"),
              shifts = c(ratios, "ses_species", "h_hat", "c_hat", "log_w_hat"))
is_test <- seq_len(n) > n_ref
weighted_quantile <- function(x, w, p) {
  o <- order(x)
  x[o][pmin(findInterval(p, cumsum(w[o]) / sum(w)) + 1, length(x))]
}
evaluate <- function(view, tol) {
  s <- as.matrix(table[, views[[view]]])
  complete <- rowSums(!is.finite(s)) == 0
  reference <- !is_test & complete
  test <- which(is_test & complete)
  rows <- lapply(test, function(i) {
    fit <- suppressWarnings(abc::abc(target = s[i, ], param = table[reference, names_parameters],
                                     sumstat = s[reference, , drop = FALSE], tol = tol,
                                     method = "loclinear", transf = rep("logit", 2),
                                     logit.bounds = limits))
    q <- apply(fit$adj.values, 2, weighted_quantile, w = fit$weights, p = c(0.05, 0.5, 0.95))
    data.frame(parameter = names_parameters, truth = unlist(table[i, names_parameters]),
               low = q[1, ], estimate = q[2, ], high = q[3, ], row = i)
  })
  out <- do.call(rbind, rows)
  out$view <- view
  out$tol <- tol
  attr(out, "used") <- c(reference = sum(reference), test = length(test))
  out
}
runs <- list()
for (view in names(views)) for (tol in tolerances) runs[[paste(view, tol)]] <- evaluate(view, tol)
estimates <- do.call(rbind, runs)

r2 <- function(d) 1 - sum((d$truth - d$estimate)^2) / sum((d$truth - mean(d$truth))^2)
summarise <- function(d) {
  range <- diff(limits[d$parameter[1], ])
  boot <- replicate(1000, r2(d[sample.int(nrow(d), replace = TRUE), ]))
  covered <- d$truth >= d$low & d$truth <= d$high
  data.frame(view = d$view[1], tol = d$tol[1], parameter = d$parameter[1], n = nrow(d),
             r2 = r2(d), r2_low = quantile(boot, 0.025), r2_high = quantile(boot, 0.975),
             bias = mean(d$estimate - d$truth) / range, coverage = mean(covered),
             coverage_se = sqrt(mean(covered) * (1 - mean(covered)) / nrow(d)),
             width = mean(d$high - d$low) / range,
             recoverable = r2(d) >= 0.5 & mean(covered) >= 0.85 & mean(covered) <= 0.95,
             too_wide = mean(covered) > 0.95)
}
set.seed(1)
results <- do.call(rbind, lapply(split(estimates,
                                       list(estimates$view, estimates$tol, estimates$parameter)),
                                 summarise))
rownames(results) <- NULL
results <- results[order(results$tol, results$parameter, match(results$view, names(views))), ]

# the raw estimates of population_shifts() against the truth, in the reference table
raw <- with(table[!is_test, ], c(cor_h = cor(local_adaptation, h_hat, use = "complete.obs"),
                                 cor_c = cor(displacement, c_hat, use = "complete.obs")))

saveRDS(list(results = results, estimates = estimates, table = table, raw = raw,
             minutes = minutes, used = lapply(runs, attr, "used"),
             settings = list(n_ref = n_ref, n_test = n_test, n_ind = n_ind)), output)
cat("minutes:", round(minutes, 1), "\n")
print(round(raw, 2))
print(results, digits = 2, row.names = FALSE)
