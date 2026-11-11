# Summaries for the article vignettes/articles/rates-assembly.Rmd, from the precomputed results
# of the validation, bridge and Januario experiments (no new simulation).
#
#   Rscript data-raw/vignette_rates_assembly.R
# Writes vignettes/articles/rates-assembly.rds.

# 1. identifiability of branch rates (data-raw/validation_rates.R, validation_bamm.R)
v <- readRDS("data-raw/validation_rates.rds")
rrphylo <- aggregate(cbind(with_change, with_rate, ceiling) ~ model, v$rrphylo, median)
bamm <- readRDS("data-raw/validation_bamm.rds")
bamm <- data.frame(replicates = nrow(bamm), spearman = stats::median(bamm$spearman),
                   prob_shift = mean(bamm$prob_shift), prob_none = mean(bamm$prob_none))

# 2. bridge experiments (data-raw/bridge_experiments.R)
b <- readRDS("data-raw/bridge_experiments.rds")$main
b$design <- paste(b$n_sites, "sites x", b$n_ind, "individuals")
# one metacommunity per i: the branch-rate correlations do not depend on h, c or the design
per_sim <- b[!duplicated(b$i), ]
turnover <- do.call(rbind, lapply(split(per_sim, per_sim[c("dispersal", "extirpation")]),
                                   function(d) data.frame(dispersal = d$dispersal[1],
                                                          extirpation = d$extirpation[1],
                                                          rho_turnover = mean(d$rho_turnover, na.rm = TRUE),
                                                          n = sum(!is.na(d$rho_turnover)),
                                                          n_sims = nrow(d))))
# metacommunities in which c could not be fitted (too few co-occurring populations; all in the
# cell with the lowest dispersal and the highest extirpation)
collapsed <- sum(tapply(b$c_free, b$i, function(x) all(is.na(x))))
h_rmse <- aggregate(cbind(error = (h_free - h)^2) ~ design, b, function(x) sqrt(mean(x)))
names(h_rmse)[2] <- "rmse_h"
# power for c = 1 against the 95% quantile of c-hat under c = 0 in the same cell
cell <- c("dispersal", "extirpation", "bandwidth", "h", "design")
power_c <- do.call(rbind, lapply(split(b, b[cell], drop = TRUE), function(d) {
  null <- d$c_free[d$c == 0]
  alt <- d$c_free[d$c == 1]
  if (sum(!is.na(null)) < 10 || sum(!is.na(alt)) < 10) return(NULL)
  q <- stats::quantile(null, 0.95, na.rm = TRUE)
  cbind(d[1, cell], power = mean(alt > q, na.rm = TRUE))
}))
power_c_design <- aggregate(power ~ design, power_c, function(x) c(mean = mean(x), min = min(x),
                                                                   max = max(x)))
power_c_design <- do.call(data.frame, power_c_design)
big <- b[b$n_sites == 100 & b$n_ind == 20, ]
tstats <- aggregate(cbind(ses_T_IP.IC, ses_T_PC.PR) ~ h + c, big, mean)

# 3. kernel width and identifiability of c (data-raw/bridge_narrow.R)
n <- readRDS("data-raw/bridge_narrow.rds")
narrow <- do.call(rbind, lapply(split(n, n$ratio), function(d) {
  qf <- stats::quantile(d$c_free[d$c == 0], 0.95, na.rm = TRUE)
  qx <- stats::quantile(d$c_fixed[d$c == 0], 0.95, na.rm = TRUE)
  do.call(rbind, lapply(c(1, 3), function(cc) {
    k <- d[d$c == cc, ]
    data.frame(ratio = d$ratio[1], c = cc, power_free = mean(k$c_free > qf, na.rm = TRUE),
               power_fixed = mean(k$c_fixed > qx, na.rm = TRUE))
  }))
}))

# 4. bootstrap of population_shifts() after PR 58 (data-raw/bridge_bootstrap.R)
bb <- readRDS("data-raw/bridge_bootstrap.rds")$boot
bb$design <- paste(bb$n_sites, "x", bb$n_ind)
bootstrap <- do.call(rbind, lapply(split(bb, bb[c("c", "design")], drop = TRUE), function(d) {
  d <- d[!is.na(d$p_displacement), ]
  data.frame(c = d$c[1], design = d$design[1], n = nrow(d),
             rejection = mean(d$p_displacement < 0.05),
             cover_c = mean(d$displacement_low <= d$c & d$displacement_high >= d$c),
             cover_h = mean(d$local_adaptation_low <= d$h & d$local_adaptation_high >= d$h),
             c_corrected = stats::median(d$displacement_corrected))
}))

# 5. power of the test of Januario et al. (2025) (data-raw/januario_power.R, januario_essim.R)
j <- readRDS("data-raw/januario_power.rds")$main
januario <- aggregate(cbind(pgls_dr = p_dr_any < 0.05, pgls_true = p_true_any < 0.05,
                            cor_dr_lambda, cor_tau_z, cor_eps_z) ~ link + b + window, j,
                      function(x) mean(x, na.rm = TRUE), na.action = na.pass)
es <- readRDS("data-raw/januario_essim.rds")$main
es$es_any <- pmin(1, 2 * pmin(es$p_es_tau, es$p_es_eps)) <= 0.05
essim <- aggregate(cbind(pgls_dr = p_dr_any < 0.05, pgls_true = p_true_any < 0.05,
                         es_tau = p_es_tau <= 0.05, es_eps = p_es_eps <= 0.05, es_any,
                         r_es_tau) ~ link + b, es, function(x) mean(x, na.rm = TRUE),
                   na.action = na.pass)
essim$n <- as.vector(table(es$link, es$b)[cbind(essim$link, as.character(essim$b))])

out <- list(rrphylo = rrphylo, bamm = bamm, turnover = turnover, collapsed = collapsed,
            n_metacommunities = nrow(per_sim), h_rmse = h_rmse, power_c = power_c,
            power_c_design = power_c_design, tstats = tstats, narrow = narrow,
            bootstrap = bootstrap, januario = januario, essim = essim)
saveRDS(out, "vignettes/articles/rates-assembly.rds")
str(out, max.level = 1)
