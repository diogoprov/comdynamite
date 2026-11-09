# Diagnostic of januario_power.R: why the PGLS of lambda_DR has no power. Same replicates
# (dispersal link, window 30), tested with PGLS, OLS and Spearman correlations.
src <- readLines("data-raw/januario_power.R")
eval(parse(text = src[1:(grep("^set.seed", src)[1] - 1)]))
set.seed(7)
res <- parallel::mclapply(rep(c(0, 1), each = 60), function(b) {
  r <- one_replicate(b, link = "dispersal", window = 30, base = 0.01)
  d <- r$data
  c(b = b, p_pgls = r$summary$p_dr_any,
    p_ols = anova(lm(log(dr) ~ eps + tau, d), lm(log(dr) ~ 1, d))$`Pr(>F)`[2],
    p_spear_tau = cor.test(d$dr, d$tau, method = "spearman", exact = FALSE)$p.value,
    p_spear_z = cor.test(d$dr, d$z, method = "spearman", exact = FALSE)$p.value)
}, mc.cores = 6)
res <- as.data.frame(do.call(rbind, res))
print(aggregate(cbind(p_pgls, p_ols, p_spear_tau, p_spear_z) ~ b, res, function(p) mean(p < 0.05)))
