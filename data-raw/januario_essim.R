# ES-sim (Harvey & Rabosky 2018) on the design of data-raw/januario_power.R: does a test of the
# correlation between tip speciation rate and persistence, with a null from Brownian traits on
# the tree, detect the causal link that the PGLS of Januario et al. (2025) misses?
#
#   Rscript data-raw/januario_essim.R
# Results in data-raw/januario_essim.rds.
#
# ES-sim, implemented from the description in Harvey & Rabosky (2018): the statistic is the
# Pearson correlation between log tip rate and the trait; its null distribution comes from traits
# simulated by Brownian motion on the tree. The Pearson correlation does not change under affine
# transformations of the trait, so the root value and rate of the Brownian model do not matter
# and the null draws are L z, with L the Cholesky factor of the tree covariance. The tip rate is
# lambda_DR (Jetz et al. 2012), which the equal-splits measure of ES-sim approximates, as in
# Januario et al. (2025). Two-sided p-values from 999 null draws, for tau and for eps.

devtools::load_all(quiet = TRUE)
# the functions of januario_power.R, without its main run
for (e in parse("data-raw/januario_power.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
      identical(e[[3]][[1]], as.name("function"))) eval(e)
}

es_sim <- function(rate, trait, tree, n_null = 999) {
  root <- t(chol(ape::vcv(tree)[names(trait), names(trait)]))
  observed <- stats::cor(log(rate), trait)
  null <- vapply(seq_len(n_null), function(k) {
    stats::cor(log(rate), as.vector(root %*% stats::rnorm(length(trait))))
  }, numeric(1))
  c(r = observed, p = (1 + sum(abs(null) >= abs(observed))) / (n_null + 1))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(2027)
runs <- rbind(
  expand.grid(b = c(0, 0.25, 0.5, 1), window = 30, link = "dispersal", base = 0.01,
              rep = 1:200, stringsAsFactors = FALSE),
  expand.grid(b = c(0, 1), window = 30, link = "extirpation", base = 0.02, rep = 1:100,
              stringsAsFactors = FALSE))
start <- Sys.time()
out <- parallel::mclapply(seq_len(nrow(runs)), function(i) {
  r <- runs[i, ]
  tryCatch({
    x <- one_replicate(r$b, link = r$link, window = r$window, base = r$base)
    d <- x$data
    tau <- es_sim(stats::setNames(d$dr, d$species), stats::setNames(d$tau, d$species), x$tree)
    eps <- es_sim(stats::setNames(d$dr, d$species), stats::setNames(d$eps, d$species), x$tree)
    cbind(i = i, rep = r$rep, base = r$base, x$summary,
          r_es_tau = tau[["r"]], p_es_tau = tau[["p"]],
          r_es_eps = eps[["r"]], p_es_eps = eps[["p"]])
  }, error = function(e) data.frame(i = i, error = conditionMessage(e)))
}, mc.cores = 6, mc.set.seed = TRUE)
ok <- vapply(out, function(o) is.null(o$error), logical(1))
results <- list(main = do.call(rbind, out[ok]), errors = do.call(rbind, out[!ok]), runs = runs,
                elapsed = difftime(Sys.time(), start, units = "mins"))
saveRDS(results, "data-raw/januario_essim.rds")
cat("done in", format(results$elapsed), "with", sum(!ok), "errors\n")
m <- results$main
m$es_any <- pmin(1, 2 * pmin(m$p_es_tau, m$p_es_eps)) < 0.05 # Bonferroni over the two traits
print(aggregate(cbind(pgls = p_dr_any < 0.05, es_tau = p_es_tau < 0.05, es_eps = p_es_eps < 0.05,
                      es_any, r_es_tau) ~ link + b, m, mean))
