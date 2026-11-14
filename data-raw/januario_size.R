# Why PGLS under Brownian motion rejects too often (about 0.15 at alpha = 0.05) under the null in
# data-raw/januario_design.R when tips are sampled at random and the trait has measurement error.
#
#   Rscript data-raw/januario_size.R          (REPS=300 by default)
# Results in data-raw/januario_size.rds.
#
# Hypothesis: a random sample of 200 from 4000 tips contains pairs that are very close in the
# pruned tree. Independent contrasts divide by the square root of the branch length, so in those
# pairs both the iid measurement error of w and the non-Brownian part of log lambda_DR (computed on
# the complete tree, where the two species need not be sisters) become large. The variance of the
# contrast residuals then grows with cw^2, and the ordinary least-squares standard error is too
# small (heteroscedasticity tied to the regressor).
#
# Each replicate (b = 0, no background shifts) records, per sampling scheme and response:
#   y = log lambda_DR, or controls: y = Brownian on the complete tree, Brownian + iid noise;
# and per reliability R (1, 0.75) of w, p-values of:
#   pic:    PGLS (regression of contrasts through the origin);
#   hc3:    the same with a heteroscedasticity-consistent (HC3) standard error;
#   lambda: PGLS with Pagel's lambda estimated (phylolm);
#   es:     ES-sim (999 Brownian traits);
# and diagnostics: Spearman correlation of cy^2 with cw^2, and the ratio of the largest to the
# median contrast standard deviation (sqrt of the contrast variance).

source("data-raw/januario_design.R")

REPS <- as.integer(Sys.getenv("REPS", 300))

tests <- function(y, w, tr) {
  cy <- ape::pic(y, tr); cwv <- ape::pic(w, tr, var.contrasts = TRUE)
  cw <- cwv[, 1]
  m <- stats::lm(cy ~ cw - 1)
  se_hc3 <- sqrt(sandwich::vcovHC(m, type = "HC3")[1, 1])
  t_hc3 <- stats::coef(m)[[1]] / se_hc3
  d <- data.frame(y = y, w = w, row.names = tr$tip.label)
  lam <- tryCatch(summary(phylolm::phylolm(y ~ w, d, tr, model = "lambda"))$coefficients["w", "p.value"],
                  error = function(e) NA_real_)
  c(pic = summary(m)$coefficients["cw", "Pr(>|t|)"],
    hc3 = 2 * stats::pt(-abs(t_hc3), df = length(cw) - 1),
    lambda = lam, es = es_sim(y, w, tr)[["p"]],
    rho_sq = suppressWarnings(stats::cor(cy^2, cw^2, method = "spearman")),
    sd_ratio = max(1 / sqrt(cwv[, 2])) / stats::median(1 / sqrt(cwv[, 2])))
}

one_null <- function() {
  s <- sim_big_tree(b = 0)
  tree <- s$tree
  responses <- list(dr = log(dr_rates(tree)[tree$tip.label]),
                    bm = ape::rTraitCont(tree, sigma = 1)[tree$tip.label])
  responses$bm_noise <- responses$bm + stats::rnorm(length(responses$bm), 0, sd(responses$bm))
  sets <- list(clade = sample_clade(tree), random = sample_random(tree),
               dispersed = sample_dispersed(tree))
  out <- list()
  for (scheme in names(sets)) {
    tr <- ape::keep.tip(tree, sets[[scheme]])
    z <- as.vector(scale(s$x[tr$tip.label]))
    for (R in c(1, 0.75)) {
      w <- stats::setNames(sqrt(R) * z + sqrt(1 - R) * stats::rnorm(length(z)), tr$tip.label)
      for (resp in names(responses)) {
        out[[length(out) + 1]] <- data.frame(scheme = scheme, reliability = R, response = resp,
                                             t(tests(responses[[resp]][tr$tip.label], w, tr)))
      }
    }
  }
  do.call(rbind, out)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(2029)
start <- Sys.time()
out <- parallel::mclapply(seq_len(REPS), function(i) {
  tryCatch(cbind(rep = i, one_null()), error = function(e) data.frame(i = i, error = conditionMessage(e)))
}, mc.cores = CORES, mc.set.seed = TRUE)
ok <- vapply(out, function(o) is.null(o$error), logical(1))
res <- list(main = do.call(rbind, out[ok]), errors = do.call(rbind, out[!ok]),
            elapsed = difftime(Sys.time(), start, units = "mins"))
saveRDS(res, Sys.getenv("OUT", "data-raw/januario_size.rds"))
m <- res$main
print(aggregate(cbind(pic = pic < 0.05, hc3 = hc3 < 0.05, lambda = lambda < 0.05, es = es <= 0.05,
                      rho_sq, sd_ratio) ~ response + reliability + scheme, m, mean, na.action = na.pass),
      digits = 2)
cat("done in", format(res$elapsed), "with", sum(!ok), "errors\n")
