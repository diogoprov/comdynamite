# The power analysis of Januario, Macedo-Rego & Rabosky (2025, Am. Nat. 205: 388-412) applied to
# trees where a trait truly drives speciation (the trees of data-raw/januario_power.R).
#
#   Rscript data-raw/januario_their_power.R
# Results in data-raw/januario_their_power.rds.
#
# Their procedure (Methods, "Power Analyses"): keep the observed tip rates and the tree fixed,
# simulate predictors under Brownian motion with a known correlation rho with the rates
# (rho ~ U(-1, 1)), count how often the PGLS slope has p < 0.05, and fit
#   P(signif) = (1 + exp(-a (1 - t))) / (1 + exp(-a (|rho| - t)))
# by maximum likelihood; t is the correlation detected with probability ~0.5.
# Their simulation code is not in the Dryad package (doi:10.5061/dryad.0p2ngf277), so two readings
# of "a Brownian trait with a known correlation" are run, with e Brownian on the tree:
#   raw:  x = rho * scale(y) + sqrt(1 - rho^2) * scale(e)   (correlation at the tips);
#   cond: x = rho * (y - root) / sigma_y + sqrt(1 - rho^2) * e   (bivariate Brownian model
#         conditional on y, sigma_y the Brownian rate of y: correlation among contrasts).
# On their tree and tip rates (Dryad files), "raw" gives thresholds close to the published ones
# (e.g. I_avg males: published t = 0.054 with DR and 0.140 with BAMM; raw 0.096 and 0.190; cond
# 0.239 for both), and only "raw" reproduces their lower thresholds with DR than with BAMM.
#
# Here the same procedure runs on trees simulated with lambda = exp(b X), X Brownian, and y =
# log lambda_DR. Alongside, the causal test: PGLS of y on the true X (a perfect measurement of
# the trait, which removes the persistence estimation step of januario_power.R).
# PGLS under Brownian motion is computed as regression through the origin of independent
# contrasts, which gives the same slope and test as nlme::gls with corBrownian.

devtools::load_all(quiet = TRUE)
src <- readLines("data-raw/januario_power.R")
eval(parse(text = src[grep("^sim_tree_trait", src):(grep("^# DR statistic", src) - 1)]))

pic_p <- function(y, x, tree) {
  cy <- ape::pic(y, tree); cx <- ape::pic(x, tree)
  m <- stats::lm(cy ~ cx - 1)
  c(p = summary(m)$coefficients["cx", "Pr(>|t|)"], r = sum(cx * cy) / sqrt(sum(cx^2) * sum(cy^2)))
}

fit_threshold <- function(rho, signif) {
  nll <- function(par) {
    a <- exp(par[1]); t <- par[2]
    p <- (1 + exp(-a * (1 - t))) / (1 + exp(-a * (abs(rho) - t)))
    p <- pmin(pmax(p, 1e-9), 1 - 1e-9)
    -sum(stats::dbinom(signif, 1, p, log = TRUE))
  }
  fits <- lapply(c(-0.5, 0.2, 0.8), function(t0) stats::optim(c(log(10), t0), nll))
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  c(a = exp(best$par[1]), t = best$par[2])
}

one_tree <- function(b, n_sim = 200) {
  s <- sim_tree_trait(n = 200, b = b)
  tree <- s$tree
  y <- log(dr_rates(tree))[tree$tip.label]
  z <- as.vector(scale(s$x)); names(z) <- tree$tip.label
  causal <- pic_p(y, z, tree)
  # Brownian rate and root value of y (GLS estimates)
  cy <- ape::pic(y, tree)
  sigma_y <- sqrt(mean(cy^2))
  root_y <- ape::ace(y, tree, method = "pic")$ace[1]
  rho <- stats::runif(n_sim, -1, 1)
  sim <- function(how) t(vapply(rho, function(r) {
    e <- ape::rTraitCont(tree, sigma = 1)[tree$tip.label]
    x <- if (how == "raw") r * as.vector(scale(y)) + sqrt(1 - r^2) * as.vector(scale(e)) else
      r * (y - root_y) / sigma_y + sqrt(1 - r^2) * e
    pic_p(y, x, tree)
  }, numeric(2)))
  raw <- sim("raw"); cond <- sim("cond")
  data.frame(b = b, p_causal = causal[["p"]], r_pic_causal = causal[["r"]],
             r_raw_causal = stats::cor(y, z),
             t_raw = fit_threshold(rho, raw[, "p"] < 0.05)[["t"]],
             t_cond = fit_threshold(rho, cond[, "p"] < 0.05)[["t"]],
             power_raw_rho_0.3 = mean(raw[abs(rho) >= 0.25 & abs(rho) <= 0.35, "p"] < 0.05),
             r_pic_vs_rho_raw = stats::cor(raw[, "r"], rho),
             r_pic_vs_rho_cond = stats::cor(cond[, "r"], rho))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(2027)
runs <- expand.grid(b = c(0, 0.25, 0.5, 1), rep = 1:200)
start <- Sys.time()
out <- parallel::mclapply(seq_len(nrow(runs)), function(i) {
  tryCatch(cbind(rep = runs$rep[i], one_tree(runs$b[i])),
           error = function(e) data.frame(i = i, error = conditionMessage(e)))
}, mc.cores = 6, mc.set.seed = TRUE)
ok <- vapply(out, function(o) is.null(o$error), logical(1))
res <- list(main = do.call(rbind, out[ok]), errors = do.call(rbind, out[!ok]),
            elapsed = difftime(Sys.time(), start, units = "mins"))
saveRDS(res, "data-raw/januario_their_power.rds")
m <- res$main
print(aggregate(cbind(r_raw_causal, r_pic_causal, t_raw, t_cond, power_raw_rho_0.3,
                      r_pic_vs_rho_raw, r_pic_vs_rho_cond) ~ b, m, median))
print(aggregate(cbind(causal_power = p_causal < 0.05) ~ b, m, mean))
cat("done in", format(res$elapsed), "with", sum(!ok), "errors\n")
