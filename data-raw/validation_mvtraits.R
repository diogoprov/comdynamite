# Multivariate traits from trait_ou() with rate matrices from rate_matrix(), against mvMORPH
# (GPL-2 | GPL-3, only called here; Clavel et al. 2015, doi:10.1111/2041-210X.12420).
# 1. Covariance of the traits at the tips over many simulations, under Brownian motion with one
#    rate matrix per regime and regimes changing within branches (a map from trait_mk()), against
#    the expected sum_k kronecker(R_k, C_k), with C_k from mvMORPH:::vcvSplit() (internal), and against the
#    covariance of the same number of simulations by mvMORPH::mvSIM(model = "BMM").
# 2. The same for an Ornstein-Uhlenbeck process with one rate matrix, against the covariance of
#    Hansen (1997) with the root fixed: R * (exp(-a (d_i + d_j - 2 s_ij)) - exp(-a (d_i + d_j))) / (2a).
# 3. Recovery of the rate matrices by mvMORPH::mvBM() with each constraint, from traits simulated
#    with matrices of that structure.
suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(mvMORPH)
})
cores <- as.integer(Sys.getenv("CORES", "6"))
RNGkind("L'Ecuyer-CMRG")
set.seed(7)

# 1. covariance under BM with two regimes ---------------------------------------------------
tree <- treeshapesim::sim_tree_shape(12)
repeat {
  hist <- trait_mk(tree, q = 0.6, states = c("a", "b"))
  if (sum(lengths(hist$maps) > 1) >= 3 && all(colSums(hist$mapped.edge) > 0)) break
}
rates_bmm <- rate_matrix(c(1, 2), correlation = list(a = 0.6, b = -0.4), constraint = "variance")
split <- mvMORPH:::vcvSplit(hist)  # internal; one matrix per regime, in the order of colnames(hist$mapped.edge)
tips <- rownames(split[[1]])
expected <- Reduce(`+`, Map(kronecker, rates_bmm[colnames(hist$mapped.edge)], split))
n_sim <- 20000
ours <- do.call(rbind, parallel::mclapply(seq_len(n_sim), function(i) {
  as.vector(trait_ou(hist, sigma2 = rates_bmm, regimes = hist)[tips, ])
}, mc.cores = cores))
theirs <- mvSIM(hist, nsim = n_sim, model = "BMM",
                param = list(sigma = unname(rates_bmm[colnames(hist$mapped.edge)]),
                             theta = c(0, 0), ntraits = 2))
theirs <- do.call(rbind, lapply(theirs, function(x) as.vector(x[tips, ])))
# Monte Carlo standard error of each covariance: sqrt((s_ii s_jj + s_ij^2) / n)
mc_se <- sqrt((outer(diag(expected), diag(expected)) + expected^2) / n_sim)
cov_check <- data.frame(
  source = c("trait_ou", "mvSIM"),
  max_abs_diff = c(max(abs(stats::cov(ours) - expected)), max(abs(stats::cov(theirs) - expected))),
  max_z = c(max(abs(stats::cov(ours) - expected) / mc_se),
            max(abs(stats::cov(theirs) - expected) / mc_se)),
  mean_abs_z = c(mean(abs(stats::cov(ours) - expected) / mc_se),
                 mean(abs(stats::cov(theirs) - expected) / mc_se)))

# 2. covariance under OU with one rate matrix ------------------------------------------------
alpha <- 0.8
rates_ou <- rate_matrix(c(1, 0.5), correlation = 0.7)
depth <- ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)]
shared <- ape::vcv(tree)[tree$tip.label, tree$tip.label]
dsum <- outer(depth, depth, `+`)
c_ou <- (exp(-alpha * (dsum - 2 * shared)) - exp(-alpha * dsum)) / (2 * alpha)
expected_ou <- kronecker(rates_ou, c_ou)
ou <- do.call(rbind, parallel::mclapply(seq_len(n_sim), function(i) {
  as.vector(trait_ou(tree, alpha = alpha, sigma2 = rates_ou)[tree$tip.label, ])
}, mc.cores = cores))
mc_se_ou <- sqrt((outer(diag(expected_ou), diag(expected_ou)) + expected_ou^2) / n_sim)
ou_check <- c(max_abs_diff = max(abs(stats::cov(ou) - expected_ou)),
              max_z = max(abs(stats::cov(ou) - expected_ou) / mc_se_ou),
              mean_abs_z = mean(abs(stats::cov(ou) - expected_ou) / mc_se_ou))

# 3. recovery by mvBM() ----------------------------------------------------------------------
big <- treeshapesim::sim_tree_shape(100)
big$edge.length <- big$edge.length / max(ape::node.depth.edgelength(big))  # height 1
repeat {
  big_hist <- trait_mk(big, q = 1, states = c("a", "b"))
  share <- colSums(big_hist$mapped.edge) / sum(big$edge.length)
  if (all(share > 0.3)) break
}
truths <- list(
  default = rate_matrix(c(1, 2), correlation = 0.6),
  diagonal = rate_matrix(c(1, 2), constraint = "diagonal"),
  equaldiagonal = rate_matrix(1.5, constraint = "equaldiagonal", n_traits = 2),
  equal = rate_matrix(1.5, correlation = 0.6, constraint = "equal", n_traits = 2),
  proportional = rate_matrix(c(1, 2), correlation = 0.6, constraint = "proportional",
                             scale = c(a = 1, b = 4)),
  shared = rate_matrix(c(1, 2), correlation = 0.6, constraint = "shared",
                       eigenvalues = list(a = c(3, 0.5), b = c(1, 2))),
  variance = rate_matrix(c(1, 2), correlation = list(a = 0.6, b = -0.6),
                         constraint = "variance"),
  correlation = rate_matrix(variance = list(a = c(1, 2), b = c(4, 0.5)), correlation = 0.6,
                            constraint = "correlation"))
n_rep <- 30
fit_one <- function(constraint, i) {
  truth <- truths[[constraint]]
  multi <- is.list(truth)
  x <- trait_ou(if (multi) big_hist else big, sigma2 = truth,
                 regimes = if (multi) big_hist)
  constraint_mv <- if (constraint == "default") FALSE else constraint
  fit <- if (multi) {
    mvBM(big_hist, x, model = "BMM", param = list(constraint = constraint_mv),
         echo = FALSE, diagnostic = FALSE)
  } else {
    mvBM(big, x, model = "BM1", param = list(constraint = constraint_mv),
         echo = FALSE, diagnostic = FALSE)
  }
  est <- fit$sigma
  if (multi) {
    est <- lapply(seq_len(dim(est)[3]), function(k) est[, , k])
    names(est) <- colnames(big_hist$mapped.edge)
  }
  data.frame(constraint = constraint, rep = i,
             regime = if (multi) rep(names(est), each = 3) else "all",
             element = c("var1", "cov12", "var2"),
             estimate = unlist(lapply(if (multi) est else list(est),
                                      function(m) c(m[1, 1], m[1, 2], m[2, 2]))),
             convergence = fit$convergence)
}
grid <- expand.grid(constraint = names(truths), rep = seq_len(n_rep), stringsAsFactors = FALSE)
fits <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), function(j) {
  fit_one(grid$constraint[j], grid$rep[j])
}, mc.cores = cores))
true_values <- do.call(rbind, lapply(names(truths), function(cn) {
  m <- truths[[cn]]
  if (!is.list(m)) m <- list(all = m)
  do.call(rbind, lapply(names(m), function(r) {
    data.frame(constraint = cn, regime = r, element = c("var1", "cov12", "var2"),
               truth = c(m[[r]][1, 1], m[[r]][1, 2], m[[r]][2, 2]))
  }))
}))
fits$converged <- fits$convergence == 0
recovery <- merge(
  stats::aggregate(cbind(median_estimate = estimate) ~ constraint + regime + element, fits,
                   stats::median),
  stats::aggregate(converged ~ constraint + regime + element, fits, mean))
recovery <- merge(true_values, recovery)
recovery <- recovery[order(match(recovery$constraint, names(truths)), recovery$regime,
                           match(recovery$element, c("var1", "cov12", "var2"))), ]

print(cov_check)
print(ou_check)
print(recovery, row.names = FALSE, digits = 3)
saveRDS(list(cov_check = cov_check, ou_check = ou_check, recovery = recovery, fits = fits,
             mvMORPH = as.character(utils::packageVersion("mvMORPH"))),
        "data-raw/validation_mvtraits.rds")

# Results (mvMORPH 1.2.3, seed 7, 20000 simulations, 30 fits per constraint on 100 tips):
# 1. BM, two regimes changing within branches (12 tips, 300 covariances of 24 tip-trait values):
#    largest difference from the expected covariance 0.037 for trait_ou() and 0.035 for mvSIM();
#    in Monte Carlo standard errors, at most 2.86 (mean 0.79) for trait_ou() and 2.76 (mean
#    0.92) for mvSIM(). Both are within simulation error of sum_k kronecker(R_k, C_k).
# 2. OU, alpha = 0.8: largest difference 0.011, at most 3.06 standard errors (mean 0.91).
# 3. Median of 30 estimates by mvBM() with the constraint of the simulation (truth in brackets):
#    default        var1 1.02 (1)  cov 0.77 (0.85)  var2 1.92 (2)
#    diagonal       var1 1.00 (1)  cov 0 (0)        var2 2.03 (2)
#    equaldiagonal  var 1.51 (1.5)
#    equal          var 1.56 (1.5)  cov 0.92 (0.90)
#    proportional   a: 0.99, 0.81, 1.94 (1, 0.85, 2); b: 3.85, 3.13, 7.73 (4, 3.39, 8)
#    shared         a: 1.02, 1.03, 2.28 (1.12, 1.08, 2.39); b: 1.71, -0.44, 1.16 (1.75, -0.43, 1.25)
#    variance       var 1.00, 1.92 (1, 2); cov a 0.83 (0.85), b -0.87 (-0.85)
#    correlation    a: 1.00, 0.89, 1.97 (1, 0.85, 2); b: 3.98, 0.91, 0.51 (4, 0.85, 0.5)
#    mvBM() reported convergence (code 0) in 90% to 93% of the BM1 fits and in all BMM fits.
