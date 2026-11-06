# Table 5 of Peres-Neto et al. (2012): variation partitioning with the logistic generator
# of Appendix C, on a balanced tree of 32 species and 31 sites. The values of the
# generator come from PGLMM_model_simulations.m of Ives & Helmus (2011).
#
# Neither the article nor the supplement says how the two vectors were obtained. The
# fractions of Table 5 agree with the vectors that generated the data, and its rejection
# rates for E and S agree with the global test of phylo_vector(), so both are reported.
# Run from the package root; writes data-raw/validation_table5.rds.
library(comdynamite)
library(parallel)

RNGkind("L'Ecuyer-CMRG")
set.seed(2012)
n_rep <- 1000
n_perm <- 999
cores <- max(1, detectCores() - 2)
tree <- ape::stree(32, "balanced")
tree$edge.length <- rep(1, nrow(tree$edge))
tree$tip.label <- paste0("sp", 1:32)

fractions <- c("abc", "ab", "bc", "a", "b", "c", "d")
one_replicate <- function(scenario) {
  g <- logistic_gradients(scenario)
  sim <- simulate_logistic(tree, g[c("x", "y")])
  keep <- rowSums(sim$abundance) > 0
  comm <- sim$abundance[keep, ]
  comm <- comm[, colSums(comm) > 0]
  x <- g$x_test[keep]
  y <- g$y_test[keep]
  known <- trait_partition(comm, sim$traits$slope[colnames(comm), ], x, y, n_perm = n_perm)
  f <- known$fractions
  rownames(f) <- f$fraction
  first <- phylo_vector(comm, tree, x, n_perm = n_perm)
  second <- phylo_vector(comm, tree, y, n_perm = n_perm)
  estimated <- cbind(first$pev_brown, second$pev_brown)
  tested <- c(a = FALSE, c = FALSE)
  if (!is.null(estimated)) {
    fit <- trait_partition(comm, estimated, x, y, n_perm = n_perm)$fractions
    rownames(fit) <- fit$fraction
    tested <- fit[c("a", "c"), "p_max"] <= 0.05
  }
  c(share = known$component_share, f[fractions, "r2"],
    reject_es = first$p_max <= 0.05 || second$p_max <= 0.05, reject_e = first$p_max <= 0.05,
    reject_s = second$p_max <= 0.05, reject_a = tested[[1]], reject_c = tested[[2]])
}

results <- t(vapply(1:7, function(scenario) {
  runs <- mclapply(seq_len(n_rep), function(i) one_replicate(scenario), mc.cores = cores)
  colMeans(do.call(rbind, runs))
}, numeric(13)))
colnames(results) <- c("share", fractions, "reject_es", "reject_e", "reject_s", "reject_a",
                       "reject_c")
rownames(results) <- as.character(as.roman(1:7))
saveRDS(results, "data-raw/validation_table5.rds", version = 2)
print(round(results, 3))
