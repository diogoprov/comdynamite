# The internal simulator of trees with heritable shifts in speciation (R/sim_shift_tree.R),
# used in the simulations of the DAMOCLES paper.
# 1. Rates by regime: one shift with ratio 1, 2 or 4, trees of 112 tips, 200 trees each. In the
#    reconstructed tree, speciation events per unit of branch length in each regime estimate a
#    pure-birth rate; its ratio between the shifted regime and the background should be close
#    to the true ratio, and 1 when the ratio is 1. With extinction (epsilon = 0.5) the
#    reconstructed rate is not lambda, but it scales with lambda, so the ratio holds.
# 2. Conditioning: share of the tree length in terminal edges with the general sampling
#    approach (gsa_factor = 2) against cutting within the first interval with n lineages only
#    (gsa_factor just above 1, close to stopping when the n-th lineage appears), 200 trees of
#    112 tips without shifts each.
suppressMessages(devtools::load_all(".", quiet = TRUE))
RNGkind("L'Ecuyer-CMRG")
set.seed(2010)
cores <- as.integer(Sys.getenv("CORES", "6"))

ratio_estimate <- function(sim) {
  tree <- sim$tree
  internal <- tree$edge[, 2] > ape::Ntip(tree)
  births <- tapply(internal, sim$regimes, sum)
  length_in <- tapply(tree$edge.length, sim$regimes, sum)
  (births[["1"]] / length_in[["1"]]) / (births[["0"]] / length_in[["0"]])
}
recovery <- do.call(rbind, lapply(c(1, 2, 4), function(r) {
  est <- unlist(parallel::mclapply(1:200, function(i) {
    ratio_estimate(sim_shift_tree(112, n_shifts = 1, ratio = r))
  }, mc.cores = cores))
  data.frame(ratio = r, median = stats::median(est), mean_log = exp(mean(log(est))),
             q05 = stats::quantile(est, 0.05), q95 = stats::quantile(est, 0.95))
}))
cat("1. estimated ratio of speciation, shifted regime against background (200 trees each)\n")
print(recovery, row.names = FALSE, digits = 3)

terminal_share <- function(tree) {
  terminal <- tree$edge[, 2] <= ape::Ntip(tree)
  sum(tree$edge.length[terminal]) / sum(tree$edge.length)
}
conditioning <- do.call(rbind, lapply(c(1.0001, 2), function(g) {
  share <- unlist(parallel::mclapply(1:200, function(i) {
    terminal_share(sim_shift_tree(112, n_shifts = 0, gsa_factor = g)$tree)
  }, mc.cores = cores))
  data.frame(gsa_factor = g, terminal_share = mean(share), sd = stats::sd(share))
}))
cat("2. share of the tree length in terminal edges, no shifts, 112 tips\n")
print(conditioning, row.names = FALSE, digits = 3)
saveRDS(list(recovery = recovery, conditioning = conditioning),
        "data-raw/validation_shift_tree.rds")
