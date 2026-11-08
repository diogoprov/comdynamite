# history_rates(model = "sympatry") against the diversity-dependent model with biogeography of
# Drury et al. (2016) in RPANDA (GPL, only called here). RPANDA was archived on CRAN; install
# the archived TESS from CRAN and then remotes::install_github("hmorlon/PANDA").
# 1. The covariance matrix of the tips implied by the rates per branch against the one RPANDA
#    builds from the same stochastic map (RPANDA:::.VCV.rescale.DDexp_geog).
# 2. Recovery of sigma2 and r by RPANDA::fit_t_comp(model = "DDexp") from traits simulated
#    with trait_rates(edge_rates = history_rates(...)).
suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(RPANDA)
})
cores <- as.integer(Sys.getenv("CORES", "6"))
RNGkind("L'Ecuyer-CMRG")
set.seed(41)

# CreateGeoObject() fails on some maps without any change of region (an RPANDA bug: object
# 'old.Dist' not found; https://github.com/hmorlon/PANDA/issues/56), so such histories are
# drawn again
one_history <- function(n, q) {
  repeat {
    tree <- ape::rcoal(n)
    tree$tip.label <- paste0("t", seq_len(n))
    h <- trait_mk(tree, q = q, states = c("A", "B", "C"))
    if (any(lengths(h$maps) > 1)) return(h)
  }
}

vcv_diff <- function(i) {
  h <- one_history(15, 0.5)
  geo <- RPANDA::CreateGeoObject(h, h)
  scaled <- h
  class(scaled) <- "phylo"
  scaled$edge.length <- h$edge.length * history_rates(h, "sympatry", sigma2 = 1.3, r = -0.4)
  ours <- ape::vcv(scaled)
  theirs <- RPANDA:::.VCV.rescale.DDexp_geog(h, 1.3, -0.4, geo, check = FALSE)
  max(abs(ours - theirs[rownames(ours), colnames(ours)]))
}
diffs <- unlist(parallel::mclapply(seq_len(30), vcv_diff, mc.cores = cores))

recover <- function(i, r) {
  h <- one_history(60, 0.3)
  geo <- RPANDA::CreateGeoObject(h, h)
  x <- trait_rates(h, edge_rates = history_rates(h, "sympatry", sigma2 = 1, r = r))
  fit <- RPANDA::fit_t_comp(ape::as.phylo(h), x$tips, model = "DDexp", geography.object = geo)
  c(true_r = r, sig2 = fit$sig2, r = fit$r)
}
grid <- rep(c(-0.3, 0), each = 40)
fits <- do.call(rbind, parallel::mclapply(seq_along(grid), function(i) recover(i, grid[i]),
                                          mc.cores = cores, mc.set.seed = TRUE))
fits <- as.data.frame(fits)

cat("largest difference between covariance matrices (30 histories):", signif(max(diffs), 2), "\n")
cat("fit_t_comp(DDexp) on 60 tips, 3 regions, true sigma2 = 1:\n")
print(stats::aggregate(cbind(sig2, r) ~ true_r, fits, stats::median))
saveRDS(list(vcv_diff = diffs, fits = fits), "data-raw/validation_sympatry.rds")
