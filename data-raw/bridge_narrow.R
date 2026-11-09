devtools::load_all(quiet = TRUE)
RNGkind("L'Ecuyer-CMRG"); set.seed(99)
one <- function(r) {
  tree <- treeshapesim::sim_tree_shape(30)
  sim <- simulate_metacommunity(landscape_grid(10, env_bandwidth = 0.5), species_traits(tree),
    competition_matrix(tree, "phylo"), dispersal = 0.01, extirpation = 0.01, n_steps = 600)
  spread <- 9.8
  out <- list()
  for (ratio in c(0.2, 0.5, 1.8)) for (cc in c(0, 1, 3)) {
    w <- ratio * spread
    ind <- sample_individuals(sim, n = 20, local_adaptation = 0.5, displacement = cc,
                              displacement_width = w, population_sd = 1, intraspecific_sd = 5)
    s <- sample(unique(ind$site), 20); ind <- ind[ind$site %in% s & ind$individual <= 5, ]
    f <- tryCatch(suppressWarnings(population_shifts(ind)), error = function(e) NULL)
    g <- tryCatch(suppressWarnings(population_shifts(ind, width = w)), error = function(e) NULL)
    if (is.null(f) || is.null(g)) next
    out[[length(out)+1]] <- data.frame(r = r, ratio = ratio, c = cc, c_free = f$displacement,
                                       w_free = f$width / w, c_fixed = g$displacement)
  }
  do.call(rbind, out)
}
res <- do.call(rbind, parallel::mclapply(1:60, one, mc.cores = 2))
saveRDS(res, "data-raw/bridge_narrow.rds")
pw <- do.call(rbind, lapply(split(res, res$ratio), function(d) {
  n <- d[d$c == 0, ]; qf <- quantile(n$c_free, 0.95); qx <- quantile(n$c_fixed, 0.95)
  do.call(rbind, lapply(c(0, 1, 3), function(cc) { k <- d[d$c == cc, ]
    data.frame(ratio = d$ratio[1], c = cc, med_free = median(k$c_free), med_fixed = median(k$c_fixed),
               power_free = mean(k$c_free > qf), power_fixed = mean(k$c_fixed > qx),
               w_ratio = median(k$w_free)) }))
}))
print(format(pw, digits = 2), row.names = FALSE)
