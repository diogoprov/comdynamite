# Spatial diagnostic for the DAMOCLES paper (dev/plano-artigo-damocles.md, S5): do the 14
# communities of Benicio et al. (2021) behave as independent realisations of DAMOCLES on the
# same tree, as the shared fit assumes?
# 1. The 112 species on the best tree of data-raw/benicio2021_trees.R (community_112), rescaled to
#    height one; one pair of rates for all communities (damocles_fit(shared = TRUE)).
# 2. Observed: Jaccard similarity between every pair of communities, geographic distance
#    (haversine, km), the regions (inland: latitude north of 23 S, seven areas; coast and Serra
#    do Mar: seven areas). Statistics: mean similarity, mean within-region minus between-region
#    similarity, and the Pearson correlation between similarity and distance (Mantel r), over all
#    pairs and over the pairs within a region (distance decay apart from the regional split).
# 3. Null: N_SIM sets of 14 communities simulated independently with the fitted rates on the same
#    tree (simulate_damocles()), with the same statistics. Under independence the within-minus-
#    between difference and Mantel r are zero on average; an excess in the data means shared
#    history or shared environment (the two are not separable here).
# 4. The likelihood-ratio test of separate rates for the two regions (2 df), which treats the
#    communities as independent.
# Data: DATA_DIR/benicio2021/ (Dryad) and DATA_DIR/benicio2021_trees.rds (not committed). Output:
# data-raw/s5_benicio_spatial.rds.
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE))
data_dir <- Sys.getenv("DATA_DIR", "data-raw")
n_sim <- as.integer(Sys.getenv("N_SIM", "999"))
RNGkind("L'Ecuyer-CMRG")
set.seed(505)

tree <- readRDS(file.path(data_dir, "benicio2021_trees.rds"))$community_112$best
height <- max(ape::node.depth.edgelength(tree))
tree$edge.length <- tree$edge.length / height
norm <- function(x) sub("_+$", "", gsub("'", "", x))
comm <- read.csv(file.path(data_dir, "benicio2021", "frog_data.csv"), row.names = 1,
                 check.names = FALSE)
taxonomy <- read.csv("data-raw/benicio2021_taxonomy.csv")
colnames(comm) <- taxonomy$name[match(norm(colnames(comm)), taxonomy$comm_name)]
comm <- as.matrix(comm[, tree$tip.label] > 0) * 1L
coords <- read.csv(file.path(data_dir, "benicio2021", "coordinates.csv"), row.names = 1)
coords <- coords[rownames(comm), ]
region <- ifelse(coords$Latitude > -23, "inland", "coast")

km <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  6371 * 2 * asin(sqrt(a))
}
dist_km <- outer(seq_len(nrow(coords)), seq_len(nrow(coords)), function(i, j)
  km(coords$Latitude[i], coords$Longitude[i], coords$Latitude[j], coords$Longitude[j]))
low <- lower.tri(dist_km)
same <- outer(region, region, "==")[low]
stats_of <- function(m) {
  both <- tcrossprod(m)
  rich <- rowSums(m)
  jac <- (both / (outer(rich, rich, "+") - both))[low]
  c(mean = mean(jac), within_minus_between = mean(jac[same]) - mean(jac[!same]),
    mantel_r = stats::cor(jac, dist_km[low]),
    mantel_r_within = stats::cor(jac[same], dist_km[low][same]))
}

null <- suppressWarnings(damocles_fit(comm, tree, shared = TRUE))
gamma <- null$gamma[1]
mu <- null$mu[1]
cat(sprintf("Shared fit: gamma = %.3f, mu = %.3f (tree height %.1f Myr rescaled to 1)\n",
            gamma, mu, height))
observed <- stats_of(comm)
sims <- t(replicate(n_sim, stats_of(simulate_damocles(tree, gamma, mu, n = nrow(comm)))))
summary <- data.frame(statistic = names(observed), observed = observed,
                      null_median = apply(sims, 2, stats::median),
                      null_lo = apply(sims, 2, stats::quantile, 0.025),
                      null_hi = apply(sims, 2, stats::quantile, 0.975),
                      p_upper = (1 + colSums(sweep(sims, 2, observed, ">="))) / (n_sim + 1),
                      p_lower = (1 + colSums(sweep(sims, 2, observed, "<="))) / (n_sim + 1),
                      row.names = NULL)
options(width = 200)
print(summary, digits = 3)

by_region <- lapply(split(seq_len(nrow(comm)), region), function(i)
  suppressWarnings(damocles_fit(comm[i, , drop = FALSE], tree, shared = TRUE)))
ll <- function(x) attr(x, "loglik_total")
lrt <- 2 * (sum(vapply(by_region, ll, numeric(1))) - ll(null))
regions <- data.frame(region = names(by_region),
                      gamma = vapply(by_region, function(x) x$gamma[1], numeric(1)),
                      mu = vapply(by_region, function(x) x$mu[1], numeric(1)),
                      richness = as.numeric(tapply(rowSums(comm), region, mean)[names(by_region)]))
print(regions, row.names = FALSE, digits = 3)
cat(sprintf("LRT regions vs shared: %.2f, p = %.3g (2 df)\n", lrt,
            stats::pchisq(lrt, 2, lower.tail = FALSE)))
saveRDS(list(gamma = gamma, mu = mu, height = height, observed = observed, sims = sims,
             summary = summary, regions = regions, lrt = lrt, region = region,
             distance = dist_km), "data-raw/s5_benicio_spatial.rds")
