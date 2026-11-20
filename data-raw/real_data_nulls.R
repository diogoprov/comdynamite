# Block A, item 1: the restricted nulls with edges applied to the two example data sets.
# - daphniids (314 lakes of New England, Leibold et al. 2010): node_test() with the
#   biogeographic events, with free permutations, with patch_permutations() for the events, and
#   with variogram and parametric surrogates of the predictors; pcps_test() and
#   phylo_beta_test() with free permutations and with both surrogates.
# - anurans (40 ponds in two cities, 20 each): pcps_test() and phylo_beta_test() with free
#   permutations, permutations within each city, and surrogates drawn within each city with
#   great-circle distances.
# Coordinates in longitude and latitude are projected to kilometres (equirectangular, at the
# mean latitude) where a function needs planar coordinates. Run from the package root with the
# package installed; takes a few minutes on 8 cores. Writes data-raw/real_data_nulls.rds.
library(comdynamite)
n_perm <- as.integer(Sys.getenv("N_PERM", "999"))
cores <- as.integer(Sys.getenv("CORES", "8"))
project <- function(lonlat) {
  lat0 <- mean(lonlat[, 2]) * pi / 180
  cbind(x = lonlat[, 1] * 111.32 * cos(lat0), y = lonlat[, 2] * 110.57)
}

# daphniids ------------------------------------------------------------------------------
d <- daphniids
km <- project(d$coords)
nodes <- node_contrasts(d$tree_leibold)
set.seed(1)
node_free <- node_test(d$comm, nodes, d$env, d$biogeography, n_perm = n_perm, cores = cores)
set.seed(1)
areas <- patch_permutations(km, d$biogeography, n_perm = n_perm)
node_patch <- node_test(d$comm, nodes, d$env, d$biogeography, n_perm = n_perm, cores = cores,
                        area_permutations = areas)
set.seed(1)
node_vgram <- node_test(d$comm, nodes, d$env, n_perm = n_perm, cores = cores,
                        null = "variogram", coords = km)
set.seed(1)
node_param <- node_test(d$comm, nodes, d$env, n_perm = n_perm, cores = cores,
                        null = "parametric", coords = km)

count_partial <- function(out) colSums(out$p_partial_site < 0.05, na.rm = TRUE)
daphniid_nodes <- rbind(free = count_partial(node_free), variogram = count_partial(node_vgram),
                        parametric = count_partial(node_param))
daphniid_events <- rbind(free = colSums(node_free$biogeography$p_partial_site < 0.05,
                                        na.rm = TRUE),
                         patches = colSums(node_patch$biogeography$p_partial_site < 0.05,
                                           na.rm = TRUE))

whole_tree <- function(comm, tree, env, nulls) {
  do.call(rbind, lapply(names(nulls), function(nm) {
    set.seed(2)
    pc <- do.call(pcps_test, c(list(comm, tree, env, n_perm = n_perm), nulls[[nm]]))
    set.seed(2)
    pb <- do.call(phylo_beta_test, c(list(comm, tree, env, n_perm = n_perm), nulls[[nm]]))
    data.frame(null = nm, test = c(pc$tests$test, "phylo_beta"),
               r2 = c(pc$tests$r2, pb$r2_adjusted),
               p_site = c(pc$tests$p_site, pb$p_site), p_max = c(pc$tests$p_max, pb$p_max))
  }))
}
occupied <- rowSums(d$comm) > 0
daphniid_tree <- whole_tree(d$comm[occupied, ], d$tree, d$env[occupied, ], list(
  free = list(),
  variogram = list(null = "variogram", coords = km[occupied, ]),
  parametric = list(null = "parametric", coords = km[occupied, ])))

# anurans --------------------------------------------------------------------------------
a <- anurans
keep <- rowSums(a$comm) > 0
comm <- a$comm[keep, ]
env <- scale(a$env[keep, ])
city <- a$city[keep]
coords <- a$coords[keep, ]
set.seed(3)
within_city <- t(replicate(n_perm, {
  rows <- seq_along(city)
  for (g in split(seq_along(city), city)) rows[g] <- g[sample.int(length(g))]
  rows
}))
set.seed(3)
vgram <- env_surrogates(env, coords, n_perm, "variogram", lonlat = TRUE, strata = city)
param <- env_surrogates(env, coords, n_perm, "parametric", lonlat = TRUE, strata = city)
anuran_tree <- whole_tree(comm, a$tree, env, list(
  free = list(), within_city = list(site_permutations = within_city),
  variogram = list(null = vgram), parametric = list(null = param)))

results <- list(daphniid_nodes = daphniid_nodes, daphniid_events = daphniid_events,
                daphniid_tree = daphniid_tree, anuran_tree = anuran_tree,
                n_nodes = ncol(nodes), n_anuran_sites = sum(keep))
saveRDS(results, "data-raw/real_data_nulls.rds", version = 2)
print(results, digits = 3)
