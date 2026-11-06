# Type I error and power of the framework of Peres-Neto et al. (2012), following the
# simulations of their Tables 1 to 5. Species respond to the environment with unimodal
# curves. Their stacked test at alpha / 2 is mean(p_site, p_taxa) <= 0.025.
library(comdynamite)
library(parallel)
n_rep <- 1000
n_perm <- 999
cores <- max(1, detectCores() - 2)

rescale <- function(x) -1 + 102 * (x - min(x)) / (max(x) - min(x))
responses <- function(env, optima, tolerance) {
  k <- length(optima)
  height <- runif(k, 0.5, 1)
  width <- rnorm(k, tolerance, 10)
  out <- t(height * exp(-outer(optima, env, "-")^2 / (2 * width^2)))
  colnames(out) <- names(optima)
  out
}
clean <- function(comm, ...) {
  keep <- rowSums(comm) > 0
  c(list(comm = comm[keep, colSums(comm[keep, ]) > 0]), lapply(list(...), function(x) x[keep]))
}
# fourth-corner test with the internal functions of phylo_vector()
fourth_corner <- function(comm, trait, env, n_perm) {
  trait <- as.matrix(trait[colnames(comm)])
  env <- as.matrix(env)
  observed <- abs(comdynamite:::corner_correlation(comm, trait, env)[1])
  nulls <- comdynamite:::corner_nulls(comm, trait, env, n_perm, NULL)
  p_site <- comdynamite:::permutation_p(abs(nulls$site), observed)
  p_taxa <- comdynamite:::permutation_p(abs(nulls$taxa), observed)
  list(p_site = p_site, p_taxa = p_taxa, p_max = max(p_site, p_taxa))
}
rejected <- function(x, alpha = 0.05) {
  c(row = x$p_site <= alpha, column = x$p_taxa <= alpha, max = x$p_max <= alpha,
    stacked = mean(c(x$p_site, x$p_taxa)) <= alpha / 2)
}

# Table 4: one trait
trait_case <- function(tolerance, case) {
  env <- runif(100, 0, 100)
  trait <- setNames(runif(50, -1, 101), paste0("sp", 1:50))
  comm <- responses(env, trait, tolerance)
  if (case == "random_environment") env <- rnorm(100)
  if (case == "random_trait") trait[] <- rnorm(50)
  d <- clean(comm, env = env)
  rejected(fourth_corner(d$comm, trait, d$env, n_perm))
}
grid <- expand.grid(tolerance = c(1, 20, 60),
                    case = c("random_trait", "random_environment", "both"),
                    stringsAsFactors = FALSE)
set.seed(2012)
RNGkind("L'Ecuyer-CMRG")
traits <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  runs <- mclapply(seq_len(n_rep), function(r) trait_case(grid$tolerance[i], grid$case[i]),
                   mc.cores = cores)
  cbind(grid[i, ], t(colMeans(do.call(rbind, runs))))
}))

saveRDS(traits, "data-raw/validation_gradient_traits.tmp.rds")

# Tables 1 to 3: phylogenetic vector; `simulated` is the Brownian trait of the species
phylo_case <- function(tolerance, case) {
  tree <- ape::rphylo(50, 1, 0)
  simulated <- rescale(ape::rTraitCont(tree))
  env <- runif(100, 0, 100)
  comm <- responses(env, simulated, tolerance)
  if (case == "random_environment") env <- rnorm(100)
  if (case == "random_phylogeny") tree$tip.label <- sample(tree$tip.label)
  d <- clean(comm, env = env)
  x <- phylo_vector(d$comm, tree, d$env, n_perm = n_perm, n_sim = 1000)
  out <- c(rejected(x), retained = length(x$selected), r_pev = NA, r_brown = NA, r_absolute = NA,
           match_pev = NA, match_brown = NA, match_absolute = NA, final_max = 0)
  if (length(x$selected) > 0) {
    # the code of the authors weights eigenvectors by absolute correlations
    absolute <- drop(x$eigenvectors[, x$selected, drop = FALSE] %*%
                       abs(x$eigenvector_correlation[x$selected]))
    names(absolute) <- names(x$pev)
    species <- names(x$pev)
    out[c("r_pev", "r_brown")] <- abs(x$correlation)
    out["r_absolute"] <- abs(comdynamite:::corner_correlation(d$comm, as.matrix(absolute),
                                                              as.matrix(d$env)))
    out["match_pev"] <- abs(cor(x$pev, simulated[species]))
    out["match_brown"] <- abs(cor(x$pev_brown, simulated[species]))
    out["match_absolute"] <- abs(cor(absolute, simulated[species]))
    best <- if (out[["r_brown"]] >= out[["r_pev"]]) x$pev_brown else x$pev
    final <- fourth_corner(d$comm, best, d$env, n_perm)
    out["final_max"] <- x$p_max <= 0.05 && final$p_max <= 0.05
  }
  out
}
grid <- expand.grid(tolerance = c(1, 20, 60),
                    case = c("random_environment", "random_phylogeny", "both"),
                    stringsAsFactors = FALSE)
phylo <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  runs <- mclapply(seq_len(n_rep), function(r) phylo_case(grid$tolerance[i], grid$case[i]),
                   mc.cores = cores)
  cbind(grid[i, ], t(colMeans(do.call(rbind, runs), na.rm = TRUE)))
}))

saveRDS(phylo, "data-raw/validation_gradient_phylo.tmp.rds")

# The demonstration script: environment and space on a transect of 100 sites (Table 5 is in
# validation_table5.R)
partition_case <- function(case) {
  tree <- ape::stree(32, "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  tree$tip.label <- paste0("sp", 1:32)
  env <- runif(100, 0, 100)
  maps <- as.matrix(adespatial::dbmem(cbind(1:100), silent = TRUE))
  space <- 100 * (maps[, 1] - min(maps[, 1])) / diff(range(maps[, 1]))
  comm <- responses(env, rescale(ape::rTraitCont(tree)), 1) *
    responses(space, rescale(ape::rTraitCont(tree)), 1)
  if (case == "neither") {
    env <- rnorm(100)
    space <- rnorm(100)
  }
  d <- clean(comm, env = env, space = space)
  # the test of the partition on its own uses the vector that is fixed in advance
  fixed <- setNames(rnorm(32), tree$tip.label)[colnames(d$comm)]
  if (case == "both") fixed <- phylo_vector(d$comm, tree, d$env, n_perm = n_perm)$pev_brown
  if (is.null(fixed)) return(rep(NA, 12))
  x <- trait_partition(d$comm, fixed, d$env, d$space, n_perm = n_perm)
  f <- x$fractions
  rownames(f) <- f$fraction
  c(share = x$component_share, f[c("abc", "ab", "bc", "a", "b", "c"), "r2_adjusted"],
    f[c("abc", "ab", "bc", "a", "c"), "p_max"] <= 0.05)
}
fields <- c("share", paste0("r2_", c("abc", "ab", "bc", "a", "b", "c")),
            paste0("reject_", c("abc", "ab", "bc", "a", "c")))
partition <- do.call(rbind, lapply(c("neither", "both"), function(case) {
  runs <- mclapply(seq_len(n_rep), function(r) {
    tryCatch(partition_case(case), error = function(e) conditionMessage(e))
  }, mc.cores = cores)
  failed <- vapply(runs, is.character, logical(1))
  print(table(unlist(runs[failed])))
  runs <- do.call(rbind, runs[!failed])
  colnames(runs) <- fields
  data.frame(case = case, failed = sum(failed), completed = mean(!is.na(runs[, 1])),
             t(colMeans(runs, na.rm = TRUE)))
}))

saveRDS(list(traits = traits, phylo = phylo, partition = partition, n_rep = n_rep, n_perm = n_perm),
        "data-raw/validation_gradient.rds")
