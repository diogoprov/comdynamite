#' Phylogenetic eigenvectors
#'
#' Principal coordinates of the patristic distances among species, the first
#' step of the phylogenetic community composition of Peres-Neto & Kembel (2015).
#'
#' @param tree A `phylo` object.
#' @param squared Use `-PD^2 / 2`, as in a classic principal coordinate analysis,
#'   instead of `-PD / 2`, the form printed by Peres-Neto & Kembel (2015)?
#' @return Matrix with one row per tip and one column per positive eigenvalue.
#' @references
#' Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient analysis:
#' environmental drivers of phylogenetic variation across ecological
#' communities. *Plant Ecology*, 216, 709--724.
#' \doi{10.1007/s11258-014-0405-0}
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' dim(phylo_eigenvectors(treeshapesim::sim_tree_shape(10)))
#' @export
phylo_eigenvectors <- function(tree, squared = FALSE) {
  stopifnot(inherits(tree, "phylo"))
  pd <- ape::cophenetic.phylo(tree)[tree$tip.label, tree$tip.label]
  b <- -0.5 * if (squared) pd^2 else pd
  b <- sweep(b, 1, rowMeans(b))
  b <- sweep(b, 2, colMeans(b))
  eig <- eigen(b, symmetric = TRUE)
  keep <- eig$values > sqrt(.Machine$double.eps) * max(eig$values)
  vectors <- eig$vectors[, keep, drop = FALSE]
  rownames(vectors) <- tree$tip.label
  vectors
}

#' Test phylogenetic community composition against the environment
#'
#' Implements the phylogenetic community composition (PCC) framework of
#' Peres-Neto & Kembel (2015). Phylogenetic eigenvectors are standardised with
#' species weights and averaged within sites, and each resulting variable is
#' regressed on the predictors with site weights. The test statistic is the
#' variance of the adjusted R2 across eigenvectors. It is tested by permuting
#' sites, which breaks the link between composition and environment, and by
#' permuting species across the tips, which breaks the link between composition
#' and phylogeny. The association requires both, so the p-value to report is
#' the larger of the two.
#'
#' The authors used presences and absences. Abundances are accepted, with
#' species and site totals as weights. Species without individuals leave the analysis.
#'
#' @param comm Matrix with one row per site and one column per species. Species are
#'   matched to the tips of `tree` by name, in any order, and those in only one of the
#'   two are dropped with a warning. Every site needs at least one individual.
#' @param tree A `phylo` object.
#' @param env Predictor, or matrix or data frame of predictors, with one row per
#'   site.
#' @param n_perm Number of permutations of sites and of species.
#' @param squared Passed to [phylo_eigenvectors()].
#' @param site_permutations Optional matrix of site indices with one permutation
#'   per row, for example from [torus_permutations()], to use instead of free
#'   permutations of sites. Empty sites are then allowed, because they are part of the
#'   landscape: they leave the analysis but their predictors can reach occupied sites.
#' @return List of class `mc_pcc` with `r2_mean`, the mean of the positive adjusted R2
#'   values, which is the value the authors recommend reporting, or 0 if none is positive;
#'   `statistic`, the
#'   variance of the adjusted R2 across eigenvectors; `p_site`, `p_species` and
#'   `p_max`; `r2_adjusted`, with one value per eigenvector; and `permutations`.
#' @references
#' Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient analysis:
#' environmental drivers of phylogenetic variation across ecological
#' communities. *Plant Ecology*, 216, 709--724.
#' \doi{10.1007/s11258-014-0405-0}
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(20)
#' traits <- species_traits(tree, lambda = 1, optima_range = c(-1, 101),
#'                          rmax_range = c(0, 30))
#' env <- runif(50, 0, 100)
#' comm <- simulate_static(traits, env, incidence = TRUE)$abundance
#' keep <- rowSums(comm) > 0
#' pcc_test(comm[keep, ], tree, env[keep], n_perm = 99)
#' @export
pcc_test <- function(comm, tree, env, n_perm = 999, squared = FALSE,
                     site_permutations = NULL) {
  matched <- present_species(check_phylo_comm(comm, tree, !is.null(site_permutations)),
                             minimum = 3)
  comm <- matched$comm
  tree <- matched$tree
  landscape_env <- check_env(env, comm)
  check_permutations(site_permutations, comm)
  draws <- site_draws(comm, n_perm, site_permutations)
  comm <- comm[draws$occupied, , drop = FALSE]
  env <- check_env(landscape_env[draws$occupied, , drop = FALSE], comm)
  stopifnot(nrow(comm) > ncol(env) + 1, n_perm >= 1)
  vectors <- phylo_eigenvectors(tree, squared = squared)
  observed <- pcc_r2(comm, vectors, env)
  statistic <- stats::var(observed)

  null_site <- unlist(progress_map(seq_len(nrow(draws$rows)), function(r) {
    stats::var(pcc_r2(comm, vectors, landscape_env[draws$rows[r, ], , drop = FALSE]))
  }, "Sites"))
  null_taxa <- unlist(progress_map(seq_len(n_perm), function(r) {
    stats::var(pcc_r2(comm[, sample.int(ncol(comm)), drop = FALSE], vectors, env))
  }, "Species"))
  p_site <- (sum(null_site >= statistic) + 1) / (length(null_site) + 1)
  p_species <- (sum(null_taxa >= statistic) + 1) / (length(null_taxa) + 1)
  structure(list(r2_mean = sum(pmax(observed, 0)) / max(sum(observed > 0), 1),
                 statistic = statistic, p_site = p_site,
                 p_species = p_species, p_max = max(p_site, p_species), r2_adjusted = observed,
                 permutations = permutation_record(n_perm, site_permutations)),
            class = "mc_pcc")
}

#' @export
print.mc_pcc <- function(x, ...) {
  cat("<mc_pcc> phylogenetic community composition\n")
  cat("mean positive adjusted R2:", format(x$r2_mean, digits = 3), "\n")
  permutation_line(x$permutations)
  cat("p (sites):", format(x$p_site, digits = 3), " p (species):",
      format(x$p_species, digits = 3), " p (max):", format(x$p_max, digits = 3), "\n")
  invisible(x)
}

# adjusted R2 of each eigenvector; the total sum of squares is the sum of the
# species weights, because the standardised eigenvectors have weighted variance 1
pcc_r2 <- function(comm, vectors, env) {
  n <- nrow(comm)
  species_weights <- colSums(comm)
  site_weights <- rowSums(comm)
  composition <- (comm %*% weighted_standardise(vectors, species_weights)) / site_weights
  predictors <- weighted_standardise(env, site_weights)
  slopes <- solve(crossprod(predictors, site_weights * predictors),
                  crossprod(predictors, site_weights * composition))
  fitted <- predictors %*% slopes
  r2 <- colSums(site_weights * fitted^2) / sum(species_weights)
  1 - (n - 1) / (n - ncol(env) - 1) * (1 - r2)
}

# weighted mean 0 and variance 1, with the sum of the weights as denominator
weighted_standardise <- function(x, w) {
  centred <- sweep(x, 2, colSums(w * x) / sum(w))
  sweep(centred, 2, sqrt(colSums(w * centred^2) / sum(w)), "/")
}
