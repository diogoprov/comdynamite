#' Principal coordinates of phylogenetic structure
#'
#' Tests whether the phylogenetic composition of sites follows the environment, with the
#' principal coordinates of phylogenetic structure of Duarte et al. (2016), the method of the
#' package PCPS (Debastiani & Duarte 2014). Phylogenetic fuzzy weighting (Pillar & Duarte
#' 2010) spreads the relative abundance of each species over its relatives, which gives the
#' matrix P; the square root of the Bray-Curtis dissimilarity among rows of P gives the
#' principal coordinates. The test runs on all that dissimilarity, as a distance-based
#' regression, and on each axis in `axes`, as a linear regression.
#'
#' Each test has two permutations, as in `PCPS::pcps.sig()` and `PCPS::matrix.p.sig()`: of
#' sites, which asks whether the environment is related to composition, and of species over
#' the tips of the tree, which asks whether that relation depends on the phylogeny. `p_max` is
#' the larger p-value, and the one to read. Eigenvectors are computed again for each
#' permutation of species. Axes are not rotated onto the observed ones, which for one axis
#' at a time changes only its sign and scale and leaves the F statistic unchanged.
#'
#' Free permutations of sites assume sites are exchangeable. When both the environment and the
#' communities are autocorrelated, without any relation between them, the test of sites rejects
#' far too often: on a grid of 10 x 10 sites, at the 5% level, free permutations rejected 29-80%
#' of the time and `p_max` 3-24%. `site_permutations` gives a restricted null model, such as
#' [torus_permutations()] on a complete grid, which the package PCPS does not offer. With fields
#' that wrap around the grid the torus rejected 4-7%, and `p_max` 1-4%. With fields that stop at
#' the edges of the grid it rejected 6-25%, and `p_max` 1-10%, the highest values for fields
#' that are broad relative to the grid, so the torus reduces this error without removing it in
#' a landscape with edges. It also costs power: with a correlation of 0.7 between the predictor
#' and the field behind the communities, the whole-dissimilarity `p_site` fell from 98-100% to
#' 23-79%.
#'
#' `data-raw/validation_pcps.R` checks the statistics against the package PCPS and gives these
#' rates.
#'
#' @param comm Matrix with one row per site and one column per species. Species are matched
#'   to the tips of `tree` by name.
#'   Or an [mc_data()] object, whose parts fill the arguments left out.
#' @param tree A `phylo` object.
#' @param env Predictor, or matrix or data frame of predictors, with one row per site.
#' @param axes Principal coordinates to test, in order of eigenvalue.
#' @param n_perm Number of permutations of sites and of species.
#' @param site_permutations Optional matrix of site indices with one permutation per row, for
#'   example from [torus_permutations()], to use instead of free permutations of sites. Empty
#'   sites are then allowed: they leave the analysis but their predictors can reach occupied
#'   sites.
#' @return List of class `mc_pcps` with `tests`, a data frame with one row for all the
#'   dissimilarity (`whole`) and one for each axis, and the columns `share`, the share of the
#'   positive eigenvalues; `r2`; `statistic`, the F statistic; `p_site`, `p_species` and
#'   `p_max`; `scores`, the sites on all the axes with positive eigenvalues; and
#'   `permutations`, the numbers of permutations and whether those of sites were free or
#'   restricted.
#' @references
#' Debastiani, V. J., & Duarte, L. D. S. (2014). PCPS -- an R-package for exploring
#' phylogenetic eigenvectors across metacommunities. *Frontiers of Biogeography*, 6.
#' \doi{10.21425/F5FBG22943}
#'
#' Duarte, L. D. S., Debastiani, V. J., Freitas, A. V. L., & Pillar, V. D. (2016). Dissecting
#' phylogenetic fuzzy weighting: theory and application in metacommunity phylogenetics.
#' *Methods in Ecology and Evolution*, 7, 937--946. \doi{10.1111/2041-210X.12547}
#'
#' Pillar, V. D., & Duarte, L. D. S. (2010). A framework for metacommunity analysis of
#' phylogenetic structure. *Ecology Letters*, 13, 587--596.
#' \doi{10.1111/j.1461-0248.2010.01456.x}
#' @seealso [pcc_test()] and [phylo_beta_test()] for other tests over the whole tree, and
#'   [torus_permutations()], [rigid_permutations()].
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(12)
#' traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100), rmax_range = c(0, 30))
#' env <- seq(0, 100, length.out = 36)
#' comm <- simulate_static(traits, env, incidence = TRUE)$abundance
#' pcps_test(comm, tree, env, n_perm = 49)
#' @export
pcps_test <- function(comm, tree, env, axes = 1:2, n_perm = 999, site_permutations = NULL) {
  if (inherits(comm, "mc_data")) {
    if (missing(tree)) tree <- data_part(comm, "tree")
    if (missing(env)) env <- data_part(comm, "env")
    comm <- comm$comm
  }
  matched <- present_species(check_phylo_comm(comm, tree, !is.null(site_permutations)),
                             minimum = 3)
  comm <- matched$comm
  tree <- matched$tree
  landscape_env <- check_env(env, comm)
  check_permutations(site_permutations, comm)
  draws <- site_draws(comm, n_perm, site_permutations)
  comm <- comm[draws$occupied, , drop = FALSE]
  env <- check_env(landscape_env[draws$occupied, , drop = FALSE], comm)
  stopifnot(nrow(comm) > ncol(env) + 1, n_perm >= 1, is.numeric(axes), length(axes) >= 1,
            axes >= 1, axes == round(axes), !anyDuplicated(axes))

  weights <- comm / rowSums(comm)
  belonging <- fuzzy_belonging(tree)
  gower <- pcps_gower(weights, belonging)
  ordination <- pcps_axes(gower)
  if (max(axes) > ncol(ordination$scores)) {
    stop("`axes` asks for axis ", max(axes), " and there are ", ncol(ordination$scores),
         " with positive eigenvalues", call. = FALSE)
  }
  scores <- ordination$scores[, axes, drop = FALSE]
  observed <- pcps_fit(gower, scores, env)
  null_site <- progress_map(seq_len(nrow(draws$rows)), function(r) {
    pcps_fit(gower, scores, landscape_env[draws$rows[r, ], , drop = FALSE])$statistic
  }, "Sites")
  null_species <- progress_map(seq_len(n_perm), function(r) {
    moved <- sample.int(ncol(comm))
    null_gower <- pcps_gower(weights, belonging[moved, moved])
    null_scores <- pcps_axes(null_gower)$scores
    present <- axes <= ncol(null_scores)
    statistic <- rep(NA_real_, length(axes) + 1)
    statistic[c(TRUE, present)] <- pcps_fit(null_gower, null_scores[, axes[present], drop = FALSE],
                                            env)$statistic
    statistic
  }, "Species")
  p_value <- function(nulls) {
    vapply(seq_along(observed$statistic), function(i) {
      permutation_p(vapply(nulls, `[`, numeric(1), i), observed$statistic[i])
    }, numeric(1))
  }
  p_site <- p_value(null_site)
  p_species <- p_value(null_species)
  tests <- data.frame(test = c("whole", paste0("pcps", axes)),
                      share = c(1, ordination$share[axes]), r2 = observed$r2,
                      statistic = observed$statistic, p_site = p_site, p_species = p_species,
                      p_max = pmax(p_site, p_species))
  structure(list(tests = tests, scores = ordination$scores,
                 permutations = permutation_record(n_perm, site_permutations)),
            class = "mc_pcps")
}

#' @export
print.mc_pcps <- function(x, ...) {
  cat("<mc_pcps> principal coordinates of phylogenetic structure\n")
  permutation_line(x$permutations)
  print(format(x$tests, digits = 3), row.names = FALSE, quote = FALSE)
  invisible(x)
}

# degree of belonging of each species to the others, from the patristic distances: rows
# sum to 1 (Pillar & Duarte 2010)
fuzzy_belonging <- function(tree) {
  distances <- ape::cophenetic.phylo(tree)[tree$tip.label, tree$tip.label]
  similarity <- 1 - distances / max(distances)
  similarity / rowSums(similarity)
}

# Gower-centred matrix of the square root of Bray-Curtis among rows of the matrix P
pcps_gower <- function(weights, belonging) {
  p <- weights %*% belonging
  sums <- rowSums(p)
  half <- -0.5 * as.matrix(stats::dist(p, method = "manhattan")) / outer(sums, sums, "+")
  half <- sweep(half, 1, rowMeans(half))
  unname(sweep(half, 2, colMeans(half)))
}

# axes with positive eigenvalues, scaled by the root of the eigenvalue, and their shares
pcps_axes <- function(gower) {
  eigenvalues <- eigen(gower, symmetric = TRUE)
  keep <- eigenvalues$values > sqrt(.Machine$double.eps) * max(eigenvalues$values)
  values <- eigenvalues$values[keep]
  scores <- sweep(eigenvalues$vectors[, keep, drop = FALSE], 2, sqrt(values), "*")
  colnames(scores) <- paste0("pcps", seq_along(values))
  list(scores = scores, share = values / sum(values))
}

# R2 and F of the distance-based regression of all the dissimilarity, then of each axis;
# NA when a permutation leaves the predictors without variance
pcps_fit <- function(gower, scores, env) {
  n <- nrow(gower)
  fit <- qr(sweep(env, 2, colMeans(env)))
  if (fit$rank == 0) return(list(r2 = NA_real_, statistic = rep(NA_real_, ncol(scores) + 1)))
  basis <- qr.Q(fit)[, seq_len(fit$rank), drop = FALSE]
  explained <- c(sum(basis * (gower %*% basis)), colSums(crossprod(basis, scores)^2))
  total <- c(sum(diag(gower)), colSums(scores^2))
  list(r2 = explained / total,
       statistic = (explained / fit$rank) / ((total - explained) / (n - fit$rank - 1)))
}
