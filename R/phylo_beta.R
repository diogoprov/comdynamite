#' Phylogenetic beta diversity between communities
#'
#' Mean phylogenetic distance between the species of each pair of communities
#' (COMDIST), the metric used by Peres-Neto & Kembel (2015). Values in `comm` are
#' the weights of the species, so incidence data give the unweighted mean and
#' abundances give the abundance-weighted mean.
#'
#' Two communities with the same species have a distance equal to their mean
#' pairwise distance, not zero, whereas the distance of a community to itself is
#' set to zero. Peres-Neto & Kembel (2015) describe the matrix as fully Euclidean,
#' but this does not always hold. If sites a and b both have species x and y, site
#' c has only x and site d has only y, then a and b lie halfway between c and d
#' and should coincide, yet their distance is not zero. The pseudo-F of
#' [phylo_beta_test()] does not require a Euclidean matrix (McArdle & Anderson
#' 2001).
#'
#' @param comm Matrix with one row per site and one column per species, named
#'   as the tips of `tree`. Every site needs at least one individual.
#' @param tree A `phylo` object with the same species as `comm`.
#' @return A `dist` object.
#' @references
#' McArdle, B. H., & Anderson, M. J. (2001). Fitting multivariate models to
#' community data: a comment on distance-based redundancy analysis. *Ecology*,
#' 82, 290--297. \doi{10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2}
#'
#' Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient analysis:
#' environmental drivers of phylogenetic variation across ecological
#' communities. *Plant Ecology*, 216, 709--724.
#' \doi{10.1007/s11258-014-0405-0}
#' @examples
#' set.seed(1)
#' tree <- ape::rcoal(5)
#' comm <- matrix(rbinom(20, 1, 0.7) + diag(1, 4, 5), 4,
#'                dimnames = list(NULL, tree$tip.label))
#' phylo_beta(comm > 0, tree)
#' @export
phylo_beta <- function(comm, tree) {
  matched <- check_phylo_comm(comm, tree)
  stats::as.dist(beta_matrix(matched$comm / rowSums(matched$comm), tip_distances(matched$tree)))
}

#' Test phylogenetic beta diversity against the environment
#'
#' Implements the phylogenetic beta-diversity framework of Peres-Neto & Kembel
#' (2015), a distance-based redundancy analysis (Legendre & Anderson 1999;
#' McArdle & Anderson 2001) of [phylo_beta()] on the predictors. The test
#' statistic is the pseudo-F. It is tested by permuting sites, which breaks the
#' link between composition and environment, and by permuting species across the
#' tips, which breaks the link between composition and phylogeny. The association
#' requires both, so the p-value to report is the larger of the two.
#'
#' The equations printed in the article omit the 1 in the adjusted R2, the
#' transpose in the hat matrix and the centring of the matrix that is regressed.
#' Here the predictors and the matrix are centred, and the residual degrees of
#' freedom are `n - p - 1`, with `n - 1` for the total, instead of the printed
#' `n - p` and `n`. The pseudo-F and R2 therefore equal those of a distance-based
#' redundancy analysis in vegan. The printed degrees of freedom would multiply
#' the pseudo-F by a constant, which does not change the permutation tests.
#'
#' The authors used presences and absences. Abundances are accepted, as weights
#' in [phylo_beta()].
#'
#' @inheritParams phylo_beta
#' @param env Predictor, or matrix or data frame of predictors, with one row per
#'   site.
#' @param n_perm Number of permutations of sites and of species.
#' @param site_permutations Optional matrix of site indices with one permutation
#'   per row, for example from [torus_permutations()], to use instead of free
#'   permutations of sites.
#' @return List of class `mc_phylo_beta` with `r2`, the adjusted R2; `statistic`,
#'   the pseudo-F; `p_site`, `p_taxa` and `p_max`; and `scores`, the sites on the
#'   principal coordinates of the predicted matrix, which the authors suggest
#'   using to order sites in a heatmap of the community matrix.
#' @references
#' Legendre, P., & Anderson, M. J. (1999). Distance-based redundancy analysis:
#' testing multispecies responses in multifactorial ecological experiments.
#' *Ecological Monographs*, 69, 1--24.
#' \doi{10.1890/0012-9615(1999)069[0001:DBRATM]2.0.CO;2}
#'
#' McArdle, B. H., & Anderson, M. J. (2001). Fitting multivariate models to
#' community data: a comment on distance-based redundancy analysis. *Ecology*,
#' 82, 290--297. \doi{10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2}
#'
#' Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient analysis:
#' environmental drivers of phylogenetic variation across ecological
#' communities. *Plant Ecology*, 216, 709--724.
#' \doi{10.1007/s11258-014-0405-0}
#' @seealso [pcc_test()], the other framework of the same article.
#' @examples
#' set.seed(1)
#' tree <- ape::rcoal(20)
#' traits <- species_traits(tree, lambda = 1, optima_range = c(-1, 101),
#'                          rmax_range = c(0, 30))
#' env <- runif(50, 0, 100)
#' comm <- simulate_static(traits, env, incidence = TRUE)$abundance
#' keep <- rowSums(comm) > 0
#' phylo_beta_test(comm[keep, ], tree, env[keep], n_perm = 99)
#' @export
phylo_beta_test <- function(comm, tree, env, n_perm = 999, site_permutations = NULL) {
  matched <- check_phylo_comm(comm, tree)
  comm <- matched$comm
  tree <- matched$tree
  env <- as.matrix(env)
  stopifnot(nrow(env) == nrow(comm), nrow(comm) > ncol(env) + 1, n_perm >= 1)
  if (!is.null(site_permutations)) stopifnot(ncol(site_permutations) == nrow(comm))
  decomposition <- qr(sweep(env, 2, colMeans(env)))
  if (decomposition$rank < ncol(env)) stop("predictors are collinear", call. = FALSE)
  basis <- qr.Q(decomposition)
  weights <- comm / rowSums(comm)
  distances <- tip_distances(tree)
  centred <- gower_centre(beta_matrix(weights, distances))
  fit <- beta_fit(centred, basis)

  sites <- if (is.null(site_permutations)) {
    t(replicate(n_perm, sample.int(nrow(comm))))
  } else {
    site_permutations
  }
  null_site <- apply(sites, 1, function(i) {
    beta_fit(centred, basis[i, , drop = FALSE])$statistic
  })
  null_taxa <- replicate(n_perm, {
    shuffled <- weights[, sample.int(ncol(weights)), drop = FALSE]
    beta_fit(gower_centre(beta_matrix(shuffled, distances)), basis)$statistic
  })
  p_site <- (sum(null_site >= fit$statistic) + 1) / (length(null_site) + 1)
  p_taxa <- (sum(null_taxa >= fit$statistic) + 1) / (length(null_taxa) + 1)
  structure(list(r2 = fit$r2, statistic = fit$statistic, p_site = p_site, p_taxa = p_taxa,
                 p_max = max(p_site, p_taxa), scores = predicted_scores(centred, basis)),
            class = "mc_phylo_beta")
}

#' @export
print.mc_phylo_beta <- function(x, ...) {
  cat("<mc_phylo_beta> phylogenetic beta diversity\n")
  cat("adjusted R2:", format(x$r2, digits = 3), " pseudo-F:",
      format(x$statistic, digits = 3), "\n")
  cat("p (sites):", format(x$p_site, digits = 3), " p (species):",
      format(x$p_taxa, digits = 3), " p (max):", format(x$p_max, digits = 3), "\n")
  invisible(x)
}

# the community and the tree with the same species, in the order of the tips
check_phylo_comm <- function(comm, tree) {
  stopifnot(inherits(tree, "phylo"))
  comm <- check_species_comm(comm, tree$tip.label, "the phylogeny")
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, colnames(comm)))
  list(comm = comm[, tree$tip.label, drop = FALSE], tree = tree)
}

# numeric community with the species in `species`, in that order, and no empty site
check_species_comm <- function(comm, species, source) {
  comm <- as.matrix(comm) * 1
  stopifnot(all(comm >= 0))
  comm <- comm[, match_species(comm, species, source), drop = FALSE]
  if (any(rowSums(comm) == 0)) stop("every site needs at least one individual", call. = FALSE)
  comm
}

# species shared by `comm` and `species`, matched by name, as match.phylo.comm() of
# picante does. Species of `comm` without individuals do not matter, so they raise no
# warning
match_species <- function(comm, species, source) {
  if (is.null(colnames(comm))) stop("`comm` needs species names in its columns", call. = FALSE)
  if (is.null(species)) stop(source, " needs species names", call. = FALSE)
  repeated <- unique(c(colnames(comm)[duplicated(colnames(comm))], species[duplicated(species)]))
  if (length(repeated) > 0) stop("repeated species: ", species_list(repeated), call. = FALSE)
  only_comm <- setdiff(colnames(comm)[colSums(comm) > 0], species)
  only_source <- setdiff(species, colnames(comm))
  if (length(only_comm) > 0) {
    warning(length(only_comm), " species of `comm` not in ", source, " dropped: ",
            species_list(only_comm), call. = FALSE)
  }
  if (length(only_source) > 0) {
    warning(length(only_source), " species of ", source, " not in `comm` dropped: ",
            species_list(only_source), call. = FALSE)
  }
  shared <- intersect(species, colnames(comm))
  if (length(shared) < 2) stop("`comm` and ", source, " share fewer than two species", call. = FALSE)
  shared
}

species_list <- function(x) {
  shown <- paste(x[seq_len(min(length(x), 10))], collapse = ", ")
  if (length(x) > 10) paste0(shown, " and ", length(x) - 10, " more") else shown
}

tip_distances <- function(tree) {
  ape::cophenetic.phylo(tree)[tree$tip.label, tree$tip.label]
}

# weights sum to 1 in each row, so each cell is the mean distance between two sites
beta_matrix <- function(weights, distances) {
  beta <- weights %*% distances %*% t(weights)
  diag(beta) <- 0
  beta
}

gower_centre <- function(beta) {
  a <- -0.5 * beta^2
  a <- sweep(a, 1, rowMeans(a))
  sweep(a, 2, colMeans(a))
}

# basis is orthonormal, so trace(HGH) = trace(Q'GQ) and the residual trace is the
# difference from the trace of G
beta_fit <- function(centred, basis) {
  n <- nrow(centred)
  p <- ncol(basis)
  total <- sum(diag(centred))
  fitted <- sum(basis * (centred %*% basis))
  residual <- total - fitted
  list(statistic = (fitted / p) / (residual / (n - p - 1)),
       r2 = 1 - (residual / (n - p - 1)) / (total / (n - 1)))
}

predicted_scores <- function(centred, basis) {
  hat <- tcrossprod(basis)
  eig <- eigen(hat %*% centred %*% hat, symmetric = TRUE)
  keep <- eig$values > sqrt(.Machine$double.eps) * max(eig$values)
  sweep(eig$vectors[, keep, drop = FALSE], 2, sqrt(eig$values[keep]), "*")
}
