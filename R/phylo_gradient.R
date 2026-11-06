#' Test the association between a trait and the environment
#'
#' Weighted correlation between a site predictor and the mean of a species trait in each
#' site, which is the fourth-corner statistic in the form of Peres-Neto et al. (2012).
#' Permuting sites tests the link between composition and the predictor, and permuting
#' species tests the link between composition and the trait. The association needs both,
#' so the p-value to report is the larger of the two (ter Braak et al. 2012).
#'
#' @param comm Matrix with one row per site and one column per species. Every site needs
#'   at least one individual.
#' @param trait Numeric vector with the value of each species. Species are matched to
#'   `comm` by name, in any order, and those in only one of the two are dropped with a
#'   warning.
#' @param env Numeric vector with the value of each site.
#' @param n_perm Number of permutations of sites and of species.
#' @param site_permutations Optional matrix of site indices with one permutation per
#'   row, for example from [torus_permutations()], to use instead of free permutations.
#' @return List of class `mc_fourth_corner` with `correlation`, `p_site`, `p_taxa` and
#'   `p_max`.
#' @section Original code:
#' The function `TraitEnvCor` of the authors keeps the observed mean and standard
#' deviation when it permutes and stacks the two null distributions, with a significance
#' level of alpha / 2. Here each permuted variable is standardised again. The stacked
#' p-value is close to the mean of the two p-values, so any test it rejects at alpha / 2
#' is also rejected by the larger p-value at alpha.
#' @references
#' Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects of spatial
#' contingency and environmental filtering on metacommunity phylogenetics. *Ecology*,
#' 93, S14--S30. \doi{10.1890/11-0494.1}
#'
#' ter Braak, C. J. F., Cormont, A., & Dray, S. (2012). Improved testing of species
#' traits--environment relationships in the fourth-corner problem. *Ecology*, 93,
#' 1525--1526. \doi{10.1890/12-0126.1}
#' @examples
#' set.seed(1)
#' env <- runif(40, 0, 100)
#' trait <- setNames(runif(15, 0, 100), paste0("sp", 1:15))
#' comm <- outer(env, trait, function(e, t) exp(-(e - t)^2 / 800))
#' colnames(comm) <- names(trait)
#' fourth_corner_test(comm, trait, env, n_perm = 99)
#' @export
fourth_corner_test <- function(comm, trait, env, n_perm = 999, site_permutations = NULL) {
  stopifnot(is.numeric(trait), length(env) == nrow(comm), n_perm >= 1)
  comm <- check_species_comm(comm, names(trait), "`trait`")
  if (!is.null(site_permutations)) stopifnot(ncol(site_permutations) == nrow(comm))
  trait <- as.matrix(trait[colnames(comm)])
  env <- as.matrix(env)
  observed <- corner_correlation(comm, trait, env)[1]
  nulls <- corner_nulls(comm, trait, env, n_perm, site_permutations)
  p_site <- permutation_p(abs(nulls$site), abs(observed))
  p_taxa <- permutation_p(abs(nulls$taxa), abs(observed))
  structure(list(correlation = observed, p_site = p_site, p_taxa = p_taxa,
                 p_max = max(p_site, p_taxa)), class = "mc_fourth_corner")
}

#' @export
print.mc_fourth_corner <- function(x, ...) {
  cat("<mc_fourth_corner> trait-environment correlation\n")
  cat("correlation:", format(x$correlation, digits = 3), "\n")
  cat("p (sites):", format(x$p_site, digits = 3), " p (species):",
      format(x$p_taxa, digits = 3), " p (max):", format(x$p_max, digits = 3), "\n")
  invisible(x)
}

#' Represent a phylogeny by one vector related to a predictor
#'
#' The phylogenetic representation of Peres-Neto et al. (2012). Eigenvectors of the
#' phylogenetic distances are extracted with species weights. A global test relates all
#' of them to the predictor. Only if it is significant, eigenvectors enter in the order of
#' their correlations until the first that is not significant, and they are summed with
#' their correlations as weights (`pev`). Brownian vectors are then added to `pev`, and the
#' sum with the largest correlation is kept (`pev_brown`).
#'
#' The p-values of the global test are the ones to report. The two vectors are chosen to
#' fit the predictor, so a test of them with [fourth_corner_test()] is not valid on its
#' own. They are meant for [trait_partition()] after a significant global test.
#'
#' Species are matched to the tips of `tree` by name, in any order, and those in only one
#' of the two are dropped with a warning. Species without individuals have no weight and
#' are left out.
#'
#' @inheritParams fourth_corner_test
#' @param tree A `phylo` object.
#' @param predictor Numeric vector with the value of each site. The authors estimate one
#'   vector for each predictor.
#' @param alpha Significance level of the global test and of the selection.
#' @param n_sim Number of Brownian vectors.
#' @return List of class `mc_phylo_vector` with `correlation`, the correlations of `pev`
#'   and `pev_brown`; `p_site`, `p_taxa` and `p_max` of the global test; `selected`, the
#'   columns of `eigenvectors` retained; `pev` and `pev_brown`, named vectors or `NULL`
#'   when the global test is not significant or no eigenvector is retained;
#'   `eigenvectors`; `eigenvalues`; and `eigenvector_correlation`.
#' @section Article and original code:
#' The code and the text differ. Here we follow the text in the first three points and
#' the code in the last:
#' \itemize{
#'   \item The text sums eigenvectors weighted by their correlations. The code uses the
#'     absolute correlations, and eigenvectors with opposite signs then cancel.
#'   \item The text uses the square root of the path length as distance. The code example
#'     uses one minus the covariance.
#'   \item The text draws 1000 Brownian vectors, and the code 100.
#'   \item The text multiplies each Brownian vector by `pev`, and the code adds them.
#' }
#' In addition, `PCoAWeigthed` of the authors is called with weights that do not sum to
#' one and returns coordinates without dividing by the square root of the weights, so its
#' centring is not the one in the text. `DoubleCenterWeigth`, which is not used in their
#' example, does both, and agrees with this function.
#' @references
#' Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects of spatial
#' contingency and environmental filtering on metacommunity phylogenetics. *Ecology*,
#' 93, S14--S30. \doi{10.1890/11-0494.1}
#' @examples
#' set.seed(1)
#' tree <- ape::rcoal(20)
#' traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100),
#'                          rmax_range = c(0, 30))
#' env <- runif(50, 0, 100)
#' comm <- simulate_static(traits, env, incidence = TRUE)$abundance
#' keep <- rowSums(comm) > 0
#' phylo_vector(comm[keep, ], tree, env[keep], n_perm = 99, n_sim = 100)
#' @export
phylo_vector <- function(comm, tree, predictor, n_perm = 999, alpha = 0.05, n_sim = 1000,
                         site_permutations = NULL) {
  matched <- check_phylo_comm(comm, tree)
  comm <- matched$comm
  tree <- matched$tree
  stopifnot(length(predictor) == nrow(comm), n_perm >= 1, n_sim >= 1, alpha > 0, alpha < 1)
  if (!is.null(site_permutations)) stopifnot(ncol(site_permutations) == nrow(comm))
  comm <- comm[, colSums(comm) > 0, drop = FALSE]
  species <- colnames(comm)
  predictor <- as.matrix(predictor)

  # with the square root of path lengths as distance, D^2 / 2 is half the path length
  pcoa <- weighted_pcoa(ape::cophenetic.phylo(tree)[species, species] / 2, colSums(comm))
  observed <- corner_correlation(comm, pcoa$vectors, predictor)[1, ]
  nulls <- corner_nulls(comm, pcoa$vectors, predictor, n_perm, site_permutations)
  global <- sum(observed^2 * pcoa$values)
  p_site <- permutation_p(nulls$site^2 %*% pcoa$values, global)
  p_taxa <- permutation_p(nulls$taxa^2 %*% pcoa$values, global)
  out <- list(correlation = c(pev = NA_real_, pev_brown = NA_real_), p_site = p_site,
              p_taxa = p_taxa, p_max = max(p_site, p_taxa), selected = integer(), pev = NULL,
              pev_brown = NULL, eigenvectors = pcoa$vectors, eigenvalues = pcoa$values,
              eigenvector_correlation = observed)
  if (out$p_max <= alpha) {
    ranked <- order(abs(observed), decreasing = TRUE)
    p_each <- pmax(vapply(ranked, function(j) permutation_p(abs(nulls$site[, j]), abs(observed[j])),
                          numeric(1)),
                   vapply(ranked, function(j) permutation_p(abs(nulls$taxa[, j]), abs(observed[j])),
                          numeric(1)))
    out$selected <- ranked[cumsum(p_each > alpha) == 0]
  }
  if (length(out$selected) > 0) {
    pev <- drop(pcoa$vectors[, out$selected, drop = FALSE] %*% observed[out$selected])
    candidates <- pev + brownian_vectors(ape::vcv(tree)[species, species], n_sim)
    fit <- corner_correlation(comm, candidates, predictor)[1, ]
    best <- which.max(abs(fit))
    out$pev <- stats::setNames(pev, species)
    out$pev_brown <- stats::setNames(candidates[, best], species)
    out$correlation <- c(pev = corner_correlation(comm, as.matrix(pev), predictor)[1],
                         pev_brown = fit[best])
  }
  structure(out, class = "mc_phylo_vector")
}

#' @export
print.mc_phylo_vector <- function(x, ...) {
  cat("<mc_phylo_vector>", ncol(x$eigenvectors), "eigenvectors,", length(x$selected), "retained\n")
  cat("global test: p (sites):", format(x$p_site, digits = 3), " p (species):",
      format(x$p_taxa, digits = 3), " p (max):", format(x$p_max, digits = 3), "\n")
  if (length(x$selected) > 0) {
    cat("correlation: PEV", format(x$correlation[["pev"]], digits = 3), " PEV-B",
        format(x$correlation[["pev_brown"]], digits = 3), "\n")
  }
  invisible(x)
}

#' Partition trait or phylogenetic variation between environment and space
#'
#' The phylogenetic gradient analysis of Peres-Neto et al. (2012). Trait values are
#' standardised with species weights, which gives each site a mean and a standard
#' deviation. The chosen component is related to `env` and `space` with site weights, and
#' the variation is split into the fractions of Peres-Neto et al. (2006): `a`, environment
#' alone; `b`, shared; `c`, space alone; and `d`, unexplained.
#'
#' The total variation includes both the mean and the dispersion within sites, so `r2`
#' answers how much of all the variation each fraction holds, and `r2_component` how much
#' of the component. Fractions `b` and `d` cannot be tested.
#'
#' @inheritParams fourth_corner_test
#' @param traits Numeric vector or matrix with one row per species, which can hold traits
#'   or the vectors of [phylo_vector()]. Species are matched to `comm` by name, in any
#'   order, and those in only one of the two are dropped with a warning.
#' @param env,space Numeric vector, or matrix or data frame, with one row per site, for
#'   example the output of [spatial_eigenvectors()] for `space`.
#' @param component `"mean"` for the composition of sites or `"variance"` for their
#'   dispersion.
#' @return List of class `mc_trait_partition` with `fractions`, a data frame with one row
#'   for each of `abc`, `ab`, `bc`, `a`, `b`, `c` and `d`, and the columns `r2`,
#'   `r2_adjusted`, `r2_component`, `p_site`, `p_taxa` and `p_max`; `component`; and
#'   `component_share`, the share of the total variation in the component.
#' @section Original code:
#' Results for the mean component agree with `PartitionRDATwoDataMatricesPhylogeny` of the
#' authors in the raw and adjusted fractions of the total. The code departs from the text
#' in the following:
#' \itemize{
#'   \item The fractions relative to the component mix the two denominators, and their
#'     adjusted values use the number of predictors of both sets in every model, which
#'     can give values above one.
#'   \item The total variation does not grow with the number of traits, as it should
#'     according to Appendix A.
#'   \item For the variance, the explained sum of squares contains the intercept, so
#'     `b` absorbs the mean dispersion. Here the intercept is removed.
#'   \item The permutations of the variance take the square root twice.
#' }
#' @references
#' Peres-Neto, P. R., Legendre, P., Dray, S., & Borcard, D. (2006). Variation partitioning
#' of species data matrices: estimation and comparison of fractions. *Ecology*, 87,
#' 2614--2625. \doi{10.1890/0012-9658(2006)87[2614:VPOSDM]2.0.CO;2}
#'
#' Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects of spatial
#' contingency and environmental filtering on metacommunity phylogenetics. *Ecology*,
#' 93, S14--S30. \doi{10.1890/11-0494.1}
#' @examples
#' set.seed(1)
#' env <- runif(40, 0, 100)
#' space <- spatial_eigenvectors(cbind(1:40))[, 1:3]
#' trait <- setNames(runif(15, 0, 100), paste0("sp", 1:15))
#' comm <- outer(env, trait, function(e, t) exp(-(e - t)^2 / 800))
#' colnames(comm) <- names(trait)
#' trait_partition(comm, trait, env, space, n_perm = 49)
#' @export
trait_partition <- function(comm, traits, env, space, component = c("mean", "variance"),
                            n_perm = 999) {
  component <- match.arg(component)
  traits <- as.matrix(traits)
  env <- as.matrix(env)
  space <- as.matrix(space)
  stopifnot(is.numeric(traits), nrow(env) == nrow(comm), nrow(space) == nrow(comm), n_perm >= 1)
  comm <- check_species_comm(comm, rownames(traits), "`traits`")
  traits <- traits[colnames(comm), , drop = FALSE]
  site_weights <- rowSums(comm)
  n <- nrow(comm)
  sizes <- c(abc = ncol(env) + ncol(space), ab = ncol(env), bc = ncol(space))
  total <- ncol(traits) * sum(comm)

  components <- function(response, env, space) {
    env <- weighted_standardise(env, site_weights)
    space <- weighted_standardise(space, site_weights)
    c(abc = explained_ss(cbind(env, space), response, site_weights),
      ab = explained_ss(env, response, site_weights),
      bc = explained_ss(space, response, site_weights))
  }
  tested <- function(ss) {
    r2 <- ss / total
    c(r2, a = r2[["abc"]] - r2[["bc"]], c = r2[["abc"]] - r2[["ab"]])
  }

  response <- trait_response(comm, traits, component)
  component_ss <- sum(site_weights * sweep(response, 2, colSums(site_weights * response) /
                                             sum(site_weights))^2)
  share <- sum(site_weights * response^2) / total
  ss <- components(response, env, space)
  observed <- tested(ss)
  null_site <- replicate(n_perm, {
    i <- sample.int(n)
    tested(components(response, env[i, , drop = FALSE], space[i, , drop = FALSE]))
  })
  null_taxa <- replicate(n_perm, {
    shuffled <- traits[sample.int(nrow(traits)), , drop = FALSE]
    tested(components(trait_response(comm, shuffled, component), env, space))
  })
  p_site <- vapply(names(observed), function(f) permutation_p(null_site[f, ], observed[[f]]),
                   numeric(1))
  p_taxa <- vapply(names(observed), function(f) permutation_p(null_taxa[f, ], observed[[f]]),
                   numeric(1))

  structure(list(fractions = partition_table(ss, total, component_ss, sizes, n, p_site, p_taxa),
                 component = component, component_share = share),
            class = "mc_trait_partition")
}

#' @export
print.mc_trait_partition <- function(x, ...) {
  cat("<mc_trait_partition>", x$component, "component,", format(x$component_share, digits = 3),
      "of the total variation\n")
  print(format(x$fractions, digits = 3), quote = FALSE)
  invisible(x)
}

#' Moran's eigenvector maps
#'
#' Spatial predictors of Dray et al. (2006). Distances above the largest edge of the
#' minimum spanning tree are truncated, and the eigenvectors of the centred similarity
#' matrix with positive Moran's I are kept, in order of decreasing spatial scale.
#'
#' Appendix B of Peres-Neto et al. (2012) recommends centring with site weights so that
#' the maps stay orthogonal in a weighted regression. The function `MEM` of the authors
#' does not, and neither does this one. [trait_partition()] standardises the maps with
#' site weights, and orthogonality is not needed to partition variation.
#'
#' @param coords Matrix of coordinates with one row per site, or a `dist` object.
#' @return Matrix with one row per site and one column per map, with their Moran's I in
#'   the attribute `moran`.
#' @references
#' Dray, S., Legendre, P., & Peres-Neto, P. R. (2006). Spatial modelling: a comprehensive
#' framework for principal coordinate analysis of neighbour matrices (PCNM). *Ecological
#' Modelling*, 196, 483--493. \doi{10.1016/j.ecolmodel.2006.02.015}
#' @examples
#' dim(spatial_eigenvectors(cbind(1:20)))
#' @export
spatial_eigenvectors <- function(coords) {
  distances <- as.matrix(if (inherits(coords, "dist")) coords else stats::dist(coords))
  n <- nrow(distances)
  stopifnot(n > 2)
  threshold <- max(distances[ape::mst(distances) > 0])
  distances[distances > threshold] <- 4 * threshold
  similarity <- 1 - (distances / max(distances))^2
  diag(similarity) <- 0
  centring <- diag(n) - 1 / n
  eig <- eigen(centring %*% similarity %*% centring, symmetric = TRUE)
  moran <- n / sum(similarity) * eig$values
  keep <- moran > 1e-14
  structure(eig$vectors[, keep, drop = FALSE], moran = moran[keep])
}

# weighted correlations between the columns of env and the site means of the columns of
# traits; species without individuals have no weight
corner_correlation <- function(comm, traits, env) {
  site_weights <- rowSums(comm)
  means <- comm %*% weighted_standardise(traits, colSums(comm)) / site_weights
  crossprod(weighted_standardise(env, site_weights), site_weights * means) / sum(site_weights)
}

# one row per permutation and one column per column of traits
corner_nulls <- function(comm, traits, env, n_perm, site_permutations) {
  sites <- if (is.null(site_permutations)) {
    t(replicate(n_perm, sample.int(nrow(comm))))
  } else {
    site_permutations
  }
  permuted <- matrix(env[as.vector(t(sites)), 1], nrow(comm))
  list(site = corner_correlation(comm, traits, permuted),
       taxa = t(replicate(n_perm, {
         corner_correlation(comm, traits[sample.int(nrow(traits)), , drop = FALSE], env)[1, ]
       })) |> matrix(n_perm))
}

permutation_p <- function(null, observed) (sum(null >= observed) + 1) / (length(null) + 1)

# principal coordinates with species weights (Appendix B of Peres-Neto et al. 2012);
# half_squared holds D^2 / 2. The coordinates have weighted mean zero and are orthogonal
# under the weights
weighted_pcoa <- function(half_squared, weights) {
  w <- weights / sum(weights)
  n <- length(w)
  projector <- diag(n) - matrix(w, n, n, byrow = TRUE)
  centred <- projector %*% half_squared %*% t(projector)
  root <- sqrt(w)
  eig <- eigen(-outer(root, root) * centred, symmetric = TRUE)
  keep <- eig$values > sqrt(.Machine$double.eps) * max(eig$values)
  list(vectors = sweep(eig$vectors[, keep, drop = FALSE], 2, sqrt(eig$values[keep]), "*") / root,
       values = eig$values[keep])
}

# one Brownian vector per column
brownian_vectors <- function(covariance, n_sim) {
  eig <- eigen(covariance, symmetric = TRUE)
  eig$vectors %*% (sqrt(pmax(eig$values, 0)) * matrix(stats::rnorm(nrow(covariance) * n_sim),
                                                      nrow(covariance)))
}

# site means or standard deviations of the standardised traits (Appendix A)
trait_response <- function(comm, traits, component) {
  site_weights <- rowSums(comm)
  standardised <- weighted_standardise(traits, colSums(comm))
  means <- comm %*% standardised / site_weights
  if (component == "mean") return(means)
  sqrt(pmax(comm %*% standardised^2 / site_weights - means^2, 0))
}

# weighted sum of squares explained by centred predictors, around the weighted mean
explained_ss <- function(predictors, response, site_weights) {
  root <- sqrt(site_weights)
  centred <- sweep(response, 2, colSums(site_weights * response) / sum(site_weights))
  sum(qr.fitted(qr(root * predictors), root * centred)^2)
}

partition_table <- function(ss, total, component_ss, sizes, n, p_site, p_taxa) {
  fractions <- function(x) {
    c(x, a = x[["abc"]] - x[["bc"]], b = x[["ab"]] + x[["bc"]] - x[["abc"]],
      c = x[["abc"]] - x[["ab"]], d = 1 - x[["abc"]])
  }
  r2 <- ss / total
  table <- data.frame(fraction = c("abc", "ab", "bc", "a", "b", "c", "d"),
                      r2 = fractions(r2),
                      r2_adjusted = fractions(1 - (n - 1) / (n - sizes - 1) * (1 - r2)),
                      r2_component = fractions(ss / component_ss), row.names = NULL)
  table$p_site <- unname(p_site[table$fraction])
  table$p_taxa <- unname(p_taxa[table$fraction])
  table$p_max <- pmax(table$p_site, table$p_taxa)
  table
}
