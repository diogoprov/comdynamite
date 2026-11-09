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
#'   `p_max`; `scores`, the sites on all the axes with positive eigenvalues; `share`, the
#'   share of each of those axes; `species`, the
#'   correlation of each species' column of the matrix P with those axes (its position in the
#'   biplot of the package PCPS); `env`, the predictors of the occupied sites; and
#'   `permutations`, the numbers of permutations and whether those of sites were free or
#'   restricted. [plot.mc_pcps()] draws the ordination.
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
#' res <- pcps_test(comm, tree, env, n_perm = 49)
#' res
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(res)
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
  composition <- weights %*% belonging
  species <- suppressWarnings(stats::cor(composition, ordination$scores))
  rownames(species) <- colnames(comm)
  rownames(ordination$scores) <- rownames(comm)
  structure(list(tests = tests, scores = ordination$scores, share = ordination$share,
                 species = species, env = env,
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

#' Plot the principal coordinates of phylogenetic structure
#'
#' Draws the sites on two principal coordinates of [pcps_test()], coloured by a predictor, with
#' the species as arrows from the correlation of their column of the matrix P with the axes,
#' as in the biplot of the package PCPS (Debastiani & Duarte 2014).
#'
#' @param x Result of [pcps_test()].
#' @param axes The two principal coordinates to draw.
#' @param env Predictor that colours the sites: a column number or name of the predictors given
#'   to [pcps_test()], or `NULL` for no colour.
#' @param species Draw the species as arrows.
#' @param species_labels Write the names of the species at the arrows.
#' @param palette Name of a palette of [grDevices::hcl.colors()] for the predictor.
#' @param point_size Size of the site points.
#' @param ... Not used.
#' @return A ggplot object.
#' @seealso [pcps_test()].
#' @references
#' Debastiani, V. J., & Duarte, L. D. S. (2014). PCPS -- an R-package for exploring
#' phylogenetic eigenvectors across metacommunities. *Frontiers of Biogeography*, 6.
#' \doi{10.21425/F5FBG22943}
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(12)
#' traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100), rmax_range = c(0, 30))
#' env <- seq(0, 100, length.out = 36)
#' comm <- simulate_static(traits, env, incidence = TRUE)$abundance
#' plot(pcps_test(comm, tree, env, n_perm = 19))
#' @export
plot.mc_pcps <- function(x, axes = c(1, 2), env = 1, species = TRUE, species_labels = TRUE,
                         palette = "viridis", point_size = 2.5, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot() of pcps_test() needs the package ggplot2", call. = FALSE) # nocov
  }
  if (is.null(x$species)) {
    stop("`x` comes from an older version of pcps_test(); run it again", call. = FALSE)
  }
  n_axes <- ncol(x$scores)
  if (!is.numeric(axes) || length(axes) != 2 || any(axes < 1) || any(axes > n_axes) ||
      any(axes != round(axes)) || axes[1] == axes[2]) {
    stop("`axes` must be two different axes between 1 and ", n_axes, call. = FALSE)
  }
  stopifnot(is.logical(species), length(species) == 1, !is.na(species),
            is.logical(species_labels), length(species_labels) == 1, !is.na(species_labels))
  axis_label <- function(a) sprintf("PCPS %d (%.1f%%)", a, 100 * x$share[a])
  sites <- data.frame(a1 = x$scores[, axes[1]], a2 = x$scores[, axes[2]])
  colour_name <- NULL
  if (!is.null(env)) {
    predictors <- as.data.frame(x$env)
    if (is.null(colnames(x$env))) names(predictors) <- paste0("env", seq_along(predictors))
    if (length(env) != 1 || !(env %in% names(predictors) ||
                              (is.numeric(env) && env >= 1 && env <= ncol(predictors)))) {
      stop("`env` must be one column of the predictors given to pcps_test(), or NULL",
           call. = FALSE)
    }
    colour_name <- if (is.numeric(env)) names(predictors)[env] else env
    sites$env <- predictors[[env]]
  }
  .data <- ggplot2::.data
  p <- ggplot2::ggplot(sites, ggplot2::aes(x = .data$a1, y = .data$a2)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85")
  if (species) {
    # arrows scaled to the spread of the sites, as in a biplot
    scale <- 0.9 * max(abs(as.matrix(sites[c("a1", "a2")]))) /
      max(abs(x$species[, axes]), na.rm = TRUE)
    sp <- data.frame(a1 = x$species[, axes[1]] * scale, a2 = x$species[, axes[2]] * scale,
                     label = rownames(x$species))
    sp <- sp[stats::complete.cases(sp), , drop = FALSE]
    p <- p + ggplot2::geom_segment(data = sp, ggplot2::aes(x = 0, y = 0, xend = .data$a1,
                                                           yend = .data$a2),
                                   colour = "grey50", linewidth = 0.4,
                                   arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")))
    if (species_labels) {
      p <- p + ggplot2::geom_text(data = sp, ggplot2::aes(label = .data$label),
                                  colour = "grey30", size = 3, vjust = -0.4)
    }
  }
  if (is.null(colour_name)) {
    p <- p + ggplot2::geom_point(size = point_size)
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$env), size = point_size) +
      ggplot2::scale_colour_gradientn(colours = grDevices::hcl.colors(30, palette),
                                      name = colour_name)
  }
  p + ggplot2::labs(x = axis_label(axes[1]), y = axis_label(axes[2])) +
    ggplot2::coord_equal() + ggplot2::theme_minimal()
}
