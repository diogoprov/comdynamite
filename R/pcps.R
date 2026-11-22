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
#' 23-79%. In landscapes with edges, `null = "variogram"` or `"parametric"` came closer to the
#' nominal rate (see [env_surrogates()]).
#'
#' `data-raw/validation_pcps.R` checks the statistics against the package PCPS and gives these
#' rates.
#'
#' @param comm Matrix with one row per site and one column per species. Species are matched
#'   to the tips of `tree` by name.
#'   Or an [mc_data()] object, whose parts fill the arguments left out.
#' @param tree A `phylo` object.
#' @param env Predictor, or matrix or data frame of predictors, with one row per site.
#'   Columns must be numeric or logical: factors and character columns stop with an error, so
#'   code a factor as dummy variables first, for example with [stats::model.matrix()].
#' @param axes Principal coordinates to test, in order of eigenvalue: a vector of distinct
#'   whole numbers, 1 or more.
#' @param n_perm Number of permutations of sites and of species. One whole number, 1 or more.
#' @param site_permutations Optional matrix of site indices with one permutation per row, for
#'   example from [torus_permutations()], to use instead of free permutations of sites. Empty
#'   sites are then allowed: they leave the analysis but their predictors can reach occupied
#'   sites.
#' @param null Null model of sites: `"permutation"` (free, or `site_permutations` when given),
#'   or surrogates of the predictors with the same spatial autocorrelation, `"variogram"` or
#'   `"parametric"`, for landscapes with edges, or a list of surrogates from
#'   [env_surrogates()], with one row per site (or per occupied site), for example drawn
#'   within strata (see [env_surrogates()], which gives their
#'   rates of type I error). Surrogates need `coords`; empty sites are then not allowed.
#' @param coords Matrix of coordinates with one row per site, without missing or infinite
#'   values, for `null = "variogram"` or `"parametric"`. Taken from an [mc_data()] object when
#'   it has them.
#' @return List of class `mc_pcps` with `tests`, a data frame with one row for all the
#'   dissimilarity (`whole`) and one for each axis, and the columns `share`, the share of the
#'   positive eigenvalues; `r2`; `statistic`, the F statistic; `p_site`, `p_species` and
#'   `p_max`; `scores`, the sites on all the axes with positive eigenvalues; `share`, the
#'   share of each of those axes; `species`, the
#'   correlation of each species' column of the matrix P with those axes (its position in the
#'   biplot of the package PCPS); `env`, the predictors of the occupied sites; and
#'   `permutations`, the numbers of permutations and whether those of sites were free or
#'   restricted; and `tree`, the tree pruned to the species present. [plot.mc_pcps()] draws
#'   the ordination.
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
#' @family the whole tree
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
pcps_test <- function(comm, tree, env, axes = 1:2, n_perm = 999, site_permutations = NULL,
                      null = c("permutation", "variogram", "parametric"), coords = NULL) {
  cite_method("pcps_test")
  null <- check_null(null, site_permutations)
  if (inherits(comm, "mc_data")) {
    if (missing(tree)) tree <- data_part(comm, "tree")
    if (missing(env)) env <- data_part(comm, "env")
    if (is.null(coords)) coords <- comm$coords
    comm <- comm$comm
  }
  check_count(n_perm, "n_perm")
  check_pcps_axes(axes)
  if (!is.null(coords)) check_coords(coords)
  matched <- present_species(check_phylo_comm(comm, tree, !is.null(site_permutations)),
                             minimum = 3)
  comm <- matched$comm
  tree <- matched$tree
  landscape_env <- check_env(env, comm)
  check_permutations(site_permutations, comm)
  draws <- site_draws(comm, n_perm, site_permutations)
  comm <- comm[draws$occupied, , drop = FALSE]
  env <- check_env(landscape_env[draws$occupied, , drop = FALSE], comm)
  check_enough_sites(comm, env)

  weights <- comm / rowSums(comm)
  belonging <- fuzzy_belonging(tree)
  gower <- pcps_gower(weights, belonging)
  ordination <- pcps_axes(gower)
  scores <- pcps_chosen_axes(ordination$scores, axes)
  observed <- pcps_fit(gower, scores, env)
  warn_perfect_fit(observed$r2, c("the whole dissimilarity", colnames(scores)),
                   "the predictors and the response")
  site <- site_null(null, draws, landscape_env, env, coords, n_perm)
  null_site <- progress_map(seq_len(site$n), function(r) {
    pcps_fit(gower, scores, site$env_at(r))$statistic
  }, "Sites")
  null_species <- pcps_species_null(weights, belonging, axes, env, n_perm)
  tests <- pcps_tests(observed, null_site, null_species, axes, ordination$share)
  species <- pcps_species_cor(weights, belonging, ordination$scores)
  rownames(ordination$scores) <- rownames(comm)
  structure(list(tests = tests, scores = ordination$scores, share = ordination$share,
                 species = species, env = env, tree = tree,
                 permutations = permutation_record(n_perm, site_permutations, null = null)),
            class = "mc_pcps")
}

#' stops unless `axes` of `pcps_test()` are distinct whole numbers, 1 or more
#' @noRd
check_pcps_axes <- function(axes) {
  if (not_finite_numbers(axes) || any(axes < 1) || any(axes != round(axes)) ||
        anyDuplicated(axes) > 0) {
    stop("`axes` must be distinct whole numbers, 1 or more", call. = FALSE)
  }
}

#' TRUE unless `x` holds one or more finite numbers
#' @noRd
not_finite_numbers <- function(x) {
  !is.numeric(x) || length(x) < 1 || !all(is.finite(x))
}

#' the scores of the `axes` asked of `pcps_test()`, which must be among those of the ordination
#' @noRd
pcps_chosen_axes <- function(scores, axes) {
  if (max(axes) > ncol(scores)) {
    stop("`axes` asks for axis ", max(axes), " and there are ", ncol(scores),
         " with positive eigenvalues", call. = FALSE)
  }
  scores[, axes, drop = FALSE]
}

#' correlations between the columns of the matrix P, one per species, and the axes
#' @noRd
pcps_species_cor <- function(weights, belonging, scores) {
  species <- suppressWarnings(stats::cor(weights %*% belonging, scores))
  rownames(species) <- colnames(weights)
  species
}

#' statistics of `pcps_test()` with the species permuted on the tree, `n_perm` times: a list of
#' vectors with the statistic of the whole dissimilarity and of each of `axes`, NA for axes
#' that a permutation does not have
#' @noRd
pcps_species_null <- function(weights, belonging, axes, env, n_perm) {
  progress_map(seq_len(n_perm), function(r) {
    moved <- sample.int(ncol(weights))
    null_gower <- pcps_gower(weights, belonging[moved, moved])
    null_scores <- pcps_axes(null_gower)$scores
    present <- axes <= ncol(null_scores)
    statistic <- rep(NA_real_, length(axes) + 1)
    statistic[c(TRUE, present)] <- pcps_fit(null_gower, null_scores[, axes[present], drop = FALSE],
                                            env)$statistic
    statistic
  }, "Species")
}

#' the table of tests of `pcps_test()`, with the p-values of the observed statistics against
#' the draws of sites and of species
#' @noRd
pcps_tests <- function(observed, null_site, null_species, axes, share) {
  p_value <- function(nulls) {
    vapply(seq_along(observed$statistic), function(i) {
      permutation_p(vapply(nulls, `[`, numeric(1), i), observed$statistic[i])
    }, numeric(1))
  }
  p_site <- p_value(null_site)
  p_species <- p_value(null_species)
  data.frame(test = c("whole", paste0("pcps", axes)),
             share = c(1, share[axes]), r2 = observed$r2,
             statistic = observed$statistic, p_site = p_site, p_species = p_species,
             p_max = pmax(p_site, p_species))
}

#' @export
print.mc_pcps <- function(x, digits = 3, ...) {
  check_count(digits, "digits")
  cat("<mc_pcps> principal coordinates of phylogenetic structure\n")
  permutation_line(x$permutations)
  print(format(x$tests, digits = digits), row.names = FALSE, quote = FALSE)
  invisible(x)
}

#' degree of belonging of each species to the others, from the patristic distances: rows
#' sum to 1 (Pillar & Duarte 2010)
#' @noRd
fuzzy_belonging <- function(tree) {
  distances <- ape::cophenetic.phylo(tree)[tree$tip.label, tree$tip.label]
  similarity <- 1 - distances / max(distances)
  similarity / rowSums(similarity)
}

#' Gower-centred matrix of the square root of Bray-Curtis among rows of the matrix P
#' @noRd
pcps_gower <- function(weights, belonging) {
  p <- weights %*% belonging
  sums <- rowSums(p)
  half <- -0.5 * as.matrix(stats::dist(p, method = "manhattan")) / outer(sums, sums, "+")
  half <- sweep(half, 1, rowMeans(half))
  unname(sweep(half, 2, colMeans(half)))
}

#' axes with positive eigenvalues, scaled by the root of the eigenvalue, and their shares
#' @noRd
pcps_axes <- function(gower) {
  eigenvalues <- eigen(gower, symmetric = TRUE)
  keep <- eigenvalues$values > sqrt(.Machine$double.eps) * max(eigenvalues$values)
  values <- eigenvalues$values[keep]
  scores <- sweep(eigenvalues$vectors[, keep, drop = FALSE], 2, sqrt(values), "*")
  colnames(scores) <- paste0("pcps", seq_along(values))
  list(scores = scores, share = values / sum(values))
}

#' R2 and F of the distance-based regression of all the dissimilarity, then of each axis;
#' NA when a permutation leaves the predictors without variance
#' @noRd
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
#' as in the biplot of the package PCPS (Debastiani & Duarte 2014). With `clades`, the species
#' are drawn as stars coloured by clade, each clade with its outline and name, as in Fig. 1 of
#' Debastiani & Duarte (2014), which shows which lineages go with which sites.
#'
#' @param x Result of [pcps_test()].
#' @param axes The two principal coordinates to draw: two different whole numbers.
#' @param env Predictor that colours the sites: a column number or name of the predictors given
#'   to [pcps_test()], or `NULL` for no colour.
#' @param species `TRUE` or `FALSE`. Draw the species as arrows.
#' @param species_labels `TRUE` or `FALSE`. Write the names of the species at the arrows or stars.
#' @param clades Optional grouping of the species into clades: one whole number k, which cuts
#'   the tree into its k deepest clades by splitting, one at a time, the clade whose crown is
#'   nearest the root; or a vector of clade names named by species, such as families, in which
#'   species left out are drawn in grey without a clade. Each clade with at least three species
#'   gets its convex hull.
#' @param palette One character string, the name of a palette of [grDevices::hcl.colors()] for
#'   the predictor.
#' @param point_size Size of the site points. One number greater than 0.
#' @param ... Not used.
#' @return A ggplot object.
#' @seealso [pcps_test()].
#' @references
#' Debastiani, V. J., & Duarte, L. D. S. (2014). PCPS -- an R-package for exploring
#' phylogenetic eigenvectors across metacommunities. *Frontiers of Biogeography*, 6.
#' \doi{10.21425/F5FBG22943}
#' @family the whole tree
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(12)
#' traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100), rmax_range = c(0, 30))
#' env <- seq(0, 100, length.out = 36)
#' comm <- simulate_static(traits, env, incidence = TRUE)$abundance
#' plot(pcps_test(comm, tree, env, n_perm = 19))
#' @export
plot.mc_pcps <- function(x, axes = c(1, 2), env = 1, species = TRUE, species_labels = TRUE,
                         clades = NULL, palette = "viridis", point_size = 2.5, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot() of pcps_test() needs the package ggplot2", call. = FALSE) # nocov
  }
  if (is.null(x$species)) {
    stop("`x` comes from an older version of pcps_test(); run it again", call. = FALSE)
  }
  check_plot_axes(axes, ncol(x$scores))
  check_flag(species, "species")
  check_flag(species_labels, "species_labels")
  check_string(palette, "palette")
  check_number(point_size, "point_size", lower = 0, open = TRUE)
  axis_label <- function(a) sprintf("PCPS %d (%.1f%%)", a, 100 * x$share[a])
  sites <- data.frame(a1 = x$scores[, axes[1]], a2 = x$scores[, axes[2]])
  colour_name <- NULL
  if (!is.null(env)) {
    predictor <- pcps_predictor(x, env)
    colour_name <- predictor$name
    sites$env <- predictor$values
  }
  .data <- ggplot2::.data
  p <- ggplot2::ggplot(sites, ggplot2::aes(x = .data$a1, y = .data$a2)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85")
  group <- NULL
  if (!is.null(clades)) {
    if (is.null(x$tree)) {
      stop("`x` comes from an older version of pcps_test(); run it again to use `clades`",
           call. = FALSE)
    }
    group <- pcps_clades(clades, x$tree, rownames(x$species))
  }
  p <- p + pcps_site_layers(colour_name, !is.null(group) && species, palette, point_size)
  if (species) p <- p + pcps_species_layers(x, axes, sites, group, species_labels)
  p + ggplot2::labs(x = axis_label(axes[1]), y = axis_label(axes[2])) +
    ggplot2::coord_equal() + ggplot2::theme_minimal()
}

#' stops unless `axes` are two different axes among the `n_axes` of `x`
#' @noRd
check_plot_axes <- function(axes, n_axes) {
  if (not_pair_within(axes, n_axes) || any(axes != round(axes)) || axes[1] == axes[2]) {
    stop("`axes` must be two different axes between 1 and ", n_axes, call. = FALSE)
  }
}

#' TRUE unless `x` holds two numbers between 1 and `n`
#' @noRd
not_pair_within <- function(x, n) {
  !is.numeric(x) || length(x) != 2 || any(x < 1) || any(x > n)
}

#' the predictor that colours the sites, from its column number or name: a list with its `name`
#' and `values`
#' @noRd
pcps_predictor <- function(x, env) {
  predictors <- as.data.frame(x$env)
  if (is.null(colnames(x$env))) names(predictors) <- paste0("env", seq_along(predictors))
  if (length(env) != 1 || !(env %in% names(predictors) ||
                            (is.numeric(env) && env >= 1 && env <= ncol(predictors)))) {
    stop("`env` must be one column of the predictors given to pcps_test(), or NULL",
         call. = FALSE)
  }
  list(name = if (is.numeric(env)) names(predictors)[env] else env, values = predictors[[env]])
}

#' the points of the sites: plain without a predictor, and otherwise coloured by it, or filled
#' with `by_fill`, when the colour scale belongs to the clades
#' @noRd
pcps_site_layers <- function(colour_name, by_fill, palette, point_size) {
  .data <- ggplot2::.data
  if (is.null(colour_name)) return(list(ggplot2::geom_point(size = point_size)))
  if (!by_fill) {
    return(list(ggplot2::geom_point(ggplot2::aes(colour = .data$env), size = point_size),
                ggplot2::scale_colour_gradientn(colours = grDevices::hcl.colors(30, palette),
                                                name = colour_name)))
  }
  # the colour scale belongs to the clades, so the predictor fills the sites
  list(ggplot2::geom_point(ggplot2::aes(fill = .data$env), shape = 21, colour = "grey30",
                           stroke = 0.3, size = point_size),
       ggplot2::scale_fill_gradientn(colours = grDevices::hcl.colors(30, palette),
                                     name = colour_name))
}

#' the species as arrows, or as stars by clade when `group` gives the clade of each, with
#' their names when `species_labels` is TRUE; a list of layers
#' @noRd
pcps_species_layers <- function(x, axes, sites, group, species_labels) {
  .data <- ggplot2::.data
  # species scaled to the spread of the sites, as in a biplot
  scale <- 0.9 * max(abs(as.matrix(sites[c("a1", "a2")]))) /
    max(abs(x$species[, axes]), na.rm = TRUE)
  sp <- data.frame(a1 = x$species[, axes[1]] * scale, a2 = x$species[, axes[2]] * scale,
                   label = rownames(x$species))
  if (!is.null(group)) sp$clade <- group
  sp <- sp[stats::complete.cases(sp[c("a1", "a2")]), , drop = FALSE]
  layers <- if (is.null(group)) {
    list(ggplot2::geom_segment(data = sp, ggplot2::aes(x = 0, y = 0, xend = .data$a1,
                                                       yend = .data$a2),
                               colour = "grey50", linewidth = 0.4,
                               arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"))))
  } else {
    pcps_clade_layers(sp)
  }
  if (species_labels) {
    layers <- c(layers, list(ggplot2::geom_text(data = sp, ggplot2::aes(label = .data$label),
                                                colour = "grey30", size = 3, vjust = -0.4)))
  }
  layers
}

#' clade of each species, from a number of clades or a vector named by species
#' @noRd
pcps_clades <- function(clades, tree, species) {
  if (is.numeric(clades) && length(clades) == 1 && is.null(names(clades))) {
    return(deepest_clades(clades, tree, species))
  }
  if (is.null(names(clades)) || length(clades) == 0) {
    stop("`clades` must be a number of clades or a vector of clade names named by species",
         call. = FALSE)
  }
  if (!any(species %in% names(clades))) {
    stop("no species of `x` is among the names of `clades`", call. = FALSE)
  }
  factor(as.character(clades)[match(species, names(clades))])
}

#' the `k` deepest clades of the tree, named by letters from the oldest crown, as a factor with
#' the clade of each of `species`
#' @noRd
deepest_clades <- function(k, tree, species) {
  n_tip <- ape::Ntip(tree)
  if (is.na(k) || k != round(k) || k < 2 || k > n_tip) {
    stop("`clades` as a number must be a whole number between 2 and ", n_tip, call. = FALSE)
  }
  depth <- ape::node.depth.edgelength(tree)
  groups <- split_deepest(tree, depth, k)
  member <- character(n_tip)
  for (g in seq_along(groups)) {
    tips <- if (groups[g] <= n_tip) tree$tip.label[groups[g]] else
      ape::extract.clade(tree, groups[g])$tip.label
    member[match(tips, tree$tip.label)] <- LETTERS[(g - 1) %% 26 + 1]
  }
  factor(member[match(species, tree$tip.label)], levels = LETTERS[seq_along(groups)])
}

#' the crown nodes (or tips) of the `k` clades left after splitting, one at a time, the clade
#' whose crown is nearest the root; ordered by the depth of the crown
#' @noRd
split_deepest <- function(tree, depth, k) {
  n_tip <- ape::Ntip(tree)
  groups <- n_tip + 1
  while (length(groups) < k) {
    internal <- groups[groups > n_tip]
    split <- internal[which.min(depth[internal])]
    groups <- c(setdiff(groups, split), tree$edge[tree$edge[, 1] == split, 2])
  }
  groups[order(depth[groups])]
}

#' stars coloured by clade, convex hulls of clades with three or more species and their names
#' @noRd
pcps_clade_layers <- function(sp) {
  .data <- ggplot2::.data
  grouped <- sp[!is.na(sp$clade), , drop = FALSE]
  hulls <- do.call(rbind, lapply(split(grouped, grouped$clade, drop = TRUE), function(d) {
    if (nrow(d) < 3) return(NULL)
    d[grDevices::chull(d$a1, d$a2), c("a1", "a2", "clade")]
  }))
  centres <- stats::aggregate(cbind(a1, a2) ~ clade, grouped, mean)
  layers <- list()
  if (!is.null(hulls)) {
    layers <- c(layers, ggplot2::geom_polygon(data = hulls,
                                              ggplot2::aes(colour = .data$clade,
                                                           group = .data$clade),
                                              fill = NA, linetype = 2, linewidth = 0.4,
                                              show.legend = FALSE))
  }
  if (anyNA(sp$clade)) {
    layers <- c(layers, ggplot2::geom_point(data = sp[is.na(sp$clade), , drop = FALSE],
                                            shape = 8, colour = "grey60", size = 2))
  }
  n_clades <- nlevels(droplevels(grouped$clade))
  # names under the stars, so a clade of one species still shows its star
  c(layers,
    ggplot2::geom_label(data = centres, ggplot2::aes(label = .data$clade,
                                                     colour = .data$clade),
                        size = 3.5, fontface = "bold", fill = grDevices::adjustcolor("white", 0.7),
                        vjust = 1.4, show.legend = FALSE),
    ggplot2::geom_point(data = grouped, ggplot2::aes(colour = .data$clade), shape = 8, size = 2),
    ggplot2::scale_colour_manual(values = grDevices::hcl.colors(max(n_clades, 2), "Dark 3"),
                                 name = "Clade"))
}
