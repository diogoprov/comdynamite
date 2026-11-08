# Classic metrics of the spacing and spread of trait values within a community, written from
# their definitions. Each returns one number, so that it can be given to community_index().

#' Classic metrics of trait spacing and spread
#'
#' Metrics of the distribution of trait values within a community used in tests of community
#' assembly: nearest-neighbour distances (Weiher et al. 1998; Stubbs & Wilson 2004), all
#' pairwise distances, the minimum spanning tree, and the range and kurtosis (Cornwell &
#' Ackerly 2009). Aiba et al. (2013) compared the statistical power of many of these metrics.
#' Each function returns one number, so that it can be given as `index` to
#' [community_index()]; with `multivariate = TRUE` there, the metrics of distances receive the
#' matrix of traits of each site.
#'
#' * `nn_distances()`: distance from each value to its nearest neighbour (a vector).
#' * `mnnd()`, `sdnnd()`, `cvnnd()`, `min_nnd()`, `max_nnd()`: mean, standard deviation,
#'   coefficient of variation (sd / mean), minimum and maximum of the nearest-neighbour
#'   distances.
#' * `mnd()`, `sdnd()`, `cvnd()`: mean, standard deviation and coefficient of variation of all
#'   pairwise distances.
#' * `mst_length()`: total length of the minimum spanning tree; `mst_ratio()`: shortest over
#'   longest edge of the tree, as in `MinMaxMST()` of the package cati.
#' * `trait_range()`: range of a single trait; `trait_kurtosis()`: excess kurtosis of a single
#'   trait (fourth standardised moment, with population moments, minus 3).
#'
#' Distances are Euclidean among the rows of a matrix, or the absolute differences of a vector;
#' a `dist` object (for instance Gower distances from `cluster::daisy()`) is used as it is.
#' Identical values are at distance zero and count as such: an individual whose nearest
#' neighbour has the same value has a nearest-neighbour distance of zero (cati discards zero
#' distances). For a single trait the nearest neighbours and the minimum spanning tree are found
#' from the sorted values, without the matrix of distances, so that large communities are fast.
#'
#' Traits on different scales should be standardised over all communities before the metrics
#' are computed (for instance with [scale()] on the whole data). Standardising within each
#' community, as cati does by default, changes the units from one community to the next and
#' between observed and null communities.
#'
#' @param x Numeric vector (one trait), matrix or data frame (traits in columns, individuals or
#'   species in rows) or a `dist` object.
#' @param relative Divide by the range of the trait (one trait) or by the largest pairwise
#'   distance (several traits), as the `div_range` option of cati.
#' @param na.rm Drop missing values (rows with any missing value of a matrix) before the
#'   computation; with `FALSE` (the default) missing values give `NA`.
#' @return A single number, or `NA` when there are too few values: two for the metrics of
#'   distances and the range, three for the standard deviation of nearest-neighbour distances,
#'   four for the kurtosis. `nn_distances()` returns a vector with one distance per value.
#' @references Aiba, M., Katabuchi, M., Takafumi, H., Matsuzaki, S. S., Sasaki, T., & Hiura, T.
#'   (2013). Robustness of trait distribution metrics for community assembly studies under the
#'   uncertainties of assembly processes. *Ecology*, 94, 2873--2885. \doi{10.1890/13-0269.1}
#'
#'   Cornwell, W. K., & Ackerly, D. D. (2009). Community assembly and shifts in plant trait
#'   distributions across an environmental gradient in coastal California. *Ecological
#'   Monographs*, 79, 109--126. \doi{10.1890/07-1134.1}
#'
#'   Stubbs, W. J., & Wilson, J. B. (2004). Evidence for limiting similarity in a sand dune
#'   community. *Journal of Ecology*, 92, 557--567. \doi{10.1111/j.0022-0477.2004.00898.x}
#'
#'   Weiher, E., Clarke, G. D. P., & Keddy, P. A. (1998). Community assembly rules,
#'   morphological dispersion, and the coexistence of plant species. *Oikos*, 81, 309--322.
#'   \doi{10.2307/3547051}
#' @seealso [community_index()].
#' @examples
#' x <- c(1, 2, 2, 5, 9)
#' nn_distances(x)
#' mnnd(x)
#' mst_length(x)
#' mnnd(cbind(x, rev(x)))
#' # as an index of communities against a null model
#' set.seed(1)
#' ind <- data.frame(site = rep(paste0("s", 1:5), each = 15),
#'                   species = sample(letters[1:4], 75, replace = TRUE), trait = rnorm(75))
#' community_index(ind, mnnd, null = "regional.ind", nperm = 99)
#' @name neighbour_metrics
NULL

#' @rdname neighbour_metrics
#' @export
nn_distances <- function(x, na.rm = FALSE) {
  x <- metric_input(x, na.rm)
  if (is.null(x)) return(NA_real_)
  if (is.numeric(x) && is.null(dim(x)) && !inherits(x, "dist")) {
    n <- length(x)
    if (n < 2) return(rep(NA_real_, n))
    o <- order(x)
    gaps <- diff(x[o])
    nearest <- pmin(c(Inf, gaps), c(gaps, Inf))
    out <- numeric(n)
    out[o] <- nearest
    return(out)
  }
  d <- as.matrix(metric_dist(x))
  if (nrow(d) < 2) return(rep(NA_real_, nrow(d)))
  diag(d) <- Inf
  apply(d, 1, min)
}

#' @rdname neighbour_metrics
#' @export
mnnd <- function(x, relative = FALSE, na.rm = FALSE) {
  nn_summary(x, mean, 2, relative, na.rm)
}

#' @rdname neighbour_metrics
#' @export
sdnnd <- function(x, relative = FALSE, na.rm = FALSE) {
  nn_summary(x, stats::sd, 3, relative, na.rm)
}

#' @rdname neighbour_metrics
#' @export
cvnnd <- function(x, na.rm = FALSE) {
  nn_summary(x, function(d) stats::sd(d) / mean(d), 3, FALSE, na.rm)
}

#' @rdname neighbour_metrics
#' @export
min_nnd <- function(x, relative = FALSE, na.rm = FALSE) {
  nn_summary(x, min, 2, relative, na.rm)
}

#' @rdname neighbour_metrics
#' @export
max_nnd <- function(x, relative = FALSE, na.rm = FALSE) {
  nn_summary(x, max, 2, relative, na.rm)
}

#' @rdname neighbour_metrics
#' @export
mnd <- function(x, relative = FALSE, na.rm = FALSE) {
  pair_summary(x, mean, relative, na.rm)
}

#' @rdname neighbour_metrics
#' @export
sdnd <- function(x, relative = FALSE, na.rm = FALSE) {
  pair_summary(x, stats::sd, relative, na.rm)
}

#' @rdname neighbour_metrics
#' @export
cvnd <- function(x, na.rm = FALSE) {
  pair_summary(x, function(d) stats::sd(d) / mean(d), FALSE, na.rm)
}

#' @rdname neighbour_metrics
#' @export
mst_length <- function(x, relative = FALSE, na.rm = FALSE) {
  edges <- mst_edges(x, na.rm)
  if (length(edges) == 0 || anyNA(edges)) return(NA_real_)
  sum(edges) / if (relative) metric_scale(x, na.rm) else 1
}

#' @rdname neighbour_metrics
#' @export
mst_ratio <- function(x, na.rm = FALSE) {
  edges <- mst_edges(x, na.rm)
  if (length(edges) == 0 || anyNA(edges)) return(NA_real_)
  if (max(edges) == 0) return(NA_real_)
  min(edges) / max(edges)
}

#' @rdname neighbour_metrics
#' @export
trait_range <- function(x, na.rm = FALSE) {
  x <- single_trait(x, na.rm)
  if (anyNA(x) || length(x) < 2) return(NA_real_)
  diff(range(x))
}

#' @rdname neighbour_metrics
#' @export
trait_kurtosis <- function(x, na.rm = FALSE) {
  x <- single_trait(x, na.rm)
  if (anyNA(x) || length(x) < 4) return(NA_real_)
  centred <- x - mean(x)
  m2 <- mean(centred^2)
  if (m2 == 0) return(NA_real_)
  mean(centred^4) / m2^2 - 3
}

# checks the input of the metrics and drops missing values when asked; returns NULL when a
# missing value remains (the metric is then NA)
metric_input <- function(x, na.rm) {
  if (inherits(x, "dist")) {
    if (anyNA(x)) {
      if (na.rm) stop("a `dist` object with missing values cannot be used", call. = FALSE)
      return(NULL)
    }
    return(x)
  }
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.numeric(x)) stop("`x` must be numeric, a numeric matrix or a dist object", call. = FALSE)
  if (is.matrix(x) && ncol(x) == 1) x <- as.vector(x)
  if (is.matrix(x)) {
    complete <- stats::complete.cases(x)
    if (!all(complete)) {
      if (!na.rm) return(NULL)
      x <- x[complete, , drop = FALSE]
    }
  } else if (anyNA(x)) {
    if (!na.rm) return(NULL)
    x <- x[!is.na(x)]
  }
  x
}

metric_dist <- function(x) if (inherits(x, "dist")) x else stats::dist(x)

n_values <- function(x) {
  if (inherits(x, "dist")) attr(x, "Size") else if (is.matrix(x)) nrow(x) else length(x)
}

# range of one trait, or the largest pairwise distance
metric_scale <- function(x, na.rm) {
  x <- metric_input(x, na.rm)
  if (is.numeric(x) && is.null(dim(x)) && !inherits(x, "dist")) diff(range(x)) else
    max(metric_dist(x))
}

nn_summary <- function(x, f, min_n, relative, na.rm) {
  clean <- metric_input(x, na.rm)
  if (is.null(clean) || n_values(clean) < min_n) return(NA_real_)
  value <- f(nn_distances(clean))
  if (relative) value <- value / metric_scale(clean, FALSE)
  finite_or_na(value)
}

pair_summary <- function(x, f, relative, na.rm) {
  clean <- metric_input(x, na.rm)
  if (is.null(clean) || n_values(clean) < 2) return(NA_real_)
  d <- as.vector(metric_dist(clean))
  value <- f(d)
  if (relative) value <- value / max(d)
  finite_or_na(value)
}

# edge lengths of the minimum spanning tree; for one trait, the gaps between sorted values
mst_edges <- function(x, na.rm) {
  clean <- metric_input(x, na.rm)
  if (is.null(clean)) return(NA_real_)
  if (n_values(clean) < 2) return(numeric(0))
  if (is.numeric(clean) && is.null(dim(clean)) && !inherits(clean, "dist")) {
    return(diff(sort(clean)))
  }
  d <- as.matrix(metric_dist(clean))
  tree <- ape::mst(d)
  d[upper.tri(d) & tree == 1]
}

single_trait <- function(x, na.rm) {
  if (inherits(x, "dist") || (is.matrix(x) && ncol(x) > 1) ||
      (is.data.frame(x) && ncol(x) > 1)) {
    stop("this metric takes a single trait", call. = FALSE)
  }
  x <- as.numeric(unlist(x))
  if (na.rm) x <- x[!is.na(x)]
  x
}
