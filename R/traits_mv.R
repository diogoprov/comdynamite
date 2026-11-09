#' Evolutionary rate matrices with restricted structures
#'
#' Builds the evolutionary rate matrix (the covariance of the changes of several traits per
#' unit time) for multivariate trait evolution with [traits_ou()], with the structures that
#' mvMORPH fits (`mvBM(param = list(constraint = ))`, Clavel et al. 2015). Variances are rates
#' of evolution of each trait, and correlations are evolutionary correlations between traits.
#'
#' One matrix, for the whole tree:
#'
#' * `"default"`: any rate matrix, with the variances in `variance` and the correlations in
#'   `correlation`. A matrix supplied by the user can also be passed directly to
#'   `traits_ou(sigma2 = )`.
#' * `"diagonal"`: traits evolve independently, each with its own rate.
#' * `"equaldiagonal"`: traits evolve independently, all with the same rate.
#' * `"equal"`: all traits with the same rate, correlated (Adams 2013).
#'
#' One matrix per selective regime, for multi-regime models (`regimes` in [traits_ou()]):
#'
#' * `"proportional"`: matrices proportional to each other, the matrix built from `variance`
#'   and `correlation` multiplied by `scale`, one positive number per regime.
#' * `"shared"`: matrices with the same eigenvectors, those of the matrix built from
#'   `variance` and `correlation`, and the eigenvalues of each regime in `eigenvalues`, in the
#'   order of [eigen()] (decreasing eigenvalues of that matrix).
#' * `"variance"`: the same variances in all regimes, and the correlations of each regime in
#'   `correlation`, a named list.
#' * `"correlation"`: the same correlations in all regimes, and the variances of each regime
#'   in `variance`, a named list.
#'
#' The regimes are named by the names of `scale`, `eigenvalues`, or the list in `correlation`
#' or `variance`, and those names must be the regimes of the tree.
#'
#' @param variance Rate of each trait (variance per unit time), recycled to the number of
#'   traits. For `"correlation"`, a named list with the rates of each regime.
#' @param correlation Correlation between traits: one number for all pairs or a correlation
#'   matrix. For `"variance"`, a named list with the correlations of each regime.
#' @param constraint Structure of the matrix or matrices; see Details.
#' @param n_traits Number of traits, when neither `variance` nor `correlation` gives it.
#' @param scale For `"proportional"`, the multiplier of each regime, named by regime.
#' @param eigenvalues For `"shared"`, a named list with the eigenvalues of each regime, one per
#'   trait.
#' @return An object of class `rate_matrix`: a symmetric matrix with one row and column per
#'   trait, or, for the multi-regime structures, a named list of them, one per regime. The
#'   attribute `constraint` holds the structure. Either can be passed to
#'   `traits_ou(sigma2 = )`; [plot()][plot.rate_matrix] draws them as ellipses.
#' @references
#'   Adams, D. C. (2013). Comparing evolutionary rates for different phenotypic traits on a
#'   phylogeny using likelihood. *Systematic Biology*, 62, 181--192.
#'   \doi{10.1093/sysbio/sys083}
#'
#'   Clavel, J., Escarguel, G., & Merceron, G. (2015). mvMORPH: an R package for fitting
#'   multivariate evolutionary models to morphometric data. *Methods in Ecology and
#'   Evolution*, 6, 1311--1319. \doi{10.1111/2041-210X.12420}
#' @seealso [traits_ou()], which simulates traits with these matrices, and [plot_traits()].
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' rate_matrix(c(1, 2), correlation = 0.5)
#' rate_matrix(1, correlation = 0.3, constraint = "equal", n_traits = 3)
#' # two habitats, the open one evolving four times faster in all directions
#' r <- rate_matrix(c(1, 2), correlation = 0.5, constraint = "proportional",
#'                  scale = c(forest = 1, open = 4))
#' r
#' plot(r)
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(30)
#' habitat <- trait_mk(tr, q = 0.5, states = c("forest", "open"))
#' x <- traits_ou(habitat, sigma2 = r, regimes = habitat)
#' plot_traits(x, habitat, regimes = habitat)
#' # several correlated traits in community assembly
#' assemble_community(x, size = 10, rule = "limiting_similarity")$community
#' @export
rate_matrix <- function(variance = 1, correlation = 0,
                        constraint = c("default", "diagonal", "equaldiagonal", "equal",
                                       "proportional", "shared", "variance", "correlation"),
                        n_traits = NULL, scale = NULL, eigenvalues = NULL) {
  constraint <- match.arg(constraint)
  first <- function(x) if (is.list(x)) x[[1]] else x
  p <- n_traits %||% {
    rho <- first(correlation)
    if (is.matrix(rho)) nrow(rho) else length(first(variance))
  }
  if (length(p) != 1 || !is.finite(p) || p < 1) {
    stop("`n_traits` must be one positive number", call. = FALSE)
  }
  p <- as.integer(p)
  if (p == 1 && is.null(n_traits) && constraint != "default") {
    stop("give `n_traits`, or one variance per trait", call. = FALSE)
  }
  regime_list <- function(x, what) {
    if (!is.list(x) || is.null(names(x)) || any(names(x) == "") || anyDuplicated(names(x))) {
      stop("for \"", constraint, "\", `", what, "` must be a list named by regime",
           call. = FALSE)
    }
    x
  }
  if (is.list(variance) && constraint != "correlation") {
    stop("`variance` can be a list only with constraint \"correlation\"", call. = FALSE)
  }
  if (is.list(correlation) && constraint != "variance") {
    stop("`correlation` can be a list only with constraint \"variance\"", call. = FALSE)
  }
  out <- switch(constraint,
    default = build_rate_matrix(variance, correlation, p),
    diagonal = {
      if (!all(correlation == 0)) stop("\"diagonal\" has no correlations", call. = FALSE)
      build_rate_matrix(variance, 0, p)
    },
    equaldiagonal = {
      if (!all(correlation == 0)) stop("\"equaldiagonal\" has no correlations", call. = FALSE)
      if (length(variance) != 1) stop("\"equaldiagonal\" has one variance", call. = FALSE)
      build_rate_matrix(variance, 0, p)
    },
    equal = {
      if (length(variance) != 1) stop("\"equal\" has one variance", call. = FALSE)
      build_rate_matrix(variance, correlation, p)
    },
    proportional = {
      if (!is.numeric(scale) || is.null(names(scale)) || anyNA(scale) || any(scale <= 0)) {
        stop("for \"proportional\", `scale` must be positive numbers named by regime",
             call. = FALSE)
      }
      base <- build_rate_matrix(variance, correlation, p)
      lapply(as.list(scale), function(s) base * s)
    },
    shared = {
      eigenvalues <- regime_list(eigenvalues, "eigenvalues")
      q <- eigen(build_rate_matrix(variance, correlation, p), symmetric = TRUE)$vectors
      lapply(eigenvalues, function(l) {
        if (!is.numeric(l) || length(l) != p || anyNA(l) || any(l < 0)) {
          stop("`eigenvalues` needs ", p, " values, zero or more, per regime", call. = FALSE)
        }
        m <- q %*% diag(l, p) %*% t(q)
        (m + t(m)) / 2
      })
    },
    variance = {
      lapply(regime_list(correlation, "correlation"),
             function(r) build_rate_matrix(variance, r, p))
    },
    correlation = {
      lapply(regime_list(variance, "variance"),
             function(v) build_rate_matrix(v, correlation, p))
    })
  names_traits <- paste0("trait", seq_len(p))
  name_it <- function(m) {
    dimnames(m) <- list(names_traits, names_traits)
    m
  }
  out <- if (is.list(out)) lapply(out, name_it) else name_it(out)
  attr(out, "constraint") <- constraint
  class(out) <- c("rate_matrix", class(out))
  out
}

#' @rdname rate_matrix
#' @param x A `rate_matrix` object.
#' @param ... Not used.
#' @export
print.rate_matrix <- function(x, ...) {
  several <- is.list(x)
  m <- if (several) x[[1]] else x
  cat("<rate_matrix>", attr(x, "constraint"), "structure,", ncol(m), "traits")
  if (several) cat(",", length(x), "regimes")
  cat("\n")
  if (several) {
    for (k in names(x)) {
      cat("\nregime ", k, ":\n", sep = "")
      print(unclass(x[[k]]))
    }
  } else {
    print(matrix(unclass(x), nrow(x), dimnames = dimnames(x)))
  }
  invisible(x)
}

# the rate matrices of sigma2 as plain matrices
strip_rate_matrix <- function(x) {
  if (!inherits(x, "rate_matrix")) return(x)
  if (is.list(x)) {
    lapply(unclass(x), unclass)
  } else {
    matrix(unclass(x), nrow(x), dimnames = dimnames(x))
  }
}

# rate matrix from the variances and a common correlation or a correlation matrix
build_rate_matrix <- function(variance, correlation, p) {
  if (!is.numeric(variance) || anyNA(variance) || any(variance < 0) ||
      !length(variance) %in% c(1, p)) {
    stop("`variance` needs 1 or ", p, " values, zero or more", call. = FALSE)
  }
  if (is.matrix(correlation)) {
    rho <- unname(correlation)
    if (any(dim(rho) != p) || anyNA(rho) || !isSymmetric(rho) || any(diag(rho) != 1) ||
        any(abs(rho) > 1)) {
      stop("`correlation` must be a symmetric ", p, " x ", p, " matrix with ones on the ",
           "diagonal and values in [-1, 1]", call. = FALSE)
    }
  } else {
    if (!is.numeric(correlation) || length(correlation) != 1 || is.na(correlation) ||
        abs(correlation) > 1) {
      stop("`correlation` must be one number in [-1, 1] or a correlation matrix",
           call. = FALSE)
    }
    rho <- matrix(correlation, p, p)
    diag(rho) <- 1
  }
  if (min(eigen(rho, symmetric = TRUE, only.values = TRUE)$values) < -1e-10) {
    stop("the correlations do not form a valid correlation matrix (not positive ",
         "semidefinite)", call. = FALSE)
  }
  s <- sqrt(rep_len(variance, p))
  rho * outer(s, s)
}
