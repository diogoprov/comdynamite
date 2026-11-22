#' Competition matrix
#'
#' Intraspecific coefficients are 1 and interspecific coefficients lie in
#' \[0, `inter_max`), all multiplied by `scale`. With the default `scale` of
#' 0.05, as in Thompson et al. (2020), a species alone at its optimum with
#' `r_max = 5` reaches an equilibrium of \eqn{(5 - 1) / 0.05 = 80} individuals,
#' and `inter_max` is the strength of interspecific relative to intraspecific
#' competition. The matrix is symmetric.
#'
#' * `"none"`: no interspecific competition.
#' * `"random"`: interspecific coefficients drawn from U(0, `inter_max`).
#' * `"phylo"`: the same uniform values, assigned to pairs of species by rank
#'   of phylogenetic similarity (1 - d / max(d), with d the patristic
#'   distance), so that close relatives compete more strongly while the
#'   distribution of coefficients stays the same as with `"random"`.
#' * `"trait"`: competition falls with the distance between the niches of the two species,
#'   \eqn{\alpha_{jk} = q \exp(-(\theta_j - \theta_k)^2 / (2(\sigma_j^2 +
#'   \sigma_k^2)))}{alpha_jk = q exp(-(theta_j - theta_k)^2 / (2 (sigma_j^2 + sigma_k^2)))},
#'   with \eqn{q}{q} = `inter_max`, \eqn{\theta}{theta} the optima and \eqn{\sigma}{sigma}
#'   the niche breadths of `traits`. The last factor is the overlap of two Gaussian niches, the
#'   kernel that [simulate_evolving_metacommunity()] uses with `niche_overlap = TRUE`, so that
#'   with `inter_max = 1` two species with the same niche compete as one species. With a list
#'   of traits (several axes), the overlap is the product over axes, or, for traits from
#'   [niche_correlation()], \eqn{\exp(-d^\top (\Omega_j + \Omega_k)^{-1} d / 2)}{exp(-d'
#'   (Omega_j + Omega_k)^-1 d / 2)}, with \eqn{d}{d} the difference between the optima and
#'   \eqn{\Omega}{Omega} the covariance of each niche. Nothing is random.
#'
#' @param species A `phylo` tree, or the number of species, one whole number, 1 or more (not for
#'   `type = "phylo"`); ignored for `type = "trait"`, where the species are those of `traits`.
#' @param type Type of interspecific competition.
#' @param traits For `type = "trait"`: a data frame from [species_traits()], or any data frame
#'   with columns `optimum` and `sigma`, or a list of them for a niche of several axes.
#' @param inter_max Maximum interspecific coefficient, relative to
#'   intraspecific competition, one finite number, 0 or more.
#' @param scale Multiplier of all coefficients, one finite number greater than 0.
#' @return Square matrix: element \[j, k\] is the effect of species k on
#'   species j.
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(5)
#' round(competition_matrix(tr, "phylo", inter_max = 0.5), 3)
#' traits <- species_traits(tr)
#' round(competition_matrix(type = "trait", traits = traits, inter_max = 1), 3)
#' @export
competition_matrix <- function(species, type = c("none", "random", "phylo", "trait"),
                               inter_max = 0.1, scale = 0.05, traits = NULL) {
  type <- match.arg(type)
  check_number(inter_max, "inter_max", lower = 0)
  check_number(scale, "scale", lower = 0, open = TRUE)
  if (type == "trait") return(scale * trait_competition(traits, inter_max))
  if (!is.null(traits)) stop("`traits` is used only with type = \"trait\"", call. = FALSE)
  is_tree <- inherits(species, "phylo")
  if (type == "phylo" && !is_tree) {
    stop("type = \"phylo\" needs a 'phylo' tree as `species`", call. = FALSE)
  }
  if (!is_tree) check_count(species, "species")
  sp <- if (is_tree) species$tip.label else paste0("sp", seq_len(species))
  k <- length(sp)
  a <- matrix(0, k, k, dimnames = list(sp, sp))
  up <- upper.tri(a)
  vals <- if (type == "none") numeric(sum(up)) else stats::runif(sum(up), 0, inter_max)
  if (type == "phylo") {
    d <- ape::cophenetic.phylo(species)[sp, sp]
    s <- (1 - d / max(d))[up]
    vals <- sort(vals)[rank(s, ties.method = "random")]
  }
  a[up] <- vals
  a <- a + t(a)
  diag(a) <- 1
  scale * a
}

#' competition by the overlap of Gaussian niches, 1 on the diagonal; product over axes
#' @noRd
trait_competition <- function(traits, inter_max) {
  if (is.null(traits)) stop("type = \"trait\" needs `traits`", call. = FALSE)
  several <- is.list(traits) && !is.data.frame(traits)
  niche_cor <- if (several) attr(traits, "niche_cor")
  axes <- if (several) traits else list(traits)
  axes <- lapply(axes, function(x) {
    if (!all(c("optimum", "sigma") %in% names(x))) {
      stop("`traits` needs columns `optimum` and `sigma`", call. = FALSE)
    }
    x
  })
  sp <- rownames(axes[[1]])
  unnamed <- is.null(sp) || (is.data.frame(axes[[1]]) && .row_names_info(axes[[1]]) < 0)
  if (unnamed) {
    sp <- paste0("sp", seq_len(nrow(axes[[1]])))
  } else if (length(axes) > 1) {
    axes <- trait_axes_aligned(axes, sp)
  }
  if (!is.null(niche_cor)) return(trait_competition_correlated(axes, sp, niche_cor, inter_max))
  log_overlap <- 0
  for (x in axes) {
    w <- outer(x$sigma^2, x$sigma^2, "+")
    d2 <- outer(x$optimum, x$optimum, "-")^2
    log_overlap <- log_overlap - ifelse(is.finite(w), niche_penalty(w, d2), 0)
  }
  a <- inter_max * exp(log_overlap)
  diag(a) <- 1
  dimnames(a) <- list(sp, sp)
  a
}

#' the tables of traits of several axes with their rows in the order of the species `sp`; stops
#' when the axes have different species
#' @noRd
trait_axes_aligned <- function(axes, sp) {
  if (!all(vapply(axes, function(x) setequal(rownames(x), sp), logical(1)))) {
    stop("sets of `traits` (one per axis) have different species", call. = FALSE)
  }
  lapply(axes, function(x) x[sp, , drop = FALSE])
}

#' competition with correlated axes of the niche: the overlap of the niches of two species is
#' exp(-d' (Omega_j + Omega_k)^-1 d / 2); returns the species x species matrix
#' @noRd
trait_competition_correlated <- function(axes, sp, niche_cor, inter_max) {
  q <- length(axes)
  cor <- niche_cor_species(niche_cor, sp, q)
  sigma <- matrix(vapply(axes, function(x) x$sigma, numeric(length(sp))), length(sp), q)
  optimum <- matrix(vapply(axes, function(x) x$optimum, numeric(length(sp))), length(sp), q)
  if (any(sigma == 0)) {
    stop("a niche breadth `sigma` of 0 is not available with correlated axes ",
         "(`niche_correlation()`)", call. = FALSE)
  }
  precision <- lapply(seq_along(sp), function(j) niche_precision(sigma[j, ], cor[[j]]))
  a <- diag(1, length(sp))
  for (j in seq_along(sp)[-1]) {
    for (k in seq_len(j - 1)) {
      d <- optimum[j, ] - optimum[k, ]
      m <- niche_overlap_matrix(precision[[j]], precision[[k]])
      a[j, k] <- a[k, j] <- inter_max * exp(-0.5 * sum((d %*% m) * d))
    }
  }
  dimnames(a) <- list(sp, sp)
  a
}

#' the penalty `d2 / (2 * w)` of the Gaussian niche overlap, from the summed variances `w`
#' and the squared distances `d2`: with `w = 0`, none at the same phenotype and an infinite
#' one elsewhere
#' @noRd
niche_penalty <- function(w, d2) {
  flat <- ifelse(d2 > 0, Inf, 0)
  ifelse(w > 0, d2 / (2 * w), flat)
}
