#' Assemble a local community from a regional pool by rules
#'
#' Static community assembly, without population dynamics, following Kraft et
#' al. (2007). Species are removed from the pool one at a time until the
#' community has `size` species (or until the optional distance end point is
#' reached):
#'
#' * `"limiting_similarity"`: the GAUSE algorithm of Colwell & Winkler (1984),
#'   direct competition followed by diffuse competition. The pair of species
#'   closest in trait space is found and, of the two, the species with the
#'   smaller sum of distances to all remaining species (the greater "diffuse
#'   similarity") is removed. The process is repeated on the remaining
#'   species. Kraft et al. (2007) used this algorithm for their limiting
#'   similarity model.
#' * `"filtering"`: habitat filtering as in Kraft et al. (2007). The species
#'   farthest from a niche optimum in trait space are removed first.
#' * `"neutral"`: species are drawn at random, without regard to traits.
#'
#' Several traits can be used at once (one column each); ties are broken at random. The
#' distance between species, and between a species and the optimum, is chosen by `distance`:
#'
#' * `"euclidean"` (the default, as in Kraft et al. 2007): distances in the units of the
#'   traits, so a trait with a larger spread weighs more.
#' * `"standardized"`: Euclidean distance on traits divided by their standard deviation in
#'   the pool, so each trait weighs the same whatever its units.
#' * `"mahalanobis"`: Euclidean distance after removing the correlations between traits,
#'   \eqn{\sqrt{(x - y)^\top S^{-1} (x - y)}}{sqrt((x - y)' S^-1 (x - y))}, with \eqn{S} the
#'   covariance matrix of the traits in the pool,
#'   or `trait_cov`. In the filter, removing species by decreasing Mahalanobis distance to
#'   the optimum is removing them by increasing value of a multivariate Gaussian niche with
#'   covariance `S` centred on the optimum (Mahalanobis 1936, reprinted 2018).
#'
#' The choice is a statement about niche axes. With correlated traits, such as those of
#' [traits_ou()] with a rate matrix from [rate_matrix()], Euclidean distances, standardised or
#' not, count the shared direction of variation more than once: two traits that measure one
#' niche axis double the weight of that axis. `"mahalanobis"` treats each independent
#' direction of variation as one axis, and is the better choice when the traits are
#' correlated because they reflect the same niche dimension; `"standardized"` is the better
#' choice when each trait is a distinct niche dimension whatever its correlation with the
#' others. With uncorrelated traits on the same scale the three give the same order of
#' removal. The pool covariance needs more species than traits; with few species, or to use
#' the covariance of the process that generated the traits (for example the rate matrix,
#' which is free of the phylogenetic dependence among species), pass it as `trait_cov`.
#'
#' @param traits Numeric vector (one trait) or matrix (one column per trait),
#'   with species names as names or row names. For example, the output of
#'   [traits_ou()].
#' @param size Number of species in the local community.
#' @param rule Assembly rule.
#' @param optimum For `rule = "filtering"`: the niche optimum, either
#'   `"ancestral"` (the root state, 0 by default in [traits_ou()]),
#'   `"derived"` (the maximum of each trait in the pool, one extreme corner of
#'   trait space), or a numeric vector with one value per
#'   trait.
#' @param min_distance For `rule = "limiting_similarity"`: optional end point.
#'   Removal stops as soon as no two remaining species are closer than this
#'   (a limit to similarity), even if more than `size` species remain. In the units of
#'   `distance`.
#' @param distance How to measure distances in trait space: `"euclidean"`,
#'   `"standardized"` or `"mahalanobis"`; see Details.
#' @param trait_cov For `distance = "mahalanobis"`, an optional covariance matrix of the
#'   traits, one row and column per trait, used instead of the covariance in the pool.
#' @return An `mc_assembly` object: a list with `community` (names of the species kept),
#'   `removed` (names of the species removed, in order of removal), the `rule`, `distance` and
#'   `optimum` (for filtering) used, the `traits` of the pool as a matrix, and `metric`, the
#'   matrix \eqn{S} of the distances, \eqn{\sqrt{(x - y)^\top S^{-1} (x - y)}}{sqrt((x - y)'
#'   S^-1 (x - y))} (the identity for Euclidean distances). [plot_assembly()] draws it in trait
#'   space.
#' @references
#'   Colwell, R. K., & Winkler, D. W. (1984). A null model for null models in
#'   biogeography. In D. R. Strong, D. Simberloff, L. G. Abele, & A. B.
#'   Thistle (Eds.), *Ecological Communities: Conceptual Issues and the
#'   Evidence* (pp. 344--359). Princeton University Press.
#'   \doi{10.1515/9781400857081.344}
#'
#'   Kraft, N. J. B., Cornwell, W. K., Webb, C. O., & Ackerly, D. D. (2007).
#'   Trait evolution, community assembly, and the phylogenetic structure of
#'   ecological communities. *The American Naturalist*, 170, 271--283.
#'   \doi{10.1086/519400}
#'
#'   Mahalanobis, P. C. (2018). Reprint of: Mahalanobis, P. C. (1936) "On the generalised
#'   distance in statistics". *Sankhya A*, 80(Suppl. 1), 1--7.
#'   \doi{10.1007/s13171-019-00164-5}
#' @family rule-based assembly
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(32)
#' x <- traits_ou(tr, n_traits = 2)
#' assemble_community(x, size = 16, rule = "limiting_similarity")$community
#' assemble_community(x, size = 16, rule = "filtering", optimum = "derived")$community
#' # correlated traits: distances that discount the correlation
#' y <- traits_ou(tr, sigma2 = rate_matrix(c(1, 4), correlation = 0.9))
#' a <- assemble_community(y, size = 16, distance = "mahalanobis")
#' a
#' plot_assembly(a)
#' @export
assemble_community <- function(traits, size,
                               rule = c("limiting_similarity", "filtering", "neutral"),
                               optimum = "ancestral", min_distance = NULL,
                               distance = c("euclidean", "standardized", "mahalanobis"),
                               trait_cov = NULL) {
  rule <- match.arg(rule)
  distance <- match.arg(distance)
  x <- as.matrix(traits)
  if (!is.numeric(x)) stop("`traits` must be numeric (one column per trait)", call. = FALSE)
  if (anyNA(x)) stop("`traits` has missing values", call. = FALSE)
  if (is.null(rownames(x))) {
    named <- if (is.null(dim(traits))) names(traits)
    if (is.null(named)) warning("`traits` has no names of species: using sp1, sp2, ...", call. = FALSE)
    rownames(x) <- named %||% paste0("sp", seq_len(nrow(x)))
  }
  sp <- rownames(x)
  n <- nrow(x)
  if (!is.numeric(size) || length(size) != 1 || is.na(size) || size < 1 || size > n) {
    stop("`size` must be one number between 1 and the ", n, " species of `traits`", call. = FALSE)
  }
  w <- trait_whitening(x, distance, trait_cov)  # Euclidean distances on x %*% w
  opt <- NULL
  removed <- switch(rule,
    neutral = sample(sp, n - size),
    filtering = {
      opt <- if (is.numeric(optimum)) optimum else switch(
        match.arg(optimum, c("ancestral", "derived")),
        ancestral = rep(0, ncol(x)),
        derived = apply(x, 2, max))
      if (length(opt) != ncol(x)) {
        stop("`optimum` needs one value per trait (", ncol(x), ")", call. = FALSE)
      }
      d <- sqrt(rowSums((sweep(x, 2, opt) %*% w)^2))
      ord <- order(-d, stats::runif(n))  # farthest first, random tie-break
      sp[ord[seq_len(n - size)]]
    },
    limiting_similarity = gause_removal(x %*% w, size, min_distance))
  structure(list(community = setdiff(sp, removed), removed = removed, rule = rule,
                 distance = distance, optimum = opt, traits = x,
                 metric = solve(tcrossprod(w))),
            class = "mc_assembly")
}

#' @rdname assemble_community
#' @param x An `mc_assembly` object from `assemble_community()`.
#' @param ... Not used.
#' @export
print.mc_assembly <- function(x, ...) {
  cat("<mc_assembly>", x$rule, "with", x$distance, "distances:", length(x$community), "of",
      nrow(x$traits), "species kept\n")
  cat("community:", paste(utils::head(x$community, 10), collapse = ", "),
      if (length(x$community) > 10) "...", "\n")
  invisible(x)
}

# matrix w such that Euclidean distances between the rows of x %*% w are the chosen distances
trait_whitening <- function(x, distance, trait_cov = NULL) {
  p <- ncol(x)
  if (!is.null(trait_cov) && distance != "mahalanobis") {
    stop("`trait_cov` is used only with distance = \"mahalanobis\"", call. = FALSE)
  }
  switch(distance,
    euclidean = diag(p),
    standardized = {
      s <- apply(x, 2, stats::sd)
      if (anyNA(s) || any(s == 0)) {
        stop("standardized distances need traits that vary in the pool", call. = FALSE)
      }
      diag(1 / s, p)
    },
    mahalanobis = {
      cv <- trait_cov %||% {
        if (nrow(x) <= p) {
          stop("the pool needs more species than traits for Mahalanobis distances; ",
               "give `trait_cov`", call. = FALSE)
        }
        stats::cov(x)
      }
      cv <- as.matrix(cv)
      if (!is.numeric(cv) || any(dim(cv) != p) || anyNA(cv) || !isSymmetric(unname(cv))) {
        stop("`trait_cov` must be a symmetric matrix with one row and column per trait",
             call. = FALSE)
      }
      u <- tryCatch(chol(cv), error = function(e) {
        stop("the covariance matrix of the traits is singular: Mahalanobis distances ",
             "need traits that are not linear combinations of each other", call. = FALSE)
      })
      backsolve(u, diag(p))  # (x %*% u^-1) has identity covariance when cv is that of x
    })
}

# GAUSE (Colwell & Winkler 1984): remove, from the closest pair, the species
# with the smaller sum of distances to the remaining species. The nearest neighbour of each
# species is kept up to date, so each step costs O(n) instead of O(n^2); the closest pairs are
# listed in the order of which(arr.ind = TRUE) on the remaining species, so the random draws,
# and the species removed, are the same as when scanning the whole matrix at each step
gause_removal <- function(x, size, min_distance = NULL) {
  d <- as.matrix(stats::dist(x))
  diag(d) <- Inf
  n <- nrow(d)
  keep <- rep(TRUE, n)
  nn <- apply(d, 1, min)  # distance to the nearest remaining species
  nn_at <- max.col(-d, ties.method = "first")
  nn_at[!is.finite(nn)] <- NA
  removed <- character(0)
  while (sum(keep) > size) {
    dmin <- min(nn[keep])
    if (!is.null(min_distance) && dmin >= min_distance) break
    rows <- which(keep & nn == dmin)
    pairs <- do.call(rbind, lapply(rows, function(r) cbind(r, which(keep & d[r, ] == dmin))))
    pairs <- pairs[order(pairs[, 2], pairs[, 1]), , drop = FALSE]
    pair <- pairs[sample.int(nrow(pairs), 1), ]
    tot <- vapply(pair, function(i) sum(d[i, keep & seq_len(n) != i]), 0)
    victim <- if (tot[1] == tot[2]) pair[sample.int(2, 1)] else pair[which.min(tot)]
    keep[victim] <- FALSE
    removed <- c(removed, rownames(d)[victim])
    d[victim, ] <- Inf
    d[, victim] <- Inf
    stale <- which(keep & nn_at == victim)
    for (i in stale) {
      nn_at[i] <- which.min(d[i, ])
      nn[i] <- d[i, nn_at[i]]
    }
  }
  removed
}
