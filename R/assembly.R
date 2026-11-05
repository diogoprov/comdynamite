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
#' Distances are Euclidean, so several traits can be used at once (one column
#' each). Ties are broken at random.
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
#'   (a limit to similarity), even if more than `size` species remain.
#' @return A list with `community` (names of the species kept) and `removed`
#'   (names of the species removed, in order of removal).
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
#' @examples
#' set.seed(1)
#' tr <- ape::rcoal(32)
#' x <- traits_ou(tr, n_traits = 2)
#' assemble_community(x, size = 16, rule = "limiting_similarity")$community
#' assemble_community(x, size = 16, rule = "filtering", optimum = "derived")$community
#' @export
assemble_community <- function(traits, size,
                               rule = c("limiting_similarity", "filtering", "neutral"),
                               optimum = "ancestral", min_distance = NULL) {
  rule <- match.arg(rule)
  x <- as.matrix(traits)
  if (is.null(rownames(x))) rownames(x) <- names(traits) %||% paste0("sp", seq_len(nrow(x)))
  sp <- rownames(x)
  n <- nrow(x)
  stopifnot(size >= 1, size <= n)
  removed <- switch(rule,
    neutral = sample(sp, n - size),
    filtering = {
      opt <- if (is.numeric(optimum)) optimum else switch(
        match.arg(optimum, c("ancestral", "derived")),
        ancestral = rep(0, ncol(x)),
        derived = apply(x, 2, max))
      stopifnot(length(opt) == ncol(x))
      d <- sqrt(colSums((t(x) - opt)^2))
      ord <- order(-d, stats::runif(n))  # farthest first, random tie-break
      sp[ord[seq_len(n - size)]]
    },
    limiting_similarity = gause_removal(x, size, min_distance))
  list(community = setdiff(sp, removed), removed = removed)
}

# GAUSE (Colwell & Winkler 1984): remove, from the closest pair, the species
# with the smaller sum of distances to the remaining species
gause_removal <- function(x, size, min_distance = NULL) {
  d <- as.matrix(stats::dist(x))
  diag(d) <- Inf
  keep <- rep(TRUE, nrow(d))
  removed <- character(0)
  while (sum(keep) > size) {
    sub <- d[keep, keep, drop = FALSE]
    dmin <- min(sub)
    if (!is.null(min_distance) && dmin >= min_distance) break
    idx <- which(keep)
    pairs <- which(sub == dmin, arr.ind = TRUE)
    pair <- idx[pairs[sample.int(nrow(pairs), 1), ]]
    tot <- vapply(pair, function(i) sum(d[i, keep & seq_along(keep) != i]), 0)
    victim <- if (tot[1] == tot[2]) pair[sample.int(2, 1)] else pair[which.min(tot)]
    keep[victim] <- FALSE
    removed <- c(removed, rownames(d)[victim])
  }
  removed
}
