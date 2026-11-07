#' Permutations that move the patches of biogeographic events
#'
#' A restricted null model for the events of [node_test()] in landscapes with edges, such as
#' refugia in a biome or islands in an archipelago. Each area is split into patches of
#' neighbouring sites. In each permutation, every patch is placed around a site drawn at
#' random, on the sites nearest to it that are still free, from the smallest patch to the
#' largest. The number and size of the patches are kept, and so is the fact that sites of a
#' patch are neighbours, while their place changes. Permutations of a torus cannot do this,
#' because a landscape with edges does not wrap around.
#'
#' Patches take a compact form near their random site, so their shape is not kept, and a
#' patch next to an edge or to another patch may take an elongated form or be split.
#' Distances are Euclidean, so coordinates should be projected. A landscape with holes or
#' with parts far apart, such as islands, may give patches that cross the gaps.
#'
#' This method was built for this package and has no published evaluation. In
#' `data-raw/validation_edges.R`, with communities that followed a hidden autocorrelated
#' environment and areas without effect, it rejected 3.0% to 6.3% of the tests at the 5% level
#' in a square, a strip and an irregular band, against 44% to 55% with free permutations (see
#' the article on analysing clades). With autocorrelation at scales up to that of the
#' landscape it rejected 2.2% to 6.0%, and free permutations 49% to 75%. Areas much larger
#' than those simulated were not evaluated.
#'
#' @param coords Matrix or data frame with the two coordinates of each site, in the order of
#'   the rows of `biogeography`.
#' @param biogeography Matrix of biogeographic events with one row per site and values of 1,
#'   -1 and 0, as given to [node_test()].
#' @param n_perm Number of permutations.
#' @param gap Largest distance between neighbouring sites of a patch. The default, 1.5 times
#'   the median distance to the nearest site, joins diagonal neighbours in a regular grid.
#' @return List with one integer matrix per event, each with one permutation per row and the
#'   names of the sites as column names, for `area_permutations` in [node_test()]. Rows only
#'   move the areas of that event, so they should not be used as `site_permutations`.
#' @seealso [torus_permutations()] for a torus. [biogeography_matrix()] codes the events.
#' @examples
#' coords <- expand.grid(x = 1:10, y = 1:6)
#' event <- ifelse(coords$x <= 3 & coords$y <= 3, 1, ifelse(coords$x >= 8, -1, 0))
#' permutations <- patch_permutations(coords, cbind(barrier = event), n_perm = 3)
#'
#' # where the first permutation places the two sides
#' matrix(event[permutations$barrier[1, ]], 10)
#' @export
patch_permutations <- function(coords, biogeography, n_perm = 999, gap = NULL) {
  check_numeric_columns(coords, "`coords`", "Keep only the two coordinates.")
  coords <- as.matrix(coords)
  biogeography <- as.matrix(biogeography)
  if (!is.numeric(coords) || ncol(coords) != 2 || anyNA(coords)) {
    stop("`coords` needs two numeric columns without missing values", call. = FALSE)
  }
  if (nrow(coords) != nrow(biogeography)) {
    stop("`coords` has ", nrow(coords), " rows and `biogeography` has ", nrow(biogeography),
         call. = FALSE)
  }
  if (!all(biogeography %in% c(-1, 0, 1))) {
    stop("`biogeography` needs values of 1, -1 and 0", call. = FALSE)
  }
  stopifnot(n_perm >= 1)
  distances <- as.matrix(stats::dist(coords))
  if (is.null(gap)) {
    nearest <- apply(distances + diag(Inf, nrow(distances)), 1, min)
    gap <- 1.5 * stats::median(nearest)
  }
  stopifnot(gap > 0)
  sites <- rownames(biogeography)
  events <- lapply(seq_len(ncol(biogeography)), function(e) {
    groups <- c(area_patches(biogeography[, e] == 1, distances, gap),
                area_patches(biogeography[, e] == -1, distances, gap))
    groups <- groups[order(lengths(groups))]
    drawn <- t(vapply(seq_len(n_perm), function(r) place_patches(groups, distances),
                      integer(nrow(distances))))
    colnames(drawn) <- sites
    drawn
  })
  names(events) <- colnames(biogeography) %||% paste0("event", seq_along(events))
  events
}

# groups of sites of an area linked by steps no longer than `gap`
area_patches <- function(member, distances, gap) {
  members <- which(member)
  linked <- distances[members, members, drop = FALSE] <= gap
  label <- integer(length(members))
  for (start in seq_along(members)) {
    if (label[start] > 0) next
    label[start] <- start
    queue <- start
    while (length(queue) > 0) {
      reached <- which(linked[queue[1], ] & label == 0)
      label[reached] <- start
      queue <- c(queue[-1], reached)
    }
  }
  unname(split(members, label))
}

# one permutation: a patch goes to the free sites nearest to a random free site
place_patches <- function(groups, distances) {
  free <- rep(TRUE, nrow(distances))
  permutation <- integer(nrow(distances))
  for (patch in groups) {
    open <- which(free)
    centre <- open[sample.int(length(open), 1)]
    closest <- order(distances[centre, open], stats::runif(length(open)))
    targets <- open[closest[seq_along(patch)]]
    permutation[targets] <- patch
    free[targets] <- FALSE
  }
  permutation[free] <- setdiff(seq_along(free), unlist(groups))
  permutation
}
