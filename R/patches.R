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
#' with parts far apart, such as islands, may give patches that cross the gaps. With sparse
#' sites, as cells with records in a large region, the default `gap` follows the spacing of
#' the sites, so that neighbouring sites of an area form a patch instead of moving one by
#' one. Results should still be compared across values of `gap`.
#'
#' This method was built for this package and has no published evaluation. In
#' `data-raw/validation_edges.R`, with communities that followed a hidden autocorrelated
#' environment and areas without effect placed at random, it rejected 3.0% to 6.3% of the
#' tests at the 5% level in a square, a strip and an irregular band, against 44% to 55% with
#' free permutations (see the article on analysing clades). With autocorrelation at scales up
#' to that of the landscape it rejected 2.2% to 6.0%, and free permutations 49% to 75%.
#'
#' These rates average over places of areas, and a real analysis has one. Patches are placed
#' independently of each other, so the distance between the two sides of an event is not kept,
#' and sides tend to end close together. With areas fixed in a square of 20 by 20 sites, the
#' rate was 0.4% to 2% for adjacent areas, and 6% to 10% for areas in opposite corners or
#' halves of the landscape. The test is therefore liberal for sides far apart, as across a
#' barrier, and conservative for neighbouring sides. With 25% to 40% of the sites of that
#' square and areas in opposite corners (`data-raw/validation_gap.R`, 150 replicates and 99
#' permutations), it rejected 15% to 17% with 1.5 times the distance to the nearest site as
#' `gap`, the former default, and 5.6% to 5.9% with the current default, which follows the
#' fourth nearest site; those two cases were used to choose it. In sparse landscapes held out
#' from that choice (irregular areas with 40% of the sites, adjacent areas with 60%, random
#' blocks with 80% and irregular areas on seven clusters of sites) it rejected 2.1% to 5.9%,
#' against 2.6% to 8.4% with the former default. It is still far better than free
#' permutations, which rejected 48% to 74%, but p-values near the threshold need caution.
#' [rigid_permutations()] keeps those distances, and is the better choice for sides close
#' enough to have room to move.
#'
#' @param coords Matrix or data frame with the two coordinates of each site, in the order of
#'   the rows of `biogeography`.
#' @param biogeography Matrix of biogeographic events with one row per site and values of 1,
#'   -1 and 0, as given to [node_test()].
#' @param n_perm Number of permutations.
#' @param gap Largest distance between neighbouring sites of a patch. The default, 1.5 times
#'   the median distance from a site to its fourth nearest site, joins diagonal neighbours in
#'   a full grid of at least 8 by 8 sites, where the four nearest sites are one step away,
#'   and grows with the spacing of sparse sites, where it follows the distance between
#'   neighbouring sites rather than the nearest pair.
#' @return List with one integer matrix per event, each with one permutation per row and the
#'   names of the sites as column names, for `area_permutations` in [node_test()]. Rows only
#'   move the areas of that event, so they should not be used as `site_permutations`.
#' @seealso [rigid_permutations()] keeps the shape of the sides and the distances between
#'   them, and [torus_permutations()] gives a torus. [biogeography_matrix()] codes the events.
#' @examples
#' coords <- expand.grid(x = 1:10, y = 1:6)
#' event <- ifelse(coords$x <= 3 & coords$y <= 3, 1, ifelse(coords$x >= 8, -1, 0))
#' permutations <- patch_permutations(coords, cbind(barrier = event), n_perm = 3)
#'
#' # where the first permutation places the two sides
#' matrix(event[permutations$barrier[1, ]], 10)
#' @export
patch_permutations <- function(coords, biogeography, n_perm = 999, gap = NULL) {
  check_event_coords(coords, biogeography)
  coords <- as.matrix(coords)
  biogeography <- as.matrix(biogeography)
  stopifnot(n_perm >= 1)
  distances <- as.matrix(stats::dist(coords))
  if (is.null(gap)) {
    gap <- 1.5 * neighbour_distance(distances, 4)
    if (gap == 0) {
      stop("most sites share their coordinates with other sites: remove repeated sites or ",
           "give `gap`", call. = FALSE)
    }
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

# coordinates and events of the permutations of events
check_event_coords <- function(coords, biogeography) {
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
