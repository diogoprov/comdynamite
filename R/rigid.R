#' Permutations that move biogeographic events as rigid bodies
#'
#' A restricted null model for the events of [node_test()] that keeps the shape of each side
#' and the distances between the sides. Each permutation moves all the sites of an event
#' together, by a rotation or a reflection and a translation. A move is kept when its sites
#' land on the landscape. Sites that land where there is no site, as in a hole of a sparse
#' landscape, go to the nearest free site within twice `gap`, and a fraction `outside` of the
#' sites of an event may land beyond the landscape and be brought in the same way. Moves are
#' never carried across the landscape, as in a torus, which would split a side. This is the
#' torus-translation test of Harms et al. (2001) without the wrap-around.
#'
#' Sites on a regular lattice get all the moves of the lattice, enumerated, so the
#' permutations are distinct placements of the event, up to the sign of the contrast, which
#' does not change the test. Other coordinates get random rotations and translations, with the
#' sites snapped to the nearest free sites. In an event without zeros where one side has up
#' to a third of the sites, only that side moves and the other is the rest of the landscape,
#' as when one area is compared with all the others.
#'
#' An event that spans the landscape, as two halves, two strips along opposite edges or areas
#' in opposite corners, has few placements, sometimes none, because a move would push a side
#' out. The function then returns fewer rows than `n_perm`, with a warning, or stops. The
#' smallest p-value in [node_test()] is 1 / (rows + 1), so with fewer than 19 rows no node can
#' reach the 5% level. This is a limit of the event, not of the function: a permutation test
#' that keeps the shape of two halves has almost nothing to permute.
#'
#' With `independent`, parts of the event farther apart than that distance move independently
#' of each other, but keep at least that distance between them. This gives placements to sides
#' far apart, as across a barrier, at the price of assuming that the community is not spatially
#' autocorrelated beyond `independent`. When that assumption fails the test is liberal, and no
#' value of `independent` can be checked from the data alone, so results should be compared
#' across values. Off a lattice, a large part in an irregular region seldom fits after a random
#' move. It then has few moves and stays nearly in place while the small parts move, which
#' gives many placements that say little, and the function warns.
#'
#' This method was built for this package and has no published evaluation. In
#' `data-raw/validation_rigid.R`, with communities that followed a hidden autocorrelated
#' environment and areas without effect, 150 replicates and 199 permutations, the default
#' rejected 3.7% to 4.7% of the tests at the 5% level with areas placed at random in a square
#' of 20 by 20 sites and in a strip of 50 by 8. With areas fixed in the square, it rejected
#' 2.2% to 3.5% for adjacent areas and for several small areas, and 4.2% to 4.8% for adjacent
#' areas, an L shape and two areas at a short distance when the hidden environment had an
#' exponential covariance. [patch_permutations()] rejected 1.4% to 6.4% in those cases, and is
#' conservative for adjacent areas. With 25% to 40% of the sites of the square and areas in
#' opposite corners, the default rejected less than 1%, with 13 to 28 placements, and
#' [patch_permutations()] 15% to 19%. Where the event had room to move, power was equal to or
#' higher than that of [patch_permutations()], up to twice for adjacent areas (24% against 13%);
#' with few placements it was lower (4% against 13% for blocks in the band).
#'
#' Sides far apart remain a problem. In an irregular band most random areas had no placement
#' or fewer than 20, so the default does not test them. With `independent = 10`, random areas
#' in the three landscapes gave 3.9% to 5.5%, but areas in opposite corners of the square gave
#' 6% to 7% with autocorrelation over 3 or 4 sites and 12% over 12 sites, strips along opposite
#' edges 8% to 9% (2% to 3% with [patch_permutations()]), areas in opposite corners with sparse
#' sites 8%, and the two farthest clusters of a clustered landscape 9%. On coordinates off a lattice, the only case evaluated, a square
#' with jitter and 60 replicates, rates were 2% to 9%, so this branch is less precise.
#'
#' @param coords Matrix or data frame with the two coordinates of each site, in the order of
#'   the rows of `biogeography`. Distances are Euclidean, so they should be projected.
#' @param biogeography Matrix of biogeographic events with one row per site and values of 1,
#'   -1 and 0, as given to [node_test()].
#' @param n_perm Largest number of permutations. Fewer are returned when fewer distinct
#'   placements of the event exist.
#' @param gap A moved site counts as inside the landscape when it is closer than `gap` to a
#'   site. The default, the 90th percentile of the distances to the nearest site, makes a
#'   site just beyond the edge of a regular grid count as outside and lets sites land in the
#'   holes of a sparse landscape. With few holes it is too small, as most sites have a
#'   neighbour one step away, and holes count as outside; `gap` should then be raised. It is
#'   not the `gap` of [patch_permutations()].
#' @param outside Largest fraction of the sites of an event that may land outside the
#'   landscape in a move. Those sites go to the nearest free site within twice `gap`, and the
#'   move is dropped when one has none.
#' @param independent Distance beyond which parts of an event move independently of each
#'   other, while staying at least that far apart. The default, `Inf`, moves the whole event
#'   together and assumes nothing about the range of autocorrelation.
#' @param max_tries Largest number of moves tried when they cannot all be enumerated, that
#'   is, with several independent parts or with coordinates that are not on a lattice.
#' @return List with one integer matrix per event, each with one permutation per row and the
#'   names of the sites as column names, for `area_permutations` in [node_test()]. Rows only
#'   move the sites of that event, so they should not be used as `site_permutations`.
#' @references Harms, K. E., Condit, R., Hubbell, S. P., & Foster, R. B. (2001). Habitat
#'   associations of trees and shrubs in a 50-ha neotropical forest plot. *Journal of
#'   Ecology*, 89, 947--959.
#' @seealso [patch_permutations()] moves each patch on its own, and [torus_permutations()]
#'   gives a torus. [biogeography_matrix()] codes the events.
#' @examples
#' coords <- expand.grid(x = 1:10, y = 1:6)
#' event <- ifelse(coords$x <= 3 & coords$y <= 3, 1, ifelse(coords$x >= 5 & coords$x <= 6, -1, 0))
#' permutations <- rigid_permutations(coords, cbind(barrier = event), n_perm = 3)
#'
#' # where the first permutation places the two sides
#' matrix(event[permutations$barrier[1, ]], 10)
#' @export
rigid_permutations <- function(coords, biogeography, n_perm = 999, gap = NULL, outside = 0.1,
                               independent = Inf, max_tries = 20 * n_perm) {
  check_event_coords(coords, biogeography)
  coords <- as.matrix(coords)
  biogeography <- as.matrix(biogeography)
  stopifnot(n_perm >= 1, outside >= 0, outside <= 1, independent > 0, max_tries >= 1)
  n <- nrow(coords)
  distances <- as.matrix(stats::dist(coords))
  if (is.null(gap)) {
    nearest <- apply(distances + diag(Inf, n), 1, min)
    gap <- stats::quantile(nearest, 0.9, names = FALSE)
    if (gap == 0) {
      stop("most sites share their coordinates with another site: remove repeated sites or ",
           "give `gap`", call. = FALSE)
    }
  }
  stopifnot(gap > 0)
  # strictly closer than `gap`, so that a site one step beyond a regular grid is outside
  gap <- gap * (1 - 1e-6)
  lattice <- lattice_of(coords, gap)
  sites <- rownames(biogeography)
  events <- colnames(biogeography) %||% paste0("event", seq_len(ncol(biogeography)))
  out <- lapply(seq_len(ncol(biogeography)), function(e) {
    event <- biogeography[, e]
    moving <- seq_len(n) %in% moving_sites(event)
    groups <- if (is.finite(independent)) {
      area_patches(moving, distances, independent)
    } else {
      list(which(moving))
    }
    drawn <- event_draws(coords, distances, lattice, groups, event, n_perm, gap, outside,
                         independent, max_tries)
    if (is.null(drawn)) {
      stop("no placement of ", events[e], " other than the observed one fits the landscape: ",
           "raise `outside` or `max_tries`, or lower `independent`", call. = FALSE)
    }
    moves <- attr(drawn, "moves")
    attr(drawn, "moves") <- NULL
    if (length(groups) > 1 && min(moves) < 20) {
      warning("a part of ", events[e], " with ", lengths(groups)[which.min(moves)], " sites has ",
              "only ", min(moves), " moves, so it stays nearly in place among the permutations: ",
              "lower `independent` or raise `max_tries`", call. = FALSE)
    }
    colnames(drawn) <- sites
    drawn
  })
  names(out) <- events
  few <- vapply(out, nrow, integer(1)) < n_perm
  if (any(few)) {
    counts <- vapply(out[few], nrow, integer(1))
    warning("fewer distinct placements than `n_perm`: ",
            paste0(names(counts), " has ", counts, collapse = ", "),
            "; the smallest attainable p-value is 1 / (placements + 1)", call. = FALSE)
  }
  out
}

# sites that move: the sites of the event, or only the smaller side when the other side is
# the rest of the landscape (an event without zeros where one side has up to a third of the
# sites), in which case the larger side is the complement of where the smaller one lands
moving_sites <- function(event) {
  if (any(event == 0)) return(which(event != 0))
  sides <- list(which(event == 1), which(event == -1))
  smaller <- sides[[which.min(lengths(sides))]]
  if (length(smaller) <= length(event) / 3) smaller else seq_along(event)
}

# a layout up to the sign of the contrast, which does not change the test
layout_key <- function(event, permutation) {
  moved <- event[permutation]
  first <- moved[moved != 0][1]
  paste(moved * first, collapse = "")
}

# a regular lattice: the step of each axis, the cell of each site, the site of each cell,
# the cells closer than `gap` to a site and the symmetries. NULL when the coordinates do not
# lie on a lattice
lattice_of <- function(coords, gap) {
  steps <- apply(coords, 2, function(x) {
    u <- sort(unique(x))
    if (length(u) < 2) return(NA_real_)
    d <- diff(u)
    step <- min(d)
    if (step <= 0 || any(abs(d / step - round(d / step)) > 1e-6)) return(NA_real_)
    step
  })
  if (anyNA(steps)) return(NULL)
  steps <- unname(steps)
  origin <- unname(apply(coords, 2, min))
  cells <- sweep(sweep(coords, 2, origin), 2, steps, "/")
  if (any(abs(cells - round(cells)) > 1e-6)) return(NULL)
  cells <- unname(round(cells))
  margin <- ceiling(gap / steps)
  extent <- apply(cells, 2, max) + 2 * margin + 1
  # rounded coordinates fall on a lattice far finer than the sites
  if (prod(extent) > 50 * nrow(coords)) return(NULL)
  index <- matrix(0L, extent[1], extent[2])
  index[cbind(cells[, 1] + margin[1] + 1, cells[, 2] + margin[2] + 1)] <- seq_len(nrow(coords))
  grid <- as.matrix(expand.grid(seq_len(extent[1]) - margin[1] - 1,
                                seq_len(extent[2]) - margin[2] - 1))
  grid_xy <- sweep(sweep(grid, 2, steps, "*"), 2, origin, "+")
  near <- matrix(FALSE, extent[1], extent[2])
  for (rows in split(seq_len(nrow(grid)), ceiling(seq_len(nrow(grid)) / 2000))) {
    d2 <- outer(grid_xy[rows, 1], coords[, 1], "-")^2 + outer(grid_xy[rows, 2], coords[, 2], "-")^2
    near[cbind(grid[rows, 1] + margin[1] + 1, grid[rows, 2] + margin[2] + 1)] <-
      apply(d2, 1, min) < gap^2
  }
  # the symmetries of the lattice: 8 when the steps are equal, else the 4 that keep the axes
  maps <- list(function(x, y) cbind(x, y), function(x, y) cbind(-x, y),
               function(x, y) cbind(x, -y), function(x, y) cbind(-x, -y))
  if (isTRUE(all.equal(steps[1], steps[2]))) {
    maps <- c(maps, list(function(x, y) cbind(y, x), function(x, y) cbind(-y, x),
                         function(x, y) cbind(y, -x), function(x, y) cbind(-y, -x)))
  }
  list(steps = steps, origin = origin, cells = cells, margin = margin, index = index,
       near = near, extent = extent, maps = maps)
}

# the moves of a lattice that keep enough of a group inside the landscape: one row per
# move, with the symmetry and the translation in cells
lattice_moves <- function(group, lattice, outside) {
  cells <- lattice$cells[group, , drop = FALSE]
  m <- nrow(cells)
  margin <- lattice$margin
  extent <- lattice$extent
  near <- lattice$near
  land <- apply(lattice$cells, 2, range)
  accepted <- list()
  for (s in seq_along(lattice$maps)) {
    image <- lattice$maps[[s]](cells[, 1], cells[, 2])
    box <- apply(image, 2, range)
    tx <- seq(land[1, 1] - margin[1] - box[2, 1], land[2, 1] + margin[1] - box[1, 1])
    ty <- seq(land[1, 2] - margin[2] - box[2, 2], land[2, 2] + margin[2] - box[1, 2])
    for (dx in tx) {
      ix <- image[, 1] + dx + margin[1] + 1
      in_x <- ix >= 1 & ix <= extent[1]
      if (sum(!in_x) > outside * m) next
      for (dy in ty) {
        iy <- image[, 2] + dy + margin[2] + 1
        inside <- in_x & iy >= 1 & iy <= extent[2]
        far <- m - sum(inside)
        if (far > outside * m) next
        far <- far + sum(!near[cbind(ix[inside], iy[inside])])
        if (far > outside * m) next
        accepted[[length(accepted) + 1]] <- c(s, dx, dy)
      }
    }
  }
  do.call(rbind, accepted)
}

# the image of a group under a move of the lattice: coordinates and the site of each image
# when it lands on one, else 0
lattice_image <- function(group, move, lattice) {
  cells <- lattice$cells[group, , drop = FALSE]
  image <- lattice$maps[[move[1]]](cells[, 1], cells[, 2])
  image <- sweep(image, 2, move[2:3], "+")
  ix <- image[, 1] + lattice$margin[1] + 1
  iy <- image[, 2] + lattice$margin[2] + 1
  inside <- ix >= 1 & ix <= lattice$extent[1] & iy >= 1 & iy <= lattice$extent[2]
  site <- integer(nrow(cells))
  site[inside] <- lattice$index[cbind(ix[inside], iy[inside])]
  list(xy = sweep(sweep(image, 2, lattice$steps, "*"), 2, lattice$origin, "+"), site = site)
}

# a pool of random rigid moves of a group of coordinates that are not on a lattice, each
# with the nearest sites of every moved site, so that placing a move is cheap. A move is
# kept when few enough of its sites land outside the landscape
continuous_pool <- function(group, coords, gap, outside, size, max_tries, k = 12) {
  from <- coords[group, , drop = FALSE]
  m <- nrow(from)
  n <- nrow(coords)
  k <- min(k, n)
  pool <- vector("list", size)
  found <- 0
  for (t in seq_len(max_tries)) {
    angle <- stats::runif(1, 0, 2 * pi)
    rotation <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2)
    if (stats::runif(1) < 0.5) rotation[, 1] <- -rotation[, 1]
    image <- sweep(from, 2, from[sample.int(m, 1), ]) %*% rotation
    image <- sweep(image, 2, coords[sample.int(n, 1), ], "+")
    d2 <- outer(image[, 1], coords[, 1], "-")^2 + outer(image[, 2], coords[, 2], "-")^2
    if (sum(apply(d2, 1, min) >= gap^2) > outside * m) next
    nearest <- apply(d2, 1, function(d) order(d)[seq_len(k)])
    nearest_d2 <- matrix(d2[cbind(rep(seq_len(m), each = k), as.vector(nearest))], k)
    found <- found + 1
    pool[[found]] <- list(nearest = nearest, d2 = nearest_d2)
    if (found == size) break
  }
  if (found == 0) return(NULL)
  pool[seq_len(found)]
}

# where the sites of a pooled move go: each to its nearest free site within `reach` among
# the sites kept for it, the nearest first; NULL when one has none
place_pooled <- function(move, free, reach) {
  m <- ncol(move$nearest)
  target <- integer(m)
  for (i in order(move$d2[1, ])) {
    sites <- move$nearest[, i]
    open <- free[sites] & move$d2[, i] <= reach^2
    if (!any(open)) return(NULL)
    site <- sites[which(open)[1]]
    target[i] <- site
    free[site] <- FALSE
  }
  target
}

# where the sites of an image go: to their site when it is free, else to the nearest free
# site within `reach`, the nearest images first; NULL when a site has none, so that a site
# is never carried to another part of the landscape
place_image <- function(image, coords, free, reach) {
  target <- image$site
  target[target > 0 & !free[pmax(target, 1)]] <- 0L
  free[target[target > 0]] <- FALSE
  loose <- which(target == 0)
  if (length(loose) > 0) {
    xy <- image$xy[loose, , drop = FALSE]
    d2 <- outer(xy[, 1], coords[, 1], "-")^2 + outer(xy[, 2], coords[, 2], "-")^2
    d2 <- d2 + stats::runif(length(d2), 0, 1e-9)
    for (row in order(apply(d2, 1, min))) {
      d <- d2[row, ]
      d[!free] <- Inf
      site <- which.min(d)
      if (d[site] > reach^2) return(NULL)
      target[loose[row]] <- site
      free[site] <- FALSE
    }
  }
  target
}

# the distinct permutations of one event, or NULL when none differs from the observed
event_draws <- function(coords, distances, lattice, groups, event, n_perm, gap, outside,
                        independent, max_tries) {
  n <- nrow(coords)
  moving <- unlist(groups)
  rest <- setdiff(seq_len(n), moving)
  observed <- layout_key(event, seq_len(n))
  if (!is.null(lattice)) {
    moves <- lapply(groups, lattice_moves, lattice, outside)
  } else {
    moves <- lapply(groups, continuous_pool, coords, gap, outside, 2 * n_perm, max_tries)
  }
  if (any(vapply(moves, is.null, logical(1)))) return(NULL)
  # one group: all its moves, or a sample of them when there are many
  enumerate <- length(groups) == 1
  if (enumerate) {
    rows <- if (!is.null(lattice)) nrow(moves[[1]]) else length(moves[[1]])
    chosen <- if (rows > 2 * n_perm) sample.int(rows, 2 * n_perm) else seq_len(rows)
    tries <- length(chosen)
  } else {
    tries <- max_tries
  }
  placement <- function(group, g, pick, free) {
    if (!is.null(lattice)) {
      place_image(lattice_image(group, moves[[g]][pick, ], lattice), coords, free, 2 * gap)
    } else {
      place_pooled(moves[[g]][[pick]], free, 2 * gap)
    }
  }
  perms <- vector("list", n_perm)
  keys <- character(n_perm)
  found <- 0
  for (t in seq_len(tries)) {
    free <- rep(TRUE, n)
    target <- integer(n)
    placed <- TRUE
    for (g in seq_along(groups)) {
      group <- groups[[g]]
      attempt <- 0
      repeat {
        attempt <- attempt + 1
        pick <- if (enumerate) chosen[t] else sample.int(NROW(moves[[g]]), 1)
        where <- placement(group, g, pick, free)
        if (!is.null(where)) {
          apart <- g == 1 || min(distances[where, target[target > 0], drop = FALSE]) > independent
          if (apart) break
        }
        if (attempt >= 30) {
          placed <- FALSE
          break
        }
      }
      if (!placed) break
      target[group] <- where
      free[where] <- FALSE
    }
    if (!placed) next
    permutation <- integer(n)
    permutation[target[moving]] <- moving
    permutation[free] <- rest[sample.int(length(rest))]
    key <- layout_key(event, permutation)
    if (key == observed || key %in% keys[seq_len(found)]) next
    found <- found + 1
    perms[[found]] <- permutation
    keys[found] <- key
    if (found == n_perm) break
  }
  if (found == 0) return(NULL)
  structure(do.call(rbind, perms[seq_len(found)]), moves = vapply(moves, NROW, integer(1)))
}
