#' Permutations that move biogeographic events as rigid bodies
#'
#' A restricted null model for the events of [node_test()] that keeps the shape of each side
#' and the distances between the sides. Each permutation moves all the sites of an event
#' together, by a rotation or a reflection and a translation. A move is kept when its sites
#' land inside the landscape: within `gap` of a site, or in a hole between sites, but not
#' beyond an edge. Sites that land where there is no site, as in a hole of a sparse
#' landscape, go to the nearest free site within twice `gap`, and a fraction `outside` of the
#' sites of an event may land beyond the landscape and be brought in the same way. Moves are
#' never carried across the landscape, as in a torus, which would split a side. This is the
#' torus-translation test of Harms et al. (2001) without the wrap-around.
#'
#' The moves are enumerated, so the permutations are distinct placements of the event, up to
#' the sign of the contrast, which does not change the test, and drawn at random when there
#' are more than `n_perm`. Sites on a regular lattice get all the moves of the lattice: its
#' symmetries with every translation in cells. Other coordinates get rotations in steps that
#' move the site farthest from the centre of the event by about `gap`, at least 4 and at most
#' 360 of them, each with and without a reflection, and the translations that put the site
#' nearest the centre on each site of the landscape; the other sites go to the nearest free
#' sites. Which moves fit is then fixed by the landscape, and only the draw among them
#' changes between runs. In an event without zeros where one side has up to a third of the
#' sites, only that side moves and the other is the rest of the landscape, as when one area
#' is compared with all the others.
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
#' across values. A large part in an irregular region may have few moves that fit. It then
#' stays nearly in place while the small parts move, which gives many placements that say
#' little, and the function warns.
#'
#' This method was built for this package and has no published evaluation. In
#' `data-raw/validation_rigid.R`, with communities that followed a hidden autocorrelated
#' environment and areas without effect, 150 replicates and 199 permutations, the default
#' rejected 3.7% to 4.7% of the tests at the 5% level with areas placed at random in a square
#' of 20 by 20 sites and in a strip of 50 by 8. With areas fixed in the square, it rejected
#' 2.2% to 3.5% for adjacent areas and for several small areas, and 4.2% to 4.8% for adjacent
#' areas, an L shape and two areas at a short distance when the hidden environment had an
#' exponential covariance. [patch_permutations()] rejected 1.4% to 6.4% in those cases, and is
#' conservative for adjacent areas. Where the event had room to move, power was equal to or
#' higher than that of [patch_permutations()], up to twice for adjacent areas (24% against 13%);
#' with few placements it was lower (4% against 13% for blocks in the band).
#'
#' With sparse sites (`data-raw/validation_gap.R`, 150 replicates and 99 permutations), the
#' default rejected 5.5% to 6.1% where the event had room to move: irregular areas with 40% of
#' the sites of the square, adjacent areas with 60% and random blocks with 80%. Areas in
#' opposite corners with 25% to 40% of the sites, and irregular areas on seven clusters of
#' sites, had no placement or fewer than five, so the default does not test them.
#'
#' Off a lattice (`data-raw/validation_offgrid.R`, 150 replicates and 99 permutations), the
#' default rejected 4.7% for adjacent areas in the square with jitter and 5.9% in 400 sites
#' placed at random, and 6.8% to 8.0% for two blocks placed at random in those landscapes, in
#' the band with jitter and in a crescent of random sites, where it is liberal.
#' [patch_permutations()] rejected 1.2% to 6.3% in those cases. Power for adjacent areas with
#' jitter was 23% against 12%. With 400 sites, the enumeration and 99 placements took 0.3 to
#' 2.8 s; with 750 sites placed at random and an event of 57 sites, 999 placements took 9 s,
#' and with 750 sites on a lattice with jitter and an event of 225 sites, 701 placements took
#' 12 min, since every site of a large event has to be moved to a free site in each try.
#'
#' Sides far apart remain a problem. In an irregular band most random areas had no placement
#' or fewer than 20, so the default does not test them. With `independent = 10`, random areas
#' in the three landscapes gave 3.9% to 5.5%, but areas in opposite corners of the square gave
#' 6% to 7% with autocorrelation over 3 or 4 sites and 12% over 12 sites, strips along opposite
#' edges 8% to 9% (2% to 3% with [patch_permutations()]), and areas in opposite corners with
#' sparse sites 10% to 12%.
#'
#' @param coords Matrix or data frame with the two coordinates of each site, in the order of
#'   the rows of `biogeography`. Distances are Euclidean, so they should be projected.
#' @param biogeography Matrix of biogeographic events with one row per site and values of 1,
#'   -1 and 0, as given to [node_test()].
#' @param n_perm Largest number of permutations. Fewer are returned when fewer distinct
#'   placements of the event exist.
#' @param gap A moved site counts as inside the landscape when it is closer than `gap` to a
#'   site, or in a hole narrower than twice `gap` between sites (four times `gap` off a
#'   lattice, where sites leave larger voids between neighbours), but not beyond an edge. A
#'   site that lands inside but not on a site goes to the nearest free site within twice
#'   `gap`. The default is the median distance from a site to its fourth nearest site: the
#'   step of a full grid, where it makes a site just beyond the edge count as outside, and
#'   the spacing of sparse sites, so that holes count as inside. It is not the `gap` of
#'   [patch_permutations()], which is 1.5 times this one by default.
#' @param outside Largest fraction of the sites of an event that may land outside the
#'   landscape in a move. Those sites go to the nearest free site within twice `gap`, and the
#'   move is dropped when one has none.
#' @param independent Distance beyond which parts of an event move independently of each
#'   other, while staying at least that far apart. The default, `Inf`, moves the whole event
#'   together and assumes nothing about the range of autocorrelation.
#' @param max_tries Largest number of moves tried: of the enumerated moves of one part, in
#'   random order, when placements fail or repeat a layout, or of the combinations of several
#'   independent parts, which cannot all be enumerated.
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
    gap <- neighbour_distance(distances, 4)
    if (gap == 0) {
      stop("most sites share their coordinates with other sites: remove repeated sites or ",
           "give `gap`", call. = FALSE)
    }
  }
  stopifnot(gap > 0)
  # strictly closer than `gap`, so that a site one step beyond a regular grid is outside
  gap <- gap * (1 - 1e-6)
  landscape <- landscape_of(coords, gap)
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
    drawn <- event_draws(coords, distances, landscape, groups, event, n_perm, gap, outside,
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

# median distance from a site to its k-th nearest site: the step of a full grid when k is at
# most 4, and the spacing of sparse sites
neighbour_distance <- function(distances, k = 4) {
  n <- nrow(distances)
  k <- min(k, n - 1)
  stats::median(apply(distances + diag(Inf, n), 1, function(d) sort(d)[k]))
}

# the geometry of the landscape: a regular lattice when the sites lie on one, else a raster
# with cells of half of `gap`. Both hold the cells inside the landscape: those near a site
# and those of holes between sites, but not the rim beyond an edge. On a lattice a cell is
# near a site when closer than `gap`; off a lattice, where sites leave larger voids between
# neighbours, when closer than twice `gap`
landscape_of <- function(coords, gap) {
  lattice_of(coords, gap) %||% raster_of(coords, gap)
}

# the cells of a raster that count as inside the landscape: those closer than `radius` to a
# site, and those of a hole narrower than twice `radius` between them. Cells beyond an edge
# are not inside. `cells` are the cells of the sites, `origin` the coordinates of the centre
# of cell (1, 1) and `steps` the size of a cell
inside_cells <- function(coords, cells, extent, origin, steps, radius) {
  reach <- ceiling(radius / steps)
  offsets <- as.matrix(expand.grid(-reach[1]:reach[1], -reach[2]:reach[2]))
  offsets <- offsets[(offsets[, 1] * steps[1])^2 + (offsets[, 2] * steps[2])^2 < radius^2, ,
                     drop = FALSE]
  # a cell is near a site when its centre is closer than `radius` to the site
  n <- nrow(coords)
  n_off <- nrow(offsets)
  candidate <- cells[rep(seq_len(n), each = n_off), , drop = FALSE] +
    offsets[rep(seq_len(n_off), n), , drop = FALSE]
  centre <- sweep(sweep(candidate - 1, 2, steps, "*"), 2, origin, "+")
  d2 <- rowSums((centre - coords[rep(seq_len(n), each = n_off), , drop = FALSE])^2)
  keep <- d2 < radius^2 & candidate[, 1] >= 1 & candidate[, 1] <= extent[1] &
    candidate[, 2] >= 1 & candidate[, 2] <= extent[2]
  near <- matrix(FALSE, extent[1], extent[2])
  near[candidate[keep, , drop = FALSE]] <- TRUE
  # a cell is inside when every cell within `radius` of it is near a site: holes between
  # sites stay inside, the rim beyond an edge does not
  inside <- near
  rows <- seq_len(extent[1])
  cols <- seq_len(extent[2])
  for (o in seq_len(n_off)) {
    from_rows <- rows + offsets[o, 1]
    from_cols <- cols + offsets[o, 2]
    ok_rows <- from_rows >= 1 & from_rows <= extent[1]
    ok_cols <- from_cols >= 1 & from_cols <= extent[2]
    shifted <- matrix(FALSE, extent[1], extent[2])
    shifted[rows[ok_rows], cols[ok_cols]] <- near[from_rows[ok_rows], from_cols[ok_cols]]
    inside <- inside & shifted
  }
  # the cells of the sites are inside, including those of isolated sites
  inside[cells] <- TRUE
  inside
}

# a regular lattice: the step of each axis, the cell of each site, the site of each cell,
# the cells inside the landscape and the symmetries. NULL when the coordinates do not lie
# on a lattice
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
  # moves may carry a site up to `gap` beyond the sites
  reach <- ceiling(gap / steps)
  margin <- reach
  extent <- apply(cells, 2, max) + 2 * margin + 1
  # rounded coordinates fall on a lattice far finer than the sites
  if (prod(extent) > 50 * nrow(coords)) return(NULL)
  cells <- sweep(cells, 2, margin + 1, "+")
  origin <- origin - margin * steps
  index <- matrix(0L, extent[1], extent[2])
  index[cells] <- seq_len(nrow(coords))
  near <- inside_cells(coords, cells, extent, origin, steps, gap)
  # the symmetries of the lattice: 8 when the steps are equal, else the 4 that keep the axes
  maps <- list(function(x, y) cbind(x, y), function(x, y) cbind(-x, y),
               function(x, y) cbind(x, -y), function(x, y) cbind(-x, -y))
  if (isTRUE(all.equal(steps[1], steps[2]))) {
    maps <- c(maps, list(function(x, y) cbind(y, x), function(x, y) cbind(-y, x),
                         function(x, y) cbind(y, -x), function(x, y) cbind(-y, -x)))
  }
  list(kind = "lattice", steps = steps, origin = origin, cells = cells, index = index,
       near = near, extent = extent, maps = maps, land = apply(cells, 2, range), reach = reach)
}

# a raster over coordinates that are not on a lattice, with square cells of half of `gap`,
# a margin of twice `gap` around the sites and the cells inside the landscape
raster_of <- function(coords, gap) {
  step <- gap / 2
  origin <- unname(apply(coords, 2, min)) - 2 * gap
  extent <- unname(floor((apply(coords, 2, max) + 2 * gap - origin) / step)) + 1
  cells <- unname(floor(sweep(sweep(coords, 2, origin), 2, step, "/"))) + 1
  near <- inside_cells(coords, cells, extent, origin, c(step, step), 2 * gap)
  list(kind = "raster", steps = c(step, step), origin = origin, near = near, extent = extent,
       gap = gap)
}

# the moves of a landscape that keep enough of a group inside it: one row per move
landscape_moves <- function(group, coords, landscape, outside) {
  if (landscape$kind == "lattice") {
    lattice_moves(group, landscape, outside)
  } else {
    raster_moves(group, coords, landscape, outside)
  }
}

# the image of a group under move `pick` of `moves`: coordinates and the site of each image
# when it lands on one, else 0
landscape_image <- function(group, moves, pick, coords, landscape) {
  if (landscape$kind == "lattice") {
    lattice_image(group, moves[pick, ], landscape)
  } else {
    raster_image(moves[pick, ], attr(moves, "orientation"), coords)
  }
}

# the moves of a lattice that keep enough of a group inside the landscape: the symmetry and
# the translation in cells
lattice_moves <- function(group, lattice, outside) {
  cells <- lattice$cells[group, , drop = FALSE]
  m <- nrow(cells)
  extent <- lattice$extent
  near <- lattice$near
  land <- lattice$land
  reach <- lattice$reach
  accepted <- list()
  for (s in seq_along(lattice$maps)) {
    image <- lattice$maps[[s]](cells[, 1], cells[, 2])
    box <- apply(image, 2, range)
    # translations that keep the group within `gap` of the sites
    tx <- seq(land[1, 1] - reach[1] - box[2, 1], land[2, 1] + reach[1] - box[1, 1])
    ty <- seq(land[1, 2] - reach[2] - box[2, 2], land[2, 2] + reach[2] - box[1, 2])
    for (dx in tx) {
      ix <- image[, 1] + dx
      in_x <- ix >= 1 & ix <= extent[1]
      if (sum(!in_x) > outside * m) next
      for (dy in ty) {
        iy <- image[, 2] + dy
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

# the image of a group under a move of the lattice
lattice_image <- function(group, move, lattice) {
  cells <- lattice$cells[group, , drop = FALSE]
  image <- lattice$maps[[move[1]]](cells[, 1], cells[, 2])
  image <- sweep(image, 2, move[2:3], "+")
  inside <- image[, 1] >= 1 & image[, 1] <= lattice$extent[1] &
    image[, 2] >= 1 & image[, 2] <= lattice$extent[2]
  site <- integer(nrow(cells))
  site[inside] <- lattice$index[image[inside, , drop = FALSE]]
  list(xy = sweep(sweep(image - 1, 2, lattice$steps, "*"), 2, lattice$origin, "+"),
       site = site)
}

# the orientations of a group off a lattice: rotations in steps that move the site farthest
# from the anchor by about `gap`, at least 4 and at most 360 of them, each with and without
# a reflection, as matrices. The anchor is the site of the group nearest its centre and the
# coordinates are relative to it. A group of one site has one orientation
group_orientations <- function(group, coords, gap) {
  from <- coords[group, , drop = FALSE]
  centre <- colMeans(from)
  anchor <- which.min(rowSums(sweep(from, 2, centre)^2))
  relative <- sweep(from, 2, from[anchor, ])
  if (nrow(from) == 1) {
    return(list(anchor = anchor, relative = relative, rotations = list(diag(2))))
  }
  radius <- sqrt(max(rowSums(relative^2)))
  n_angles <- min(360, max(4, 4 * ceiling(2 * pi * radius / gap / 4)))
  angles <- 2 * pi * (seq_len(n_angles) - 1) / n_angles
  rotations <- lapply(angles, function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2))
  reflected <- lapply(rotations, function(r) r * c(-1, -1, 1, 1))
  list(anchor = anchor, relative = relative, rotations = c(rotations, reflected))
}

# the moves of a group off a lattice that keep enough of it inside the landscape: each
# orientation with the translations that put the anchor on each site, as the index of the
# orientation and the site, with the orientations as an attribute
raster_moves <- function(group, coords, raster, outside) {
  orientation <- group_orientations(group, coords, raster$gap)
  m <- length(group)
  extent <- raster$extent
  accepted <- lapply(seq_along(orientation$rotations), function(o) {
    image <- orientation$relative %*% orientation$rotations[[o]]
    ix <- floor((outer(image[, 1], coords[, 1], "+") - raster$origin[1]) / raster$steps[1]) + 1
    iy <- floor((outer(image[, 2], coords[, 2], "+") - raster$origin[2]) / raster$steps[2]) + 1
    inside <- ix >= 1 & ix <= extent[1] & iy >= 1 & iy <= extent[2]
    inside[inside] <- raster$near[cbind(ix[inside], iy[inside])]
    far <- m - colSums(inside)
    sites <- which(far <= outside * m)
    if (length(sites) == 0) return(NULL)
    cbind(o, sites)
  })
  # the move that keeps the group in place always fits, so there is at least one
  structure(do.call(rbind, accepted), orientation = orientation)
}

# the image of a group under a move off a lattice: only the anchor lands on a site
raster_image <- function(move, orientation, coords) {
  image <- orientation$relative %*% orientation$rotations[[move[1]]]
  site <- integer(nrow(image))
  site[orientation$anchor] <- move[2]
  list(xy = sweep(image, 2, coords[move[2], ], "+"), site = site)
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
    # only the free sites within `reach` of the loose images can take them
    near <- which(free & coords[, 1] >= min(xy[, 1]) - reach & coords[, 1] <= max(xy[, 1]) + reach &
                    coords[, 2] >= min(xy[, 2]) - reach & coords[, 2] <= max(xy[, 2]) + reach)
    if (length(near) < length(loose)) return(NULL)
    d2 <- outer(xy[, 1], coords[near, 1], "-")^2 + outer(xy[, 2], coords[near, 2], "-")^2
    # ties between sites, as in the holes of a lattice, are broken at random
    d2 <- d2 + rep(stats::runif(length(near), 0, 1e-9), each = length(loose))
    d2[d2 > reach^2] <- Inf
    open <- rep(TRUE, length(near))
    for (row in order(d2[cbind(seq_along(loose), max.col(-d2, "first"))])) {
      d <- d2[row, ]
      d[!open] <- Inf
      site <- which.min(d)
      if (d[site] == Inf) return(NULL)
      target[loose[row]] <- near[site]
      open[site] <- FALSE
    }
  }
  target
}

# the distinct permutations of one event, or NULL when none differs from the observed
event_draws <- function(coords, distances, landscape, groups, event, n_perm, gap, outside,
                        independent, max_tries) {
  n <- nrow(coords)
  moving <- unlist(groups)
  rest <- setdiff(seq_len(n), moving)
  observed <- layout_key(event, seq_len(n))
  moves <- lapply(groups, landscape_moves, coords, landscape, outside)
  perms <- vector("list", n_perm)
  keys <- character(n_perm)
  found <- 0
  # one try: the groups placed, or NULL when one does not fit after the others. With one
  # group `pick` is its move; with several, each group draws its own until they stay apart
  one_try <- function(pick) {
    free <- rep(TRUE, n)
    target <- integer(n)
    for (g in seq_along(groups)) {
      group <- groups[[g]]
      attempt <- 0
      repeat {
        attempt <- attempt + 1
        if (length(groups) > 1) pick <- sample.int(nrow(moves[[g]]), 1)
        image <- landscape_image(group, moves[[g]], pick, coords, landscape)
        where <- place_image(image, coords, free, 2 * gap)
        if (!is.null(where)) {
          apart <- g == 1 || min(distances[where, target[target > 0], drop = FALSE]) > independent
          if (apart) break
        }
        if (attempt >= 30) return(NULL)
      }
      target[group] <- where
      free[where] <- FALSE
    }
    permutation <- integer(n)
    permutation[target[moving]] <- moving
    permutation[free] <- rest[sample.int(length(rest))]
    permutation
  }
  # the distinct layouts among the tries, until `n_perm` are found
  draw <- function(picks) {
    for (pick in picks) {
      permutation <- one_try(pick)
      if (is.null(permutation)) next
      key <- layout_key(event, permutation)
      if (key == observed || key %in% keys[seq_len(found)]) next
      found <<- found + 1
      perms[[found]] <<- permutation
      keys[found] <<- key
      if (found == n_perm) break
    }
  }
  if (length(groups) > 1) {
    draw(rep(NA_integer_, max_tries))
  } else {
    # the moves in random order: twice `n_perm` of them, then the others up to `max_tries`
    # when placements fail or repeat a layout
    rows <- nrow(moves[[1]])
    first <- sample.int(rows, min(rows, 2 * n_perm, max_tries))
    draw(first)
    if (found < n_perm && rows > length(first)) {
      others <- setdiff(seq_len(rows), first)
      others <- others[sample.int(length(others), min(length(others), max_tries - length(first)))]
      draw(others)
    }
  }
  if (found == 0) return(NULL)
  structure(do.call(rbind, perms[seq_len(found)]), moves = vapply(moves, nrow, integer(1)))
}
