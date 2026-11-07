square <- expand.grid(x = 1:12, y = 1:12)
area <- function(xs, ys) square$x %in% xs & square$y %in% ys
apart <- function(event) {
  sort(round(as.vector(dist(as.matrix(square)[event != 0, ])), 6))
}

test_that("rigid_permutations() keeps shape and distances, without repeats", {
  set.seed(1)
  event <- (area(2:4, 2:3) | area(2, 4:5)) * 1 - area(7:8, 7:8) * 1
  drawn <- rigid_permutations(square, cbind(a = event), n_perm = 40, outside = 0)
  expect_named(drawn, "a")
  expect_equal(dim(drawn$a), c(40, 144))
  expect_true(all(apply(drawn$a, 1, function(row) all(sort(row) == 1:144))))
  moved <- apply(drawn$a, 1, function(row) event[row])
  expect_true(all(apply(moved, 2, function(m) identical(apart(m), apart(event)))))
  # layouts are distinct up to the sign, and none is the observed
  keys <- apply(moved, 2, function(m) paste(m * m[m != 0][1], collapse = ""))
  expect_equal(anyDuplicated(keys), 0)
  expect_false(paste(event * event[event != 0][1], collapse = "") %in% keys)
})

test_that("few placements give fewer rows and a warning, and none an error", {
  halves <- ifelse(square$x <= 6, 1, -1)
  expect_warning(drawn <- rigid_permutations(square, cbind(h = halves), n_perm = 9),
                 "h has 1")
  expect_equal(nrow(drawn$h), 1)
  row <- expand.grid(x = 1:9, y = 1)
  ends <- ifelse(row$x <= 4, 1, ifelse(row$x >= 6, -1, 0))
  expect_error(rigid_permutations(row, cbind(e = ends), n_perm = 9, outside = 0),
               "no placement of e")
})

test_that("the larger side is what is left when one area faces all others", {
  set.seed(2)
  event <- ifelse(area(3:5, 3:5), 1, -1)
  drawn <- rigid_permutations(square, cbind(e = event), n_perm = 20)$e
  expect_true(all(apply(drawn, 1, function(row) sum(event[row] == 1)) == 9))
  expect_true(all(apply(drawn, 1, function(row) sum(event[row] == -1)) == 135))
})

test_that("independent parts move on their own and stay apart", {
  set.seed(3)
  event <- area(1:3, 1:3) * 1 - area(10:12, 10:12) * 1
  expect_lt(nrow(suppressWarnings(rigid_permutations(square, cbind(e = event), 30))$e), 30)
  drawn <- rigid_permutations(square, cbind(e = event), n_perm = 30, independent = 5)$e
  expect_equal(nrow(drawn), 30)
  xy <- as.matrix(square)
  gaps <- apply(drawn, 1, function(row) {
    m <- event[row]
    min(as.matrix(dist(xy))[m == 1, m == -1])
  })
  expect_true(all(gaps > 5))
})

test_that("a part with few moves gives a warning", {
  event <- area(1:12, 1:9) * 1 - area(1:3, 12) * 1
  # steps differ between the axes, so this lattice has four symmetries
  wide <- transform(square, x = 2 * x)
  expect_warning(rigid_permutations(wide, cbind(e = event), n_perm = 5, independent = 2.5,
                                    outside = 0),
                 "with 108 sites has only 16 moves")
})

test_that("rigid_permutations() rejects wrong input", {
  event <- cbind(rep(c(1, 0, -1), 3))
  coords <- expand.grid(1:3, 1:3)
  expect_error(rigid_permutations(coords[-1, ], event), "rows")
  expect_error(rigid_permutations(cbind(coords, 1), event), "two numeric")
  expect_error(rigid_permutations(coords, event * 2), "values of 1")
  expect_error(rigid_permutations(cbind(rep(1, 9), rep(1, 9)), event), "share their coordinates")
})

test_that("node_test() reports the permutations of each event", {
  set.seed(6)
  tree <- ape::rcoal(6)
  comm <- matrix(rpois(864, 3), 144, dimnames = list(paste0("s", 1:144), tree$tip.label))
  events <- cbind(a = area(2:4, 2:4) * 1 - area(6:8, 2:4) * 1,
                  b = ifelse(square$x <= 6, 1, -1))
  rownames(events) <- rownames(comm)
  drawn <- suppressWarnings(rigid_permutations(square, events, n_perm = 9))
  out <- node_test(comm, node_contrasts(tree), cbind(x = rnorm(144)), events, n_perm = 9,
                   area_permutations = drawn)
  expect_equal(out$permutations$areas, c(a = 9L, b = 1L))
  expect_output(print(out), "of areas a 9, b 1")
})

test_that("off a lattice, the moves are enumerated and the placements keep the sides", {
  set.seed(4)
  coords <- as.matrix(square) + matrix(runif(288, -0.2, 0.2), ncol = 2)
  event <- area(2:4, 2:4) * 1 - area(6:8, 2:4) * 1
  drawn <- rigid_permutations(coords, cbind(e = event), n_perm = 10)$e
  expect_equal(nrow(drawn), 10)
  expect_true(all(apply(drawn, 1, function(row) all(sort(row) == 1:144))))
  expect_true(all(apply(drawn, 1, function(row) sum(event[row] == 1)) == 9))
  # each moved site is within twice `gap` of where the rigid move put it, so the pairwise
  # distances of the sides change by at most four times `gap`
  d <- as.matrix(dist(coords))
  gap <- neighbour_distance(d, 4)
  observed <- sort(d[event != 0, event != 0][upper.tri(matrix(0, 18, 18))])
  for (row in seq_len(10)) {
    moved <- event[drawn[row, ]]
    placed <- sort(d[moved != 0, moved != 0][upper.tri(matrix(0, 18, 18))])
    expect_lt(max(abs(placed - observed)), 4 * gap)
  }
  # the moves do not depend on the seed: the same landscape gives the same count
  landscape <- landscape_of(coords, gap)
  expect_equal(landscape$kind, "raster")
  moves <- landscape_moves(which(event != 0), coords, landscape, 0.1)
  expect_identical(moves, landscape_moves(which(event != 0), coords, landscape, 0.1))
  expect_gt(nrow(moves), 10)
  # the orientations: a multiple of four rotations, each with a reflection
  orientation <- attr(moves, "orientation")
  expect_equal(length(orientation$rotations) %% 8, 0)
  expect_equal(orientation$rotations[[1]], diag(2))
  # a group of one site has one orientation
  expect_length(group_orientations(1, coords, gap)$rotations, 1)
  expect_false(is.null(lattice_of(as.matrix(square), 1)))
  rounded <- round(cbind(runif(50, 0, 40), runif(50, 0, 40)), 3)
  expect_null(lattice_of(rounded, 1))
})

test_that("a sparse lattice places sites in the holes within reach", {
  set.seed(5)
  keep <- sort(sample(144, 90))
  event <- ifelse(square$x[keep] <= 4 & square$y[keep] <= 4, 1,
                  ifelse(square$x[keep] %in% 6:9 & square$y[keep] <= 4, -1, 0))
  # the default `gap` follows the fourth nearest site, so the holes count as inside
  drawn <- rigid_permutations(square[keep, ], cbind(e = event), n_perm = 10)$e
  expect_equal(nrow(drawn), 10)
  expect_true(all(apply(drawn, 1, function(row) all(sort(row) == 1:90))))
  # the rim beyond the edge does not count as inside, the holes do
  coords <- as.matrix(square[keep, ])
  lattice <- lattice_of(coords, 1.5)
  expect_equal(lattice$kind, "lattice")
  rim <- cbind(lattice$cells[, 1], 1)
  expect_false(any(lattice$near[rim]))
  holes <- lattice$near & lattice$index == 0
  expect_gt(sum(holes), 0)
})

test_that("the inside of a landscape keeps holes and drops the rim", {
  # a lattice of 5 by 5 with one hole: the hole is inside, the rim is not
  coords <- as.matrix(expand.grid(1:5, 1:5))[-13, ]
  lattice <- lattice_of(coords, 1.5 * (1 - 1e-6))
  expect_true(lattice$near[lattice$cells[1, 1] + 2, lattice$cells[1, 2] + 2])
  expect_equal(sum(lattice$near), 25)
  # an isolated site stays inside off a lattice
  coords <- rbind(cbind(runif(20, 0, 5), runif(20, 0, 5)), c(50, 50))
  raster <- raster_of(coords, 1)
  cell <- floor((coords[21, ] - raster$origin) / raster$steps) + 1
  expect_true(raster$near[cell[1], cell[2]])
})

test_that("images without a free site within reach are dropped", {
  image <- list(xy = rbind(c(0, 0), c(10, 10)), site = c(0L, 0L))
  coords <- rbind(c(0, 0), c(0.5, 0), c(10, 10))
  expect_null(place_image(image, coords, c(TRUE, TRUE, FALSE), 1))
  expect_equal(place_image(image, coords, c(FALSE, TRUE, TRUE), 1), c(2L, 3L))
  expect_null(place_image(image, coords, c(FALSE, FALSE, TRUE), 1))
})

test_that("coordinates that only look like a lattice are not treated as one", {
  # steps within the tolerance, but drifting away from the lattice
  drift <- cumsum(c(0, rep(c(1, 1 + 5e-7), 6)))
  expect_null(lattice_of(cbind(drift, drift), 1))
  # a true lattice, but far finer than its sites
  expect_null(lattice_of(cbind(c(0, 0.5, 100, 200), c(0, 0.5, 100, 200)), 1))
})

test_that("an orientation that fits nowhere has no moves", {
  set.seed(10)
  strip <- as.matrix(expand.grid(1:15, 1:3)) + matrix(runif(90, -0.2, 0.2), ncol = 2)
  event <- ifelse(strip[, 1] <= 4, 1, ifelse(strip[, 1] >= 7 & strip[, 1] <= 9, -1, 0))
  landscape <- landscape_of(strip, 1)
  moves <- landscape_moves(which(event != 0), strip, landscape, 0)
  # the rotations by a quarter turn do not fit a strip of three rows
  expect_lt(length(unique(moves[, 1])), length(attr(moves, "orientation")$rotations))
  drawn <- suppressWarnings(rigid_permutations(strip, cbind(e = event), n_perm = 10, outside = 0))$e
  expect_true(all(apply(drawn, 1, function(row) sum(event[row] == 1)) == sum(event == 1)))
})

test_that("off a lattice, an event with no other placement stops", {
  set.seed(7)
  row <- cbind(1:9, runif(9, -0.1, 0.1))
  ends <- ifelse(row[, 1] <= 4, 1, ifelse(row[, 1] >= 6, -1, 0))
  expect_error(rigid_permutations(row, cbind(e = ends), n_perm = 5, outside = 0),
               "no placement of e")
})

test_that("more moves are tried when placements fail or repeat", {
  set.seed(8)
  # a crowded event on sparse sites: many moves land on holes beyond reach
  keep <- sort(sample(144, 100))
  coords <- as.matrix(square[keep, ])
  event <- ifelse(square$x[keep] <= 6 & square$y[keep] <= 6, 1,
                  ifelse(square$x[keep] >= 7 & square$y[keep] <= 6, -1, 0))
  few <- suppressWarnings(rigid_permutations(coords, cbind(e = event), n_perm = 40,
                                             max_tries = 80))$e
  more <- suppressWarnings(rigid_permutations(coords, cbind(e = event), n_perm = 40))$e
  expect_lte(nrow(few), nrow(more))
  expect_gt(nrow(more), 1)
})

test_that("patch_permutations() follows the fourth nearest site with sparse sites", {
  set.seed(9)
  keep <- sort(sample(400, 100))
  coords <- as.matrix(expand.grid(1:20, 1:20))[keep, ]
  d <- as.matrix(dist(coords))
  gap <- 1.5 * neighbour_distance(d, 4)
  expect_gt(gap, 2)
  event <- ifelse(coords[, 1] <= 6 & coords[, 2] <= 6, 1, ifelse(coords[, 1] >= 15 & coords[, 2] >= 15, -1, 0))
  drawn <- patch_permutations(coords, cbind(e = event), n_perm = 20)$e
  # patches are those of the default gap: sites of a side stay within reach of each other
  n_patches <- length(area_patches(event == 1, d, gap)) + length(area_patches(event == -1, d, gap))
  moved <- apply(drawn, 1, function(row) event[row])
  expect_true(all(colSums(moved == 1) == sum(event == 1)))
  expect_lt(n_patches, sum(event != 0))
  # on a full grid of at least 8 by 8 the default is still 1.5; with fewer than five sites
  # the nearest that exists is used
  expect_equal(1.5 * neighbour_distance(as.matrix(dist(expand.grid(1:8, 1:8))), 4), 1.5)
  expect_equal(neighbour_distance(as.matrix(dist(cbind(c(0, 1, 3), 0))), 4), 3)
})
