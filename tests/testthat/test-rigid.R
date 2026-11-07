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

test_that("coordinates off a lattice work, and rounding does not make a huge lattice", {
  set.seed(4)
  coords <- as.matrix(square) + matrix(runif(288, -0.2, 0.2), ncol = 2)
  event <- area(2:4, 2:4) * 1 - area(6:8, 2:4) * 1
  drawn <- rigid_permutations(coords, cbind(e = event), n_perm = 10)$e
  expect_equal(nrow(drawn), 10)
  expect_true(all(apply(drawn, 1, function(row) all(sort(row) == 1:144))))
  expect_true(all(apply(drawn, 1, function(row) sum(event[row] == 1)) == 9))
  expect_false(is.null(lattice_of(as.matrix(square), 1)))
  rounded <- round(cbind(runif(50, 0, 40), runif(50, 0, 40)), 3)
  expect_null(lattice_of(rounded, 1))
})

test_that("a sparse lattice places sites in the holes within reach", {
  set.seed(5)
  keep <- sort(sample(144, 90))
  event <- ifelse(square$x[keep] <= 4 & square$y[keep] <= 4, 1,
                  ifelse(square$x[keep] %in% 6:9 & square$y[keep] <= 4, -1, 0))
  # with few holes the default `gap` counts them as outside
  expect_error(rigid_permutations(square[keep, ], cbind(e = event), n_perm = 10), "no placement")
  drawn <- suppressWarnings(rigid_permutations(square[keep, ], cbind(e = event), n_perm = 10,
                                               gap = 1.5))$e
  expect_true(all(apply(drawn, 1, function(row) all(sort(row) == 1:90))))
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
