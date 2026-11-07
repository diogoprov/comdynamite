test_that("GAUSE removes, from the closest pair, the species with smaller summed distance", {
  # a and b are the closest pair; b is closer to the others (c, d), so it goes first
  x <- c(a = 0, b = 0.5, c = 3, d = 4)
  out <- assemble_community(x, size = 3, rule = "limiting_similarity")
  expect_equal(out$removed, "b")
  expect_setequal(out$community, c("a", "c", "d"))
  # next removal: c and d are now the closest pair; c has the smaller sum
  out2 <- assemble_community(x, size = 2, rule = "limiting_similarity")
  expect_equal(out2$removed, c("b", "c"))
})

test_that("limiting similarity spaces species more evenly than random draws", {
  set.seed(1)
  x <- stats::setNames(stats::runif(40), paste0("s", 1:40))
  nn <- function(v) min(diff(sort(v)))
  ls <- replicate(50, nn(x[assemble_community(x, 10, "limiting_similarity")$community]))
  nt <- replicate(50, nn(x[assemble_community(x, 10, "neutral")$community]))
  expect_gt(mean(ls), mean(nt))
})

test_that("min_distance stops removal at a limit to similarity", {
  x <- c(a = 0, b = 0.1, c = 1, d = 2.05, e = 3)
  out <- assemble_community(x, size = 1, rule = "limiting_similarity", min_distance = 0.5)
  expect_equal(out$removed, c("b"))
  expect_true(min(dist(x[out$community])) >= 0.5)
})

test_that("filtering keeps the species nearest the optimum", {
  x <- cbind(t1 = c(a = -2, b = -0.1, c = 0.2, d = 3), t2 = 0)
  expect_setequal(assemble_community(x, 2, "filtering")$community, c("b", "c"))
  expect_setequal(assemble_community(x, 2, "filtering", optimum = "derived")$community,
                  c("c", "d"))
  expect_setequal(assemble_community(x, 1, "filtering", optimum = c(-2, 0))$community, "a")
})

test_that("multivariate traits and output structure", {
  set.seed(2)
  tr <- ape::rcoal(30)
  x <- traits_ou(tr, n_traits = 3)
  expect_equal(dim(x), c(30, 3))
  for (rule in c("limiting_similarity", "filtering", "neutral")) {
    out <- assemble_community(x, 12, rule)
    expect_length(out$community, 12)
    expect_length(out$removed, 18)
    expect_setequal(c(out$community, out$removed), tr$tip.label)
  }
  expect_error(assemble_community(x, 31))
})

test_that("traits without row names get default species names", {
  set.seed(1)
  x <- matrix(rnorm(12), 6, 2)
  expect_warning(out <- assemble_community(x, size = 3, rule = "neutral"), "no names")
  expect_length(out$community, 3)
  expect_true(all(out$community %in% paste0("sp", 1:6)))
})
