test_that("nearest-neighbour and pairwise metrics follow their definitions", {
  x <- c(1, 2, 2, 5, 9)
  expect_equal(nn_distances(x), c(1, 0, 0, 3, 4))
  expect_equal(mnnd(x), 1.6)
  expect_equal(sdnnd(x), stats::sd(c(1, 0, 0, 3, 4)))
  expect_equal(cvnnd(x), stats::sd(c(1, 0, 0, 3, 4)) / 1.6)
  expect_equal(min_nnd(x), 0)
  expect_equal(max_nnd(x), 4)
  expect_equal(mnd(x), mean(dist(x)))
  expect_equal(sdnd(x), stats::sd(dist(x)))
  expect_equal(cvnd(x), stats::sd(dist(x)) / mean(dist(x)))
  expect_equal(mnnd(x, relative = TRUE), 1.6 / 8)
  expect_equal(mnd(x, relative = TRUE), mean(dist(x)) / 8)
  expect_equal(trait_range(x), 8)
  m <- x - mean(x)
  expect_equal(trait_kurtosis(x), mean(m^4) / mean(m^2)^2 - 3)
})

test_that("the minimum spanning tree agrees with a brute-force tree", {
  prim <- function(d) {
    d <- as.matrix(d)
    inside <- 1
    edges <- numeric(0)
    while (length(inside) < nrow(d)) {
      sub <- d[inside, -inside, drop = FALSE]
      k <- which(sub == min(sub), arr.ind = TRUE)[1, ]
      edges <- c(edges, min(sub))
      inside <- c(inside, setdiff(seq_len(nrow(d)), inside)[k[2]])
    }
    edges
  }
  set.seed(1)
  for (i in 1:20) {
    m <- matrix(round(rnorm(30), 1), 15)
    expect_equal(mst_length(m), sum(prim(dist(m))))
    expect_equal(mst_length(m[, 1]), sum(prim(dist(m[, 1]))))
    expect_equal(mst_ratio(m), min(prim(dist(m))) / max(prim(dist(m))))
  }
  expect_equal(mst_length(c(3, 1, 2), relative = TRUE), 1)
})

test_that("matrices, data frames and dist objects give the same distances", {
  set.seed(2)
  m <- matrix(rnorm(40), 20)
  d <- nn_distances(m)
  full <- as.matrix(dist(m))
  diag(full) <- Inf
  expect_equal(d, apply(full, 1, min))
  expect_equal(mnnd(as.data.frame(m)), mean(d))
  expect_equal(mnnd(dist(m)), mean(d))
  expect_equal(mnnd(m, relative = TRUE), mean(d) / max(dist(m)))
  expect_equal(mnnd(matrix(1:5, ncol = 1)), mnnd(1:5))
})

test_that("metrics handle missing values, small samples and bad input", {
  expect_true(is.na(mnnd(c(1, NA, 3))))
  expect_equal(mnnd(c(1, NA, 3), na.rm = TRUE), 2)
  expect_true(is.na(mnnd(cbind(c(1, NA, 3), 1:3))))
  expect_equal(mnnd(cbind(c(1, NA, 4), c(0, 0, 0)), na.rm = TRUE), 3)
  expect_true(is.na(mnnd(5)))
  expect_true(is.na(sdnnd(c(1, 2))))
  expect_true(is.na(trait_kurtosis(1:3)))
  expect_true(is.na(trait_kurtosis(rep(1, 5))))
  expect_true(is.na(mst_ratio(c(2, 2, 2))))
  expect_true(is.na(mst_length(c(1, NA))))
  expect_true(is.na(mst_length(5)))
  expect_true(is.na(cvnnd(c(1, 1, 1))))
  expect_equal(nn_distances(4), NA_real_)
  expect_true(all(is.na(nn_distances(c(1, NA)))))
  expect_error(trait_range(cbind(1:3, 1:3)), "single trait")
  expect_error(mnnd("a"), "numeric")
  d <- dist(c(1, 2, 4))
  d[1] <- NA
  expect_true(is.na(mnnd(d)))
  expect_error(mnnd(d, na.rm = TRUE), "missing")
  expect_equal(nn_distances(dist(5)), NA_real_, ignore_attr = TRUE)
})

test_that("metrics work as indices of community_index()", {
  set.seed(3)
  ind <- data.frame(site = rep(paste0("s", 1:4), each = 12),
                    species = sample(letters[1:3], 48, replace = TRUE),
                    a = rnorm(48), b = rnorm(48))
  res <- community_index(ind, list(MNND = mnnd, MST = mst_length, range = trait_range),
                         traits = "a", null = "regional.ind", nperm = 9)
  expect_equal(unique(res$index), c("MNND", "MST", "range"))
  s1 <- ind$a[ind$site == "s1"]
  expect_equal(res$observed[res$index == "MNND" & res$site == "s1"], mnnd(s1))
  expect_equal(res$observed[res$index == "range" & res$site == "s1"], diff(range(s1)))
  expect_equal(dim(attr(res, "null")), c(9, 12))
  mv <- community_index(ind, mnnd, traits = c("a", "b"), multivariate = TRUE, nperm = 0)
  expect_equal(mv$observed[1], mnnd(as.matrix(ind[ind$site == "s1", c("a", "b")])))
  expect_equal(mv$index[1], "mnnd")
})

test_that("audit edge cases of the neighbour metrics", {
  expect_equal(nn_distances(c(1, NA, 3)), rep(NA_real_, 3))
  expect_equal(nn_distances(data.frame(a = c(1, NA, 3), b = 1:3)), rep(NA_real_, 3))
  expect_error(trait_range(letters[1:5]), "numeric")
  expect_error(trait_kurtosis(letters[1:5]), "numeric")
  expect_equal(trait_range(matrix(c(1, 4, 2))), 3)
  expect_error(mnnd(array(1:8, c(2, 2, 2))), "vector, a matrix")
  expect_error(trait_range(array(1:8, c(2, 2, 2))))
})

test_that("metrics return NA for too few or missing values", {
  expect_true(is.na(mst_ratio(c(1, NA, 2))))
  expect_true(is.na(mst_ratio(3)))
  expect_true(is.na(trait_range(c(1, NA))))
  expect_equal(trait_range(c(1, NA, 4), na.rm = TRUE), 3)
  expect_true(is.na(mnd(c(1, NA, 4))))
  expect_true(is.na(mnd(1)))
})
