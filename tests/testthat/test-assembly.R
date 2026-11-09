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

test_that("distance changes trait space as documented", {
  set.seed(10)
  x <- cbind(a = stats::rnorm(30), b = stats::rnorm(30, sd = 10))
  rownames(x) <- paste0("s", 1:30)
  # standardized: the same as Euclidean on traits divided by their standard deviations
  z <- sweep(x, 2, apply(x, 2, stats::sd), "/")
  set.seed(1)
  a <- assemble_community(x, 10, distance = "standardized")
  set.seed(1)
  b <- assemble_community(z, 10)
  expect_identical(a[c("community", "removed")], b[c("community", "removed")])
  # mahalanobis: the same as Euclidean on whitened traits
  w <- x %*% solve(chol(stats::cov(x)))
  set.seed(2)
  a <- assemble_community(x, 10, rule = "filtering", optimum = c(0, 0),
                          distance = "mahalanobis")
  set.seed(2)
  b <- assemble_community(w, 10, rule = "filtering", optimum = c(0, 0))
  expect_identical(a[c("community", "removed")], b[c("community", "removed")])
  # a given covariance replaces the one of the pool
  set.seed(3)
  a <- assemble_community(x, 10, distance = "mahalanobis",
                          trait_cov = diag(apply(x, 2, stats::var)))
  set.seed(3)
  b <- assemble_community(x, 10, distance = "standardized")
  expect_equal(a[c("community", "removed", "metric")], b[c("community", "removed", "metric")])
})

test_that("Mahalanobis filtering keeps species inside the niche ellipse", {
  # two traits along one direction: Euclidean distance favours the short axis, Mahalanobis not
  x <- rbind(along = c(2, 2), across = c(1, -1), rbind(matrix(c(1, 1, -1, -1, 3, 3, -3, -3,
                                                                0.2, -0.2, -0.2, 0.2),
                                                              ncol = 2, byrow = TRUE)))
  rownames(x) <- c("along", "across", paste0("s", 1:6))
  s <- matrix(c(1, 0.95, 0.95, 1), 2)
  euc <- assemble_community(x, 7, rule = "filtering", optimum = c(0, 0))
  mah <- assemble_community(x, 7, rule = "filtering", optimum = c(0, 0),
                            distance = "mahalanobis", trait_cov = s)
  expect_true("across" %in% euc$community)
  expect_false("across" %in% mah$community)
})

test_that("distance checks its inputs", {
  x <- matrix(stats::rnorm(6), 3, dimnames = list(c("a", "b", "c"), NULL))
  expect_error(assemble_community(x[1:2, ], 1, distance = "mahalanobis"),
               "more species than traits")
  expect_error(assemble_community(x, 2, trait_cov = diag(2)), "only with")
  expect_error(assemble_community(x, 2, distance = "mahalanobis", trait_cov = diag(3)),
               "one row and column per trait")
  expect_error(assemble_community(x, 2, distance = "mahalanobis", trait_cov = matrix(1, 2, 2)),
               "singular")
  expect_error(assemble_community(cbind(x, 1), 2, distance = "standardized"), "vary")
})

test_that("missing or non-numeric traits are rejected", {
  x <- matrix(rnorm(20), 10, dimnames = list(paste0("s", 1:10), NULL))
  x[1, 2] <- NA
  for (rule in c("filtering", "limiting_similarity", "neutral")) {
    expect_error(assemble_community(x, 5, rule), "missing values")
  }
  y <- data.frame(species = letters[1:4], t = 1:4)
  expect_error(assemble_community(y, 2), "must be numeric")
  expect_error(assemble_community(matrix(1:6, 3, dimnames = list(letters[1:3], NULL)), 2,
                                  "filtering", optimum = 1), "one value per trait")
})

test_that("assemble_community() returns an mc_assembly with the metric", {
  set.seed(2)
  x <- matrix(rnorm(60), 30, dimnames = list(paste0("s", 1:30), NULL))
  a <- assemble_community(x, 10, distance = "mahalanobis")
  expect_s3_class(a, "mc_assembly")
  expect_equal(a$metric, stats::cov(x))
  expect_equal(assemble_community(x, 10)$metric, diag(2))
  f <- assemble_community(x, 10, "filtering", optimum = c(1, 1), distance = "standardized")
  expect_equal(f$optimum, c(1, 1))
  expect_equal(unname(diag(f$metric)), unname(apply(x, 2, var)))
  expect_output(print(a), "limiting_similarity with mahalanobis distances: 10 of 30")
})

test_that("the GAUSE removal matches a scan of the whole matrix, ties included", {
  scan <- function(x, size, min_distance = NULL) {
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
  for (s in 1:12) {
    set.seed(s)
    n <- sample(5:40, 1)
    x <- if (s %% 2) matrix(sample(0:3, 2 * n, TRUE), n) else matrix(rnorm(2 * n), n)
    rownames(x) <- paste0("s", seq_len(n))
    md <- if (s %% 3 == 0) 0.8 else NULL
    set.seed(50 + s)
    a <- scan(x, 2, md)
    set.seed(50 + s)
    expect_identical(gause_removal(x, 2, md), a)
  }
})
