test_that("competition matrices: symmetric, scaled diagonal, phylo ordered by relatedness", {
  set.seed(1)
  tr <- ape::rcoal(30)
  for (type in c("none", "random", "phylo")) {
    a <- competition_matrix(tr, type)
    expect_true(isSymmetric(unname(a)))
    expect_equal(unname(diag(a)), rep(0.05, 30))
  }
  set.seed(3)
  ar <- competition_matrix(tr, "random")
  set.seed(3)
  ap <- competition_matrix(tr, "phylo")
  up <- upper.tri(ar)
  expect_equal(sort(ar[up]), sort(ap[up]))
  d <- ape::cophenetic.phylo(tr)[tr$tip.label, tr$tip.label]
  s <- (1 - d / max(d))
  rng <- tapply(ap[up], s[up], range)
  lo <- vapply(rng, `[`, 0, 1)
  hi <- vapply(rng, `[`, 0, 2)
  expect_true(all(lo[-1] >= hi[-length(hi)]))
  expect_error(competition_matrix(5, "phylo"))
})

test_that("grid landscapes: torus permutations keep distances, dispersal rows sum to 1", {
  set.seed(2)
  ls <- landscape_grid(side = 6)
  p <- torus_permutations(6)
  expect_equal(nrow(p), 6^2 * 8 - 1)
  expect_true(all(apply(p, 1, function(q) all(sort(q) == seq_len(36)))))
  expect_true(all(apply(p, 1, function(q) isTRUE(all.equal(ls$dist[q, q], ls$dist)))))
  expect_equal(unname(rowSums(ls$dispersal)), rep(1, 36))
  expect_equal(unname(diag(ls$dispersal)), rep(0, 36))
  expect_true(all(ls$env >= 0 & ls$env <= 100))
  flat <- landscape_grid(side = 6, torus = FALSE)
  expect_gt(max(flat$dist), max(ls$dist))
  expect_output(print(ls), "36 sites")
})

test_that("traits have the requested ranges and need a tree for signal", {
  set.seed(6)
  tr <- ape::rcoal(25)
  x <- species_traits(tr, lambda = 1, sigma_signal = TRUE, rmax_signal = TRUE)
  expect_equal(rownames(x), tr$tip.label)
  expect_equal(range(x$optimum), c(0, 100))
  expect_true(all(x$sigma > 0))
  expect_true(all(x$r_max >= 2 & x$r_max <= 8))
  u <- species_traits(tr, uniform_optima = TRUE)$optimum
  expect_equal(sort(u), seq(0, 100, length.out = 25))
  expect_error(species_traits(10, sigma_signal = TRUE))
  expect_equal(nrow(species_traits(10)), 10)
})

test_that("growth is maximal at the optimum and has SD sigma", {
  tr <- data.frame(optimum = 50, sigma = 10, r_max = 4, row.names = "a")
  g <- growth_rates(c(50, 60), tr)
  expect_equal(unname(g[1, 1]), 4)
  expect_equal(unname(g[2, 1]), 4 * exp(-0.5))
})

test_that("env_over_time has the requested autocorrelation", {
  set.seed(7)
  e <- env_over_time(0, 20000, rho = 0.8, sd = 2)
  expect_equal(stats::acf(e[1, ], plot = FALSE)$acf[2], 0.8, tolerance = 0.05)
  expect_equal(stats::sd(e[1, ]), 2, tolerance = 0.1)
})
