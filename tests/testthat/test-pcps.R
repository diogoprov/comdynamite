sim_pcps <- function(seed = 1, n = 30, k = 10) {
  set.seed(seed)
  tree <- ape::rcoal(k)
  optima <- stats::setNames(as.vector(ape::cophenetic.phylo(tree)[, 1]), tree$tip.label)
  env <- seq(0, 1, length.out = n)
  rate <- exp(-outer(env, optima / max(optima), "-")^2 / 0.05)
  comm <- matrix(rpois(n * k, 5 * rate), n, dimnames = list(paste0("s", seq_len(n)), names(optima)))
  keep <- rowSums(comm) > 0
  list(comm = comm[keep, ], tree = tree, env = cbind(x = env[keep]))
}

test_that("fuzzy weighting and the axes follow their definitions", {
  d <- sim_pcps()
  q <- fuzzy_belonging(d$tree)
  expect_equal(unname(rowSums(q)), rep(1, nrow(q)))
  w <- d$comm[, d$tree$tip.label] / rowSums(d$comm)
  p <- w %*% q
  bray <- as.matrix(dist(p, "manhattan")) / outer(rowSums(p), rowSums(p), "+")
  g <- pcps_gower(w, q)
  ax <- pcps_axes(g)
  n <- nrow(g)
  centre <- diag(n) - 1 / n
  expect_equal(g, unname(centre %*% (-0.5 * bray) %*% centre))
  expect_equal(sum(ax$share), 1)
  expect_true(all(diff(ax$share) <= 0))
  expect_equal(colSums(ax$scores^2)[1] / sum(colSums(ax$scores^2)), ax$share[1],
               ignore_attr = TRUE)
})

test_that("statistics of pcps_test() are those of lm() and of a distance-based regression", {
  d <- sim_pcps(2)
  out <- pcps_test(d$comm, d$tree, d$env, n_perm = 9)
  expect_s3_class(out, "mc_pcps")
  expect_named(out$tests, c("test", "share", "r2", "statistic", "p_site", "p_species", "p_max"))
  expect_equal(out$tests$test, c("whole", "pcps1", "pcps2"))
  expect_equal(out$tests$statistic[2],
               unname(summary(lm(out$scores[, 1] ~ d$env))$fstatistic[1]))
  expect_equal(out$tests$r2[3], summary(lm(out$scores[, 2] ~ d$env))$r.squared)
  expect_equal(out$tests$p_max, pmax(out$tests$p_site, out$tests$p_species))
  expect_output(print(out), "mc_pcps")
  expect_identical(summary(out), out$tests)
  expect_identical(as.data.frame(out), out$tests)
})

test_that("pcps_test() matches values of the package PCPS", {
  # PCPS 1.0.9: pcps.sig() with FUN.GLM and matrix.p.sig() with FUN.ADONIS2.global
  comm <- matrix(c(4, 0, 1, 0, 2, 3, 1, 0, 0, 5, 2, 1, 0, 1, 3, 4, 2, 0, 0, 1, 5, 0, 3, 1, 1,
                   0, 2, 6, 0, 2), 6, byrow = TRUE,
                 dimnames = list(paste0("s", 1:6), paste0("t", 1:5)))
  tree <- ape::read.tree(text = "((t1:1,t2:1):2,(t3:2,(t4:1,t5:1):1):1);")
  env <- c(1, 2, 4, 7, 8, 12)
  out <- pcps_test(comm, tree, env, n_perm = 3)
  expect_equal(out$tests$statistic, c(1.1641872889, 0.3314600482, 8.9845658671), tolerance = 1e-8)
  expect_equal(out$tests$share[2:3], c(0.4909659, 0.2503570), tolerance = 1e-6)
})

test_that("site_permutations replace the shuffle of sites and allow empty sites", {
  d <- sim_pcps(3, n = 36)
  d$comm <- matrix(rpois(36 * 10, 3) + 1, 36, dimnames = list(paste0("s", 1:36), d$tree$tip.label))
  d$env <- cbind(x = seq_len(36))
  moves <- torus_permutations(6)[1:19, ]
  out <- pcps_test(d$comm, d$tree, d$env, n_perm = 5, site_permutations = moves)
  expect_equal(out$permutations$sites, 19)
  expect_equal(out$permutations$null, "restricted")
  empty <- d$comm
  empty[3, ] <- 0
  env <- d$env
  out <- pcps_test(empty, d$tree, env, n_perm = 5, site_permutations = moves)
  expect_equal(nrow(out$scores), 35)
  expect_error(pcps_test(empty, d$tree, env, n_perm = 5), "at least one individual")
})

test_that("a null without variance does not count, and wrong axes stop", {
  fit <- pcps_fit(diag(4) - 0.25, matrix(rnorm(8), 4), cbind(rep(1, 4)))
  expect_true(all(is.na(fit$statistic)))
  d <- sim_pcps(4)
  expect_error(pcps_test(d$comm, d$tree, d$env, axes = 500, n_perm = 3), "positive eigenvalues")
  expect_error(pcps_test(d$comm, d$tree, d$env, axes = 0, n_perm = 3))
  expect_equal(pcps_test(d$comm, d$tree, d$env, axes = 3, n_perm = 3)$tests$test,
               c("whole", "pcps3"))
})

test_that("torus_permutations() handles rectangles, symmetries and subsets", {
  rectangle <- as.matrix(expand.grid(x = 1:6, y = 1:4))
  full <- torus_permutations(c(6, 4))
  expect_equal(nrow(full), 4 * 24 - 1)
  torus <- function(xy) {
    dx <- abs(outer(xy[, 1], xy[, 1], "-"))
    dy <- abs(outer(xy[, 2], xy[, 2], "-"))
    pmin(dx, 6 - dx)^2 + pmin(dy, 4 - dy)^2
  }
  expect_true(all(apply(full, 1, function(p) all(torus(rectangle[p, ]) == torus(rectangle)))))
  expect_equal(nrow(torus_permutations(c(6, 4), symmetries = FALSE)), 23)
  expect_equal(nrow(torus_permutations(5, symmetries = FALSE)), 24)
  expect_equal(nrow(torus_permutations(5)), 8 * 25 - 1)
  expect_equal(nrow(torus_permutations(5, n_perm = 10)), 10)
  expect_equal(nrow(torus_permutations(3, symmetries = FALSE, n_perm = 100)), 8)
  expect_error(torus_permutations(c(2, 3, 4)))
  expect_error(torus_permutations(1))
})
