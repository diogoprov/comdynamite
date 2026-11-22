# rOpenSci statistical standards tested in this file (see the srr tags in R/)
#' @srrstats {EA2.0, EA6.0a, EA6.0b, EA6.0c} Tested in this file, as described in
#'   R/srr-stats-eda.R.
NULL

tree <- ape::read.tree(text = "((t1:1,t2:1):2,(t3:2,(t4:1,t5:1):1):1);")
comm <- matrix(c(3, 0, 1, 2, 0,
                 1, 2, 0, 0, 4,
                 0, 1, 3, 1, 1,
                 2, 2, 0, 3, 0,
                 1, 0, 2, 0, 3,
                 0, 3, 1, 1, 2), 6, byrow = TRUE,
               dimnames = list(paste0("s", 1:6), paste0("t", 1:5)))
traits <- data.frame(size = c(5, 4, 3, 2, 1), row.names = paste0("t", 5:1))
env <- data.frame(x = c(1, 3, 2, 6, 5, 8), row.names = paste0("s", 1:6))

test_that("mc_data() matches names once and keeps everything aligned", {
  d <- mc_data(comm, tree, traits, env[6:1, , drop = FALSE], cbind(1:6, 0), letters[c(1, 1, 2, 2, 3, 3)])
  expect_s3_class(d, "mc_data")
  expect_identical(rownames(d$env), rownames(comm))
  expect_identical(d$env$x, env$x)
  expect_identical(rownames(d$traits), colnames(comm))
  expect_equal(d$traits$size, 1:5)
  expect_identical(sort(d$tree$tip.label), colnames(d$comm))
  expect_identical(unname(d$areas), letters[c(1, 1, 2, 2, 3, 3)])
  expect_identical(rownames(d$coords), rownames(comm))
  expect_output(print(d), "6 sites x 5 species")
  small <- d[2:4, c("t1", "t2", "t5")]
  expect_identical(dim(small$comm), c(3L, 3L))
  expect_identical(small$env$x, env$x[2:4])
  expect_identical(sort(small$tree$tip.label), c("t1", "t2", "t5"))
  expect_identical(rownames(small$traits), c("t1", "t2", "t5"))
  expect_identical(unname(small$areas), c("a", "b", "b"))
  expect_identical(nrow(small$coords), 3L)
  expect_error(d[, 1], "two species")
  expect_identical(as_mc_data(d), d)
})

test_that("species in only some parts are dropped with a warning", {
  extra <- cbind(comm, t9 = 1)
  d <- expect_warning_value(mc_data(extra, tree), "t9")
  expect_identical(colnames(d$comm), paste0("t", 1:5))
  d <- expect_warning_value(mc_data(comm, tree, traits[1:4, , drop = FALSE]), "t1")
  expect_identical(ape::Ntip(d$tree), 4L)
  unnamed <- comm
  rownames(unnamed) <- NULL
  expect_identical(rownames(mc_data(unnamed)$comm), paste0("site", 1:6))
  one_trait <- mc_data(comm, traits = c(t1 = 1, t2 = 2, t3 = 3, t4 = 4, t5 = 5))
  expect_identical(colnames(one_trait$traits), "trait")
})

test_that("wrong inputs stop", {
  expect_error(mc_data(unname(comm)), "species names")
  repeated <- comm
  rownames(repeated)[2] <- "s1"
  expect_error(mc_data(repeated), "repeated site names")
  expect_error(mc_data(comm, tree = "tree"), "phylo")
  expect_error(mc_data(comm, traits = 1:5), "species names")
  expect_error(mc_data(comm, traits = matrix(1:5)), "species names in its rows")
  expect_error(mc_data(comm, env = 1:3), "6 sites")
  expect_error(mc_data(comm, env = c(s1 = 1, s2 = 2, s3 = 3, s4 = 4, s5 = 5, s9 = 6)), "s6")
  expect_error(mc_data(comm, coords = data.frame(x = letters[1:6], y = 1:6)), "not numeric")
  expect_error(mc_data(comm, coords = matrix(letters[1:12], 6)), "numbers")
  expect_error(mc_data(comm, series = comm), "array")
  expect_error(mc_data(comm, series = array(1, c(6, 5, 2))), "species")
  bad <- array(1, c(5, 5, 2), dimnames = list(NULL, colnames(comm), NULL))
  expect_error(mc_data(comm, series = bad), "sites")
})

test_that("series, long tables and summaries", {
  series <- array(c(comm, comm * 2), c(6, 5, 2),
                  dimnames = list(rownames(comm), colnames(comm), c("step5", "step10")))
  d <- mc_data(comm, tree, traits, env, series = series)
  long <- as.data.frame(d)
  expect_named(long, c("site", "species", "step", "abundance"))
  expect_setequal(unique(long$step), c(5, 10))
  expect_identical(sum(long$abundance), 3 * sum(comm))
  expect_true(all(is.na(as.data.frame(mc_data(comm))$step)))
  expect_identical(dim(d[1:3, 1:2]$series), c(3L, 2L, 2L))
  listed <- structure(list(abundance = comm, env = list(a = 1:6, b = 6:1)), class = "mc_sim")
  expect_named(as_mc_data(listed)$env, c("a", "b"))
  unnamed <- series
  dimnames(unnamed)[[3]] <- NULL
  expect_identical(unique(as.data.frame(mc_data(comm, series = unnamed))$step), 1:2)
  set.seed(1)
  out <- summary(d, n_null = 9)
  set.seed(1)
  expected <- summary_stats(comm, tree, stats::setNames(traits$size, rownames(traits)), env$x,
                            n_null = 9)[names(out)]
  over_time <- c("cv_population", "cv_site", "cv_regional", "turnover")
  expected <- replace(expected, over_time, summary_stats(series)[over_time])
  expect_identical(out, expected)
  two <- mc_data(comm, env = data.frame(x = 1:6, y = 6:1))
  expect_error(summary_stats(two), "2 numeric columns")
  expect_true(is.na(summary_stats(mc_data(comm))[["trait_env"]]))
})

test_that("analyses take their arguments from an mc_data object", {
  d <- mc_data(comm, tree, traits, env, coords = cbind(1:6, 0))
  set.seed(2)
  a <- pcps_test(d, n_perm = 5)
  set.seed(2)
  expect_identical(a, pcps_test(comm, tree, env, n_perm = 5))
  set.seed(3)
  expect_identical(pcc_test(d, n_perm = 5)$tests, {
    set.seed(3)
    pcc_test(comm, tree, env, n_perm = 5)$tests
  })
  set.seed(4)
  expect_identical(phylo_beta_test(d, n_perm = 5)$tests, {
    set.seed(4)
    phylo_beta_test(comm, tree, env, n_perm = 5)$tests
  })
  expect_identical(phylo_beta(d), phylo_beta(comm, tree))
  set.seed(5)
  expect_identical(node_test(d, n_perm = 5)$nodes, {
    set.seed(5)
    nodes <- node_contrasts(d$tree)
    node_test(comm, nodes, env, n_perm = 5)$nodes
  })
  set.seed(6)
  expect_identical(phylo_vector(d, n_perm = 5, n_sim = 5)$tests, {
    set.seed(6)
    phylo_vector(comm, tree, env$x, n_perm = 5, n_sim = 5)$tests
  })
  space <- cbind(m1 = c(1, -1, 1, -1, 1, -1), m2 = c(1, 1, 0, 0, -1, -1))
  set.seed(7)
  expect_identical(trait_partition(d, space = space, n_perm = 5)$fractions, {
    set.seed(7)
    trait_partition(comm, d$traits, env, space, n_perm = 5)$fractions
  })
  expect_error(pcps_test(mc_data(comm, tree)), "no `env`")
  expect_error(phylo_vector(mc_data(comm, tree, env = data.frame(x = 1:6, y = 6:1))),
               "2 numeric columns")
})

test_that("simulations keep their tree and convert", {
  set.seed(8)
  traits <- species_traits(tree)
  expect_identical(attr(traits, "tree"), tree)
  expect_null(attr(species_traits(4), "tree"))
  land <- landscape_grid(side = 3)
  sim <- simulate_metacommunity(land, traits, competition_matrix(tree, "phylo"), n_steps = 10,
                                seed_steps = 1:3, record_steps = c(5, 10))
  expect_identical(sim$tree, tree)
  d <- as_mc_data(sim)
  expect_identical(d$env$env, land$env)
  expect_identical(unname(d$coords), unname(as.matrix(land$coords)))
  expect_identical(dim(d$series), c(9L, 5L, 2L))
  expect_identical(colnames(d$traits), c("optimum", "sigma", "r_max"))
  static <- simulate_static(traits, seq(0, 100, length.out = 4))
  expect_identical(as_mc_data(static)$tree, tree)
  expect_identical(nrow(as_mc_data(static)$env), 4L)
  logistic <- simulate_logistic(tree, logistic_gradients(1, n_sites = 7))
  expect_identical(logistic$tree, tree)
  expect_null(as_mc_data(logistic)$traits)
  expect_identical(nrow(as_mc_data(logistic)$env), 7L)
  expect_s3_class(as_mc_data(comm, tree = tree), "mc_data")
})
