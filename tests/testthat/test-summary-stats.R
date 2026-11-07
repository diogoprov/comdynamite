small <- matrix(c(4, 1, 0, 0,
                  0, 2, 2, 0,
                  1, 0, 0, 5,
                  0, 0, 0, 0), 4, byrow = TRUE,
                dimnames = list(paste0("s", 1:4), paste0("t", 1:4)))
tree <- ape::read.tree(text = "((t1:1,t2:1):2,(t3:2,t4:2):1);")

test_that("community statistics follow their definitions", {
  out <- summary_stats(small)
  expect_equal(out[["alpha"]], 6 / 4)
  expect_equal(out[["gamma"]], 4)
  expect_equal(out[["beta"]], 4 / 1.5)
  expect_equal(out[["occupancy"]], mean(c(2, 2, 1, 1) / 4))
  expect_equal(out[["occupancy_sd"]], sd(c(2, 2, 1, 1) / 4))
  expect_equal(out[["abundance"]], 15 / 4)
  shannon <- function(a) -sum(a / sum(a) * log(a / sum(a)))
  expect_equal(out[["evenness"]], mean(c(shannon(c(4, 1)), shannon(c(2, 2)), shannon(c(1, 5)))) / log(2))
  expect_equal(out[["evenness_regional"]], shannon(c(5, 3, 2, 5)) / log(4))
  expect_equal(out[["sad_sd"]], sd(log(c(5, 3, 2, 5))))
  expect_true(all(is.na(out[c("ses_mpd", "ses_mntd", "decay", "trait_env", "trait_sd", "cv_population",
                              "cv_site", "cv_regional", "turnover")])))
})

test_that("names do not change with the inputs, and an empty community gives no errors", {
  full <- summary_stats(small, tree, c(t1 = 1, t2 = 2, t3 = 3, t4 = 4), 1:4, cbind(1:4, 0), 5)
  expect_named(summary_stats(small), names(full))
  empty <- summary_stats(small * 0, tree, c(t1 = 1, t2 = 2, t3 = 3, t4 = 4), 1:4, cbind(1:4, 0), 5)
  expect_named(empty, names(full))
  expect_equal(empty[["alpha"]], 0)
  expect_true(is.na(empty[["beta"]]))
  expect_true(is.na(empty[["ses_mpd"]]))
  expect_true(is.na(summary_stats(small[, 1, drop = FALSE])[["occupancy_sd"]]))
  expect_true(is.na(summary_stats(small * 0)[["occupancy"]]))
})

test_that("phylogenetic statistics match picante", {
  skip_if_not_installed("picante")
  set.seed(1)
  big_tree <- ape::rcoal(12)
  comm <- matrix(rbinom(10 * 12, 1, 0.4), 10, dimnames = list(NULL, big_tree$tip.label))
  comm[1, ] <- 0
  comm[1, 1] <- 1
  distances <- ape::cophenetic.phylo(big_tree)[colnames(comm), colnames(comm)]
  expect_equal(mean_pairwise(comm, distances)[-1], picante::mpd(comm, distances)[-1])
  expect_equal(unname(mean_nearest(comm, distances))[-1], picante::mntd(comm, distances)[-1])
  ours <- summary_stats(comm, big_tree, n_null = 1999)
  theirs <- picante::ses.mpd(comm, distances, null.model = "taxa.labels", runs = 1999)
  expect_lt(abs(ours[["ses_mpd"]] - mean(theirs$mpd.obs.z, na.rm = TRUE)), 0.05)
  theirs <- picante::ses.mntd(comm, distances, null.model = "taxa.labels", runs = 1999)
  expect_lt(abs(ours[["ses_mntd"]] - mean(theirs$mntd.obs.z, na.rm = TRUE)), 0.05)
})

test_that("decay and traits follow their definitions", {
  coords <- cbind(c(0, 1, 3, 9), 0)
  out <- summary_stats(small, coords = coords)
  occupied <- small[1:3, ]
  bray <- as.matrix(dist(occupied, "manhattan")) / outer(rowSums(occupied), rowSums(occupied), "+")
  expect_equal(out[["decay"]], cor(as.vector(as.dist(bray)), as.vector(dist(coords[1:3, ]))))
  skip_if_not_installed("vegan")
  expect_equal(bray_curtis(occupied), as.vector(vegan::vegdist(occupied)))
  trait <- c(t4 = 40, t3 = 30, t2 = 20, t1 = 10)
  out <- summary_stats(small, trait = trait, env = c(1, 2, 3, 4))
  weighted <- c(sum(c(4, 1) * c(10, 20)) / 5, sum(c(2, 2) * c(20, 30)) / 4, sum(c(1, 5) * c(10, 40)) / 6)
  expect_equal(out[["trait_env"]], cor(weighted, 1:3))
  individuals <- rep(c(10, 20, 30, 40), c(5, 3, 2, 5))
  expect_equal(out[["trait_sd"]], sqrt(mean((individuals - mean(individuals))^2)))
  expect_true(is.na(summary_stats(small, trait = trait)[["trait_env"]]))
  expect_true(is.na(summary_stats(small[c(1, 4, 4), ], coords = cbind(1:3, 0))[["decay"]]))
  expect_true(is.na(summary_stats(small[c(1, 1, 1), ], coords = cbind(1:3, 0), )[["decay"]]))
})

test_that("temporal statistics follow their definitions", {
  series <- array(c(small, small * 2, small * 3), c(4, 4, 3),
                  dimnames = list(rownames(small), colnames(small), NULL))
  series[1, 1, 3] <- 0
  out <- summary_stats(series)
  expect_equal(out[["abundance"]], sum(series[, , 3]) / 4)
  populations <- matrix(series, ncol = 3)
  populations <- populations[rowSums(populations) > 0, ]
  expect_equal(out[["cv_population"]], mean(apply(populations, 1, function(v) sd(v) / mean(v))))
  sites <- apply(series, c(1, 3), sum)[1:3, ]
  expect_equal(out[["cv_site"]], mean(apply(sites, 1, function(v) sd(v) / mean(v))))
  regional <- apply(series, c(2, 3), sum)
  expect_equal(out[["cv_regional"]], mean(apply(regional, 1, function(v) sd(v) / mean(v))))
  still <- array(small, c(4, 4, 2), dimnames = list(rownames(small), colnames(small), NULL))
  expect_equal(summary_stats(still)[["turnover"]], 0)
  expect_equal(summary_stats(still)[["cv_population"]], 0)
  expect_true(is.na(summary_stats(series[, , 1, drop = FALSE])[["turnover"]]))
})

test_that("an mc_sim gives its own environment, coordinates, trait and history", {
  set.seed(2)
  traits <- species_traits(tree)
  land <- landscape_grid(side = 3)
  sim <- simulate_metacommunity(land, traits, competition_matrix(tree, "phylo"), n_steps = 12,
                                seed_steps = 1:3, record_steps = c(8, 10, 12))
  out <- summary_stats(sim, tree, n_null = 5)
  expect_false(is.na(out[["decay"]]))
  expect_false(is.na(out[["trait_env"]]))
  expect_false(is.na(out[["turnover"]]))
  same <- summary_stats(sim$abundance, trait = stats::setNames(traits$optimum, rownames(traits)),
                        env = land$env, coords = land$coords)
  expect_equal(out[c("alpha", "decay", "trait_env")], same[c("alpha", "decay", "trait_env")])
  static <- simulate_static(traits, seq(0, 100, length.out = 6))
  expect_false(is.na(summary_stats(static)[["trait_env"]]))
})

test_that("wrong inputs stop", {
  expect_error(summary_stats(small, env = 1:3), "one number per site")
  expect_error(summary_stats(small, env = letters[1:4]), "one number per site")
  expect_error(summary_stats(small, coords = cbind(1:3, 0)), "one row per site")
  expect_error(summary_stats(small, trait = 1:4), "named by species")
  expect_error(summary_stats(small, trait = c(t1 = 1)), "no value")
  expect_error(summary_stats(small, n_null = 1))
})
