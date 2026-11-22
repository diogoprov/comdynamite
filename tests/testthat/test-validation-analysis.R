# rOpenSci statistical standards tested in this file (see the srr tags in R/)
#' @srrstats {BS2.5, BS2.6} Tested in this file, as described in R/srr-stats-bayesian.R.
#' @srrstats {G2.13, G5.8c} Tested in this file, as described in R/srr-stats-general.R.
NULL

# Input validation of the analysis functions: scalar arguments, Inf and coordinates

d <- daphniids
nodes <- node_contrasts(d$tree_leibold)
tr4 <- ape::read.tree(text = "((s1:1,s2:1):1,(s3:1,s4:1):1);")
toy <- data.frame(site = c("A", "A", "B", "B", "B", "C"),
                  species = c("s1", "s2", "s1", "s3", "s2", "s3"),
                  trait = c(1, 3, 2, 6, 4, 8), abundance = c(1, 1, 2, 1, 1, 1))

test_that("metrics of trait spacing check their flags and refuse Inf", {
  expect_error(mnnd(1:3, relative = NA), "`relative`")
  expect_error(mnd(1:3, relative = "yes"), "`relative`")
  expect_error(mst_length(1:3, relative = c(TRUE, FALSE)), "`relative`")
  expect_error(mnnd(1:3, na.rm = 1), "`na.rm`")
  expect_error(trait_range(1:3, na.rm = NA), "`na.rm`")
  expect_error(mnnd(c(1, Inf, 3)), "`x` has infinite values")
  expect_error(trait_kurtosis(c(1, 2, 3, -Inf)), "`x` has infinite values")
  expect_error(mnd(stats::as.dist(matrix(c(0, Inf, Inf, 0), 2))), "infinite distances")
  expect_equal(mnnd(c(1, NA, 3), na.rm = TRUE), 2)
})

test_that("node_contrasts() and biogeography_matrix() check their scalars", {
  expect_error(node_contrasts(tr4, cherries = NA), "`cherries`")
  expect_error(node_contrasts(tr4, branch_lengths = "yes"), "`branch_lengths`")
  area <- rep(c("a", "b"), each = 6)
  expect_error(biogeography_matrix(area, min_sites = -1), "`min_sites`")
  expect_error(biogeography_matrix(area, min_sites = c(1, 2)), "`min_sites`")
  expect_silent(biogeography_matrix(area, min_sites = 0))
})

test_that("node_test() checks n_perm, cores, nodes and coords", {
  expect_error(node_test(d$comm, nodes, d$env, n_perm = 0), "`n_perm`")
  expect_error(node_test(d$comm, nodes, d$env, n_perm = NA), "`n_perm`")
  expect_error(node_test(d$comm, nodes, d$env, cores = 1.5), "`cores`")
  expect_error(node_test(d$comm, nodes, d$env, cores = "2"), "`cores`")
  bad_nodes <- nodes
  bad_nodes[1, 1] <- Inf
  expect_error(node_test(d$comm, bad_nodes, d$env, n_perm = 1), "`nodes`")
  coords <- d$coords
  coords[1, 1] <- NA
  expect_error(node_test(d$comm, nodes, d$env, coords = coords), "`coords` has missing")
  coords[1, 1] <- Inf
  expect_error(node_test(d$comm, nodes, d$env, coords = coords), "`coords` has infinite")
})

test_that("tests over the whole tree check their scalars and coordinates", {
  expect_error(phylo_eigenvectors(tr4, squared = NA), "`squared`")
  expect_error(pcc_test(d$comm, d$tree, d$env, n_perm = 0), "`n_perm`")
  expect_error(pcc_test(d$comm, d$tree, d$env, squared = "no"), "`squared`")
  expect_error(phylo_beta_test(d$comm, d$tree, d$env, n_perm = "9"), "`n_perm`")
  expect_error(phylo_beta_test(d$comm, d$tree, d$env, coords = cbind(1:3, c(1, NA, 3))),
               "`coords` has missing")
  expect_error(pcps_test(d$comm, d$tree, d$env, n_perm = 2.5), "`n_perm`")
  expect_error(pcps_test(d$comm, d$tree, d$env, axes = 0), "`axes`")
  expect_error(pcps_test(d$comm, d$tree, d$env, axes = c(1, 1)), "`axes`")
  expect_error(pcps_test(d$comm, d$tree, d$env, axes = c(1, NA)), "`axes`")
  expect_error(pcps_test(d$comm, d$tree, d$env, coords = cbind(1:3, c(1, Inf, 3))),
               "`coords` has infinite")
  expect_error(check_enough_sites(matrix(1, 3, 2), matrix(1, 3, 2)), "more occupied sites")
})

test_that("phylo_vector() and trait_partition() check their scalars", {
  expect_error(phylo_vector(d$comm, d$tree, d$env[, 1], n_perm = -1), "`n_perm`")
  expect_error(phylo_vector(d$comm, d$tree, d$env[, 1], n_sim = 0), "`n_sim`")
  expect_error(phylo_vector(d$comm, d$tree, d$env[, 1], level = 1), "`level`")
  expect_error(phylo_vector(d$comm, d$tree, d$env[, 1], level = NA), "`level`")
  comm <- matrix(1, 4, 3, dimnames = list(NULL, c("a", "b", "c")))
  traits <- c(a = 1, b = 2, c = 3)
  env <- 1:4
  space <- c(2, 1, 4, 3)
  expect_error(trait_partition(comm, traits, env, space, n_perm = 0), "`n_perm`")
  expect_error(trait_partition(comm, traits, env, space, level = 0), "`level`")
  expect_error(trait_partition(comm, c(a = 1, b = Inf, c = 3), env, space),
               "`traits` has infinite")
  expect_error(trait_partition(comm, traits, env, space, coords = cbind(1:4, NA)),
               "`coords` has missing")
})

test_that("mc_data() refuses missing or infinite coordinates", {
  comm <- matrix(1, 3, 2, dimnames = list(NULL, c("a", "b")))
  expect_error(mc_data(comm, coords = cbind(1:3, c(1, NA, 3))), "`coords` has missing")
  expect_error(mc_data(comm, coords = cbind(1:3, c(1, -Inf, 3))), "`coords` has infinite")
  expect_s3_class(mc_data(comm, coords = cbind(1:3, 1:3)), "mc_data")
})

test_that("trait_gradient() and its methods check their scalars", {
  expect_error(trait_gradient(toy, trait = c("trait", "site")), "`trait`")
  expect_error(trait_gradient(toy, site = NA_character_), "`site`")
  expect_error(trait_gradient(toy, species = 3), "`species`")
  expect_error(trait_gradient(toy, abundance = TRUE), "`abundance`")
  expect_error(trait_gradient(toy, leave_one_out = "no"), "`leave_one_out`")
  inf_trait <- toy
  inf_trait$trait[1] <- Inf
  expect_error(trait_gradient(inf_trait), "`trait` has infinite")
  inf_abundance <- toy
  inf_abundance$abundance[1] <- Inf
  expect_error(trait_gradient(inf_abundance, abundance = "abundance"), "`abundance`")
  tg <- trait_gradient(toy)
  expect_error(trait_gradient_test(tg, n_null = 9, null = "phylo", tree = tr4, phylo_k = 0.5),
               "`phylo_k`")
  expect_error(trait_gradient_test(tg, n_null = 9, null = "phylo", tree = tr4, phylo_k = Inf),
               "`phylo_k`")
  expect_error(trait_gradient_bayes(tg, tr4, ou_alpha = 0), "`ou_alpha`")
  expect_error(trait_gradient_bayes(tg, tr4, chains = 0), "`chains`")
  expect_error(trait_gradient_bayes(tg, tr4, parallel_chains = NA), "`parallel_chains`")
  expect_error(trait_gradient_bayes(tg, tr4, iter_warmup = -1), "`iter_warmup`")
  expect_error(trait_gradient_bayes(tg, tr4, iter_sampling = 0), "`iter_sampling`")
  expect_error(trait_gradient_bayes(tg, tr4, adapt_delta = 1), "`adapt_delta`")
  expect_error(trait_gradient_bayes(tg, tr4, seed = "a"), "`seed`")
  expect_error(trait_gradient_bayes(tg, tr4, lambda_prior = c(1, Inf)), "`lambda_prior`")
  expect_error(pp_check.mc_trait_gradient_bayes(list(), ndraws = 0), "`ndraws`")
})

test_that("plots check their scalar arguments", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_rates(tr4, rep(1, 6), log = NA), "`log`")
  expect_error(plot_rates(tr4, rep(1, 6), tip_labels = "yes"), "`tip_labels`")
  expect_error(plot_rates(tr4, rep(1, 6), palette = 1), "`palette`")
  expect_error(plot_rates(tr4, rep(1, 6), linewidth = 0), "`linewidth`")
  expect_error(plot_rates(tr4, c(1, 1, 1, 1, 1, Inf)), "`rates` has infinite")
  expect_error(plot_node_test(NULL, tr4, level = 1.5), "`level`")
  expect_error(plot_node_test(NULL, tr4, text_size = -1), "`text_size`")
  expect_error(plot_node_heatmap(NULL, level = NA), "`level`")
  expect_error(plot_node_heatmap(NULL, text_size = "big"), "`text_size`")
  expect_error(plot_trait_partition(NULL, level = -0.1), "`level`")
  expect_error(plot_trait_partition(NULL, text_size = c(1, 2)), "`text_size`")
  expect_error(plot_phylo_beta(NULL, NULL, tr4, axis = 0), "`axis`")
  expect_error(plot_phylo_beta(NULL, NULL, tr4, text_size = 0), "`text_size`")
  fake_pcps <- structure(list(species = matrix(0.5, 2, 3, dimnames = list(c("a", "b"), NULL)),
                              scores = matrix(1:12, 4, 3), share = c(0.5, 0.3, 0.2),
                              env = matrix(1:4)), class = "mc_pcps")
  expect_error(plot(fake_pcps, species = NA), "`species`")
  expect_error(plot(fake_pcps, species_labels = 1), "`species_labels`")
  expect_error(plot(fake_pcps, palette = NULL), "`palette`")
  expect_error(plot(fake_pcps, point_size = -2), "`point_size`")
  expect_error(plot_traits(NULL, tip_labels = NA), "`tip_labels`")
  expect_error(plot_traits(NULL, point_size = 0), "`point_size`")
  expect_error(plot_traits(NULL, linewidth = Inf), "`linewidth`")
  r <- rate_matrix(c(1, 2))
  expect_error(plot(r, level = Inf), "`level`")
  expect_error(plot(r, linewidth = -1), "`linewidth`")
  expect_error(plot_assembly(NULL, species_labels = "no"), "`species_labels`")
  expect_error(plot_assembly(NULL, point_size = NA), "`point_size`")
  expect_error(plot.mc_trait_gradient(NULL, species_labels = NA), "`species_labels`")
  expect_error(plot.mc_trait_gradient(NULL, point_size = 0), "`point_size`")
})
