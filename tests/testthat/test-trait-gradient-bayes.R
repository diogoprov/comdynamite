tgb_sim <- function(weights = "occurrence") {
  set.seed(1)
  tree <- treeshapesim::sim_tree_shape(12)
  sim <- simulate_metacommunity(landscape_implicit(10), species_traits(tree, lambda = 1),
                                competition_matrix(tree, "none"), n_steps = 100,
                                seed_steps = 1:5)
  ind <- sample_individuals(sim, n = 3, intraspecific_sd = 1)
  list(tree = tree, tg = trait_gradient(ind, abundance = "abundance",
                                        species_weights = weights))
}

# posterior draws with the variables the Stan model gives back
fake_draws <- function(n_species, n = 200) {
  set.seed(9)
  d <- data.frame(lambda_beta = runif(n, 0.7, 0.9), lambda_alpha = 0.5, lambda_mean = 0.3,
                  cor_phylo_beta_mean = 0.6, cor_phylo_beta_alpha = -0.2,
                  var_phylo_beta = 4, var_species_beta = 1, var_phylo_alpha = 2,
                  var_species_alpha = 2, var_phylo_mean = 3, var_species_mean = 1)
  for (i in seq_len(n_species)) {
    d[[paste0("beta[", i, "]")]] <- i
    d[[paste0("mean_trait[", i, "]")]] <- i + 1
    d[[paste0("alpha[", i, "]")]] <- 1
  }
  d
}

test_that("the data of the model are standardised and indexed by species", {
  d <- tgb_sim()
  m <- tg_bayes_data(d$tg, d$tree)
  sd <- m$stan_data
  expect_true(all(sd$w_t == 1)) # occurrence weights
  expect_equal(mean(sd$y_t), 0, tolerance = 1e-12)
  expect_equal(stats::sd(sd$y_t), 1, tolerance = 1e-12)
  expect_equal(sd$N_t, sum(m$data$measured))
  expect_equal(sd$N_p, sum(!is.na(m$data$site_mean)))
  expect_identical(m$species[sd$sp_t], m$data$species[m$data$measured])
  expect_equal(sd$L_A %*% t(sd$L_A), m$A, ignore_attr = TRUE)
  expect_equal(unname(diag(m$A)), rep(1, nrow(m$A)))
  expect_equal(as.vector(sd$scale * sd$y_t + sd$center), m$data$trait[m$data$measured])
  expect_equal(m$kish, 1)
  pruned <- ape::drop.tip(d$tree, m$species[1])
  expect_message(tg_bayes_data(d$tg, pruned), "left out")
  expect_error(tg_bayes_data(d$tg, ape::keep.tip(d$tree, m$species[1:3])), "fewer than 4")
  expect_warning(tg_bayes_data(d$tg, d$tree, ou_alpha = 1e4), "cannot be told apart")
  expect_lt(tg_bayes_data(tgb_sim("abundance")$tg, d$tree)$kish, 1)
  flat <- d$tg
  flat$populations$trait <- 3
  expect_error(tg_bayes_data(flat, d$tree), "does not vary")
})

test_that("a singular correlation matrix is refused", {
  d <- tgb_sim()
  twin <- ape::keep.tip(d$tree, unique(d$tg$populations$species))
  twin$edge.length[twin$edge[, 2] <= ape::Ntip(twin)] <- 0
  expect_error(tg_bayes_data(d$tg, twin), "singular")
})

test_that("the OU matrix tends to Brownian motion and shrinks distant correlations", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1.5,d:1.5):0.5);")
  sp <- c("a", "b", "c", "d")
  bm <- phylo_correlation(tree, sp)
  expect_equal(bm, stats::cov2cor(ape::vcv(tree)[sp, sp]))
  expect_equal(phylo_correlation(tree, sp, 1e-6), bm, tolerance = 1e-5)
  strong <- phylo_correlation(tree, sp, 3)
  expect_lt(strong["a", "b"], bm["a", "b"])
  expect_equal(unname(diag(strong)), rep(1, 4))
  expect_error(phylo_correlation(tree, sp, -1), "`ou_alpha`")
  expect_error(phylo_correlation(tree, sp, c(1, 2)), "`ou_alpha`")
})

test_that("the summary reads lambda, variances and species values from the draws", {
  sp <- c("x", "y", "z", "w")
  s <- tg_bayes_summary(fake_draws(4), sp)
  expect_equal(s$signal$quantity, c("lambda_beta", "lambda_alpha", "lambda_mean",
                                    "cor_phylo_beta_mean", "cor_phylo_beta_alpha"))
  expect_equal(s$signal$estimate[2:5], c(0.5, 0.3, 0.6, -0.2))
  expect_gt(s$signal$lower[1], 0.7)
  expect_equal(s$variances$component, c("phylo_beta", "species_beta", "phylo_alpha",
                                        "species_alpha", "phylo_mean", "species_mean"))
  expect_equal(s$variances$estimate[1], 4)
  expect_equal(s$species$beta, 1:4)
  expect_equal(s$species$mean, 2:5)
  expect_equal(s$species$alpha, rep(1, 4))
  expect_error(tg_bayes_summary(fake_draws(4)[, -2], sp), "no parameter")
  alone <- tg_bayes_summary(fake_draws(4), sp, no_beta = "y")
  expect_true(all(is.na(alone$species[2, 2:7])))
  expect_false(anyNA(alone$species[2, 8:10]))
  expect_false(anyNA(alone$species$beta[-2]))
})

test_that("diagnostics flag poor convergence and divergences", {
  expect_true(tg_diagnostics_table(1.005, 900, 800, 0)$ok)
  expect_false(tg_diagnostics_table(1.02, 900, 800, 0)$ok)
  expect_false(tg_diagnostics_table(1.005, 300, 800, 0)$ok)
  expect_false(tg_diagnostics_table(1.005, 900, 800, 3)$ok)
  skip_if_not_installed("posterior")
  set.seed(3)
  arr <- array(rnorm(1000 * 4 * 2), c(1000, 4, 2),
               dimnames = list(NULL, NULL, c("lambda[1]", "lp__")))
  draws <- posterior::as_draws_array(arr)
  expect_true(tg_diagnostics_from(draws, 0)$ok)
  expect_equal(tg_diagnostics_from(draws, 2)$divergent, 2L)
  arr[, 1:2, 1] <- arr[, 1:2, 1] + 3 # two chains elsewhere
  expect_gt(tg_diagnostics_from(posterior::as_draws_array(arr), 0)$max_rhat, 1.1)
})

test_that("trait_gradient_bayes() samples and builds the result", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  d <- tgb_sim()
  n_sp <- length(tg_bayes_data(d$tg, d$tree)$species)
  seen <- NULL
  fake_run <- function(draws) {
    function(data, ...) {
      seen <<- c(list(data = data), list(...))
      list(fit = "fit", draws = posterior::as_draws_array(draws), divergent = 0)
    }
  }
  local_mocked_bindings(run_stan = fake_run(fake_draws(n_sp)))
  # constant fake draws have no R-hat, so the fit does not pass
  expect_warning(res <- trait_gradient_bayes(d$tg, d$tree, chains = 2, seed = 4),
                 "not converged")
  expect_s3_class(res, "mc_trait_gradient_bayes")
  expect_equal(seen$chains, 2)
  expect_equal(seen$adapt_delta, 0.95)
  expect_equal(seen$seed, 4)
  expect_equal(seen$data$lambda_prior, c(1, 1))
  expect_equal(nrow(res$species), n_sp)
  expect_identical(res$fit, "fit")
  expect_output(print(res), "Brownian motion")
  expect_output(print(res), "CHECK")
  ou <- suppressWarnings(trait_gradient_bayes(d$tg, d$tree, ou_alpha = 0.5,
                                              lambda_prior = c(2, 2)))
  expect_equal(seen$data$lambda_prior, c(2, 2))
  expect_output(print(ou), "Ornstein-Uhlenbeck")
  expect_error(trait_gradient_bayes(d$tg, d$tree, lambda_prior = c(1, -1)), "`lambda_prior`")
  expect_error(trait_gradient_bayes(list(), d$tree), "come from trait_gradient")
  expect_error(trait_gradient_bayes(d$tg, "x"), "`tree`")
  # abundance weights warn; a posterior without the expected variables keeps the fit
  local_mocked_bindings(run_stan = fake_run(fake_draws(n_sp)[, -1]))
  w <- character()
  bad <- withCallingHandlers(trait_gradient_bayes(tgb_sim("abundance")$tg, d$tree),
                             warning = function(e) {
                               w <<- c(w, conditionMessage(e))
                               invokeRestart("muffleWarning")
                             })
  expect_true(any(grepl("effective sample", w)))
  expect_true(any(grepl("could not summarise", w)))
  expect_null(bad$signal)
  expect_output(print(bad), "could not be summarised")
})

test_that("pp_check refuses an unknown response", {
  expect_error(pp_check.mc_trait_gradient_bayes(list(), resp = "x"))
})
