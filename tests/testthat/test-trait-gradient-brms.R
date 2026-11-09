tgb_sim <- function() {
  set.seed(1)
  tree <- treeshapesim::sim_tree_shape(12)
  sim <- simulate_metacommunity(landscape_implicit(10), species_traits(tree, lambda = 1),
                                competition_matrix(tree, "none"), n_steps = 100,
                                seed_steps = 1:5)
  ind <- sample_individuals(sim, n = 3, intraspecific_sd = 1)
  list(tree = tree, tg = trait_gradient(ind, abundance = "abundance",
                                        species_weights = "occurrence"))
}

# posterior draws with the columns brms gives, from known values
fake_draws <- function(species, n = 200, sp = c(1, 1), ss = c(0.5, 0.5), rp = 0.6, rs = 0) {
  set.seed(9)
  d <- data.frame(b_sitemean_Intercept = rnorm(n, 0, 0.01), b_trait_Intercept = rnorm(n, 0, 0.01),
                  sd_phylo__sitemean_Intercept = sp[1], sd_phylo__trait_Intercept = sp[2],
                  sd_species__sitemean_Intercept = ss[1], sd_species__trait_Intercept = ss[2],
                  cor_phylo__sitemean_Intercept__trait_Intercept = rp,
                  cor_species__sitemean_Intercept__trait_Intercept = rs, check.names = FALSE)
  for (s in species) for (r in c("sitemean", "trait")) {
    d[[paste0("r_phylo__", r, "[", s, ",Intercept]")]] <- if (r == "trait") 1 else 0.5
    d[[paste0("r_species__", r, "[", s, ",Intercept]")]] <- 0
  }
  d
}

test_that("the data of the model are standardised and the matrix is a correlation", {
  skip_if_not_installed("brms")
  d <- tgb_sim()
  m <- tg_brms_data(d$tg, d$tree)
  expect_equal(mean(m$data$w), 1)
  measured <- m$data$measured
  expect_equal(mean(m$data$trait[measured]), 0, tolerance = 1e-12)
  expect_equal(stats::sd(m$data$trait[measured]), 1, tolerance = 1e-12)
  expect_equal(unname(diag(m$A)), rep(1, nrow(m$A)))
  expect_identical(rownames(m$A), m$species)
  expect_s3_class(m$formula, "mvbrmsformula")
  expect_true(all(m$data$sitemean[!m$data$has_p] == 0))
  pruned <- ape::drop.tip(d$tree, m$species[1])
  expect_message(tg_brms_data(d$tg, pruned), "left out")
  expect_error(tg_brms_data(d$tg, ape::keep.tip(d$tree, m$species[1:3])), "fewer than 4")
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

test_that("priors are explicit for every class", {
  skip_if_not_installed("brms")
  p <- tg_brms_priors()
  expect_s3_class(p, "brmsprior")
  expect_setequal(unique(p$class), c("Intercept", "sd", "sigma", "cor"))
})

test_that("the summary gives lambda, correlations and species values from the draws", {
  sp <- c("x", "y", "z", "w")
  s <- tg_brms_summary(fake_draws(sp), sp, center = 10, scale = 2)
  expect_equal(s$signal$quantity, c("lambda_beta", "lambda_alpha", "lambda_mean",
                                    "cor_phylo_beta_mean", "cor_phylo_beta_alpha"))
  # sd 1 phylo and 0.5 species: lambda = 1 / (1 + 0.25)
  expect_equal(s$signal$estimate[1], 0.8)
  expect_equal(s$signal$estimate[3], 0.8)
  # alpha: phylo variance 1 + 1 - 2 * 0.6 = 0.8, species 0.25 + 0.25 = 0.5
  expect_equal(s$signal$estimate[2], 0.8 / 1.3)
  expect_equal(s$signal$estimate[4], 0.6)
  expect_equal(s$signal$estimate[5], (0.6 - 1) / sqrt(0.8))
  expect_equal(s$variances$estimate[1], 4)
  expect_equal(s$species$beta, rep(10 + 2 * 0.5, 4), tolerance = 0.01)
  expect_equal(s$species$mean, rep(12, 4), tolerance = 0.01)
  expect_equal(s$species$alpha, rep(1, 4), tolerance = 0.01)
  expect_error(tg_brms_summary(fake_draws(sp)[, -3], sp, 0, 1), "no parameter")
})

test_that("diagnostics flag poor convergence", {
  expect_true(tg_diagnostics_table(1.005, 900, 800, 0)$ok)
  expect_false(tg_diagnostics_table(1.02, 900, 800, 0)$ok)
  expect_false(tg_diagnostics_table(1.005, 300, 800, 0)$ok)
  expect_false(tg_diagnostics_table(1.005, 900, 800, 3)$ok)
})

test_that("trait_gradient_brms() calls brm and builds the result", {
  skip_if_not_installed("brms")
  skip_if_not_installed("posterior")
  d <- tgb_sim()
  sp <- tg_brms_data(d$tg, d$tree)$species
  seen <- NULL
  local_mocked_bindings(
    run_brm = function(...) {
      seen <<- list(...)
      posterior::as_draws_df(fake_draws(sp))
    },
    tg_brms_diagnostics = function(fit) tg_diagnostics_table(1.001, 1000, 1000, 0))
  res <- trait_gradient_brms(d$tg, d$tree, chains = 2, iter = 500)
  expect_s3_class(res, "mc_trait_gradient_brms")
  expect_equal(seen$chains, 2)
  expect_equal(seen$control$adapt_delta, 0.95)
  expect_identical(seen$data2$A, res$A)
  expect_equal(nrow(res$species), length(sp))
  expect_output(print(res), "Brownian motion")
  ou <- trait_gradient_brms(d$tg, d$tree, ou_alpha = 0.5, control = list(adapt_delta = 0.9))
  expect_equal(seen$control$adapt_delta, 0.9)
  expect_output(print(ou), "Ornstein-Uhlenbeck")
  local_mocked_bindings(tg_brms_diagnostics = function(fit) tg_diagnostics_table(1.2, 10, 10, 5))
  expect_warning(bad <- trait_gradient_brms(d$tg, d$tree), "not converged")
  expect_output(print(bad), "CHECK")
  expect_error(trait_gradient_brms(list(), d$tree), "come from trait_gradient")
  expect_error(trait_gradient_brms(d$tg, "x"), "`tree`")
})

test_that("diagnostics are read from the draws and the sampler", {
  skip_if_not_installed("posterior")
  set.seed(3)
  arr <- array(rnorm(1000 * 4 * 2), c(1000, 4, 2),
               dimnames = list(NULL, NULL, c("b_trait_Intercept", "lp__")))
  draws <- posterior::as_draws_array(arr)
  np <- data.frame(Parameter = c("divergent__", "treedepth__"), Value = c(0, 5))
  d <- tg_diagnostics_from(draws, np)
  expect_true(d$ok)
  expect_lt(d$max_rhat, 1.01)
  np$Value[1] <- 2
  expect_equal(tg_diagnostics_from(draws, np)$divergent, 2L)
})

test_that("species names with spaces are found as brms writes them", {
  sp <- c("Hyla one", "Hyla two", "c", "d")
  draws <- fake_draws(gsub(" ", ".", sp))
  s <- tg_brms_summary(draws, sp, center = 0, scale = 1)
  expect_identical(s$species$species, sp)
  expect_false(anyNA(s$species$beta))
})

test_that("a species alone in all its sites has no beta", {
  sp <- c("x", "y", "z", "w")
  s <- tg_brms_summary(fake_draws(sp), sp, 0, 1, no_beta = "y")
  expect_true(all(is.na(s$species[s$species$species == "y", 2:7])))
  expect_false(anyNA(s$species[s$species$species == "y", 8:10]))
  expect_false(anyNA(s$species$beta[-2]))
})

test_that("unequal weights and a strong OU warn, and a failed summary keeps the fit", {
  skip_if_not_installed("brms")
  skip_if_not_installed("posterior")
  d <- tgb_sim()
  m <- tg_brms_data(d$tg, d$tree)
  expect_lte(m$kish, 1)
  expect_gt(m$kish, 0)
  expect_warning(tg_brms_data(d$tg, d$tree, ou_alpha = 1e4), "cannot be told apart")
  sp <- m$species
  local_mocked_bindings(
    run_brm = function(...) posterior::as_draws_df(fake_draws(sp)[, -3]),
    tg_brms_diagnostics = function(fit) tg_diagnostics_table(1.001, 1000, 1000, 0),
    tg_brms_data = function(...) {
      out <- m
      out$kish <- 0.3
      out
    })
  w <- character()
  res <- withCallingHandlers(trait_gradient_brms(d$tg, d$tree),
                             warning = function(e) {
                               w <<- c(w, conditionMessage(e))
                               invokeRestart("muffleWarning")
                             })
  expect_true(any(grepl("effective sample", w)))
  expect_true(any(grepl("could not summarise", w)))
  expect_null(res$signal)
  expect_false(is.null(res$fit))
  expect_output(print(res), "could not be summarised")
})

test_that("diagnostics include the phylogenetic shares", {
  skip_if_not_installed("posterior")
  set.seed(4)
  n <- 500
  names <- c("sd_phylo__sitemean_Intercept", "sd_species__sitemean_Intercept", "lp__")
  arr <- array(abs(rnorm(n * 4 * 3)) + 0.1, c(n, 4, 3), dimnames = list(NULL, NULL, names))
  np <- data.frame(Parameter = "divergent__", Value = 0)
  good <- tg_diagnostics_from(posterior::as_draws_array(arr), np)
  expect_lt(good$max_rhat, 1.01)
  # chains stuck at opposite ends of lambda
  arr[, 1:2, 1] <- arr[, 1:2, 1] * 0.2
  arr[, 3:4, 2] <- arr[, 3:4, 2] * 0.2
  bad <- tg_diagnostics_from(posterior::as_draws_array(arr), np)
  expect_gt(bad$max_rhat, 1.01)
})
