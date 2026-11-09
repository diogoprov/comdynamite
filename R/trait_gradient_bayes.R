#' Bayesian phylogenetic model of the trait gradient, with latent beta
#'
#' Fits in Stan a bivariate phylogenetic mixed model in which the position of each species
#' along the gradient (\eqn{\beta}) and its mean trait (\eqn{\bar t}) are latent species
#' effects, estimated from every site of the species, and returns the phylogenetic signal of
#' \eqn{\beta}, \eqn{\alpha} and \eqn{\bar t} with their credible intervals. Unlike
#' [trait_gradient_phylo()], which first computes \eqn{\beta} and then models it with a
#' standard error, here the uncertainty of \eqn{\beta} enters the model directly.
#'
#' @details
#' Each population (species \eqn{i} in site \eqn{j}) contributes two observations: the mean of
#' its site, \eqn{p_{ij}} (without the focal species, as chosen in [trait_gradient()]), and its
#' trait value \eqn{t_{ij}}. The model is
#' \deqn{p_{ij} = \mu_p + u^p_i + v^p_i + e^p_{ij}, \qquad
#'       t_{ij} = \mu_t + u^t_i + v^t_i + e^t_{ij},}
#' where \eqn{(u^p, u^t)} are phylogenetic effects with covariance \eqn{\Sigma_u \otimes A}
#' (\eqn{A} the phylogenetic correlation matrix), \eqn{(v^p, v^t)} are species effects without
#' phylogeny with covariance \eqn{\Sigma_v \otimes I}, and the residuals are independent. So
#' \eqn{\beta_i = \mu_p + u^p_i + v^p_i}, \eqn{\bar t_i = \mu_t + u^t_i + v^t_i} and
#' \eqn{\alpha_i = \bar t_i - \beta_i}, computed for each draw of the posterior.
#'
#' The variance among species of each response is written as a total, \eqn{\tau^2}, and a
#' phylogenetic share \eqn{\lambda}: \eqn{\sigma_u = \tau \sqrt{\lambda}} and
#' \eqn{\sigma_v = \tau \sqrt{1 - \lambda}}. \eqn{\lambda} is between 0 (no phylogenetic
#' signal) and 1 (the structure of `A` only: Brownian motion, or the OU process given), and its
#' prior is set directly, `beta(lambda_prior[1], lambda_prior[2])`, uniform by default. (A model
#' with independent priors on \eqn{\sigma_u} and \eqn{\sigma_v}, as in brms, puts an implicit
#' U-shaped prior on \eqn{\lambda}.) For \eqn{\alpha}, \eqn{\sigma^2_u(\alpha) = \sigma^2_u(t) +
#' \sigma^2_u(p) - 2 \rho_u \sigma_u(t) \sigma_u(p)}, and likewise for \eqn{\sigma^2_v(\alpha)};
#' its \eqn{\lambda} follows. The phylogenetic correlation of \eqn{\beta} with \eqn{\bar t} is
#' \eqn{\rho_u}, and with \eqn{\alpha} it follows from the same terms. Both effects are
#' non-centred, so the sampler works when either share is near 0.
#'
#' With `ou_alpha`, \eqn{A} is the correlation matrix of an Ornstein-Uhlenbeck process with
#' strength \eqn{a} = `ou_alpha` and a fixed root, as in [trait_gradient_phylo()]:
#' \eqn{C_{ij} = e^{-a(d_i + d_j - 2 s_{ij})} (1 - e^{-2 a s_{ij}}) / (2a)}, scaled to
#' correlations. \eqn{a} is taken as known, not estimated: take it from the OU fit of
#' [trait_gradient_phylo()] for \eqn{\beta} (`response == "beta"`, `model == "OU"`, column
#' `a`), which is in the same units. The same \eqn{a} applies to \eqn{\beta} and \eqn{\bar t},
#' because their phylogenetic correlation needs one matrix. As \eqn{a} goes to 0 the matrix is
#' that of Brownian motion; at large \eqn{a} it approaches the identity, and the phylogenetic
#' and species effects can no longer be told apart (the function warns). On a tree that is not
#' ultrametric, scaling the Brownian covariance to correlations is not Brownian motion, because
#' the tips have different variances.
#'
#' Both responses are centred on the mean trait value and divided by its standard deviation, so
#' that \eqn{\alpha} keeps its meaning, and the results are returned in trait units. Weights
#' (from `species_weights` in [trait_gradient()]) multiply the log-likelihood of each
#' population and are divided by their mean, so that they change the relative influence of
#' populations without changing their number. A weighted observation then counts as that many
#' full observations, so unequal weights, such as abundances, make the posterior too narrow:
#' the effective sample is \eqn{(\sum w)^2 / \sum w^2} rather than the number of populations.
#' The function warns when this is below 70% of the populations; `species_weights =
#' "occurrence"` avoids the problem.
#'
#' **Priors**, on the standardised scale: `normal(0, 1)` for the intercepts, `exponential(1)`
#' for the total standard deviations among species and for the residuals, `beta(1, 1)` for
#' both \eqn{\lambda} (change with `lambda_prior`), and `lkj(2)` for the two correlations.
#' Check convergence (`diagnostics`: R-hat below 1.01, bulk and tail effective sample sizes
#' above 400 and no divergent transitions, over the intercepts, standard deviations,
#' \eqn{\lambda}, correlations, species values and `lp__`) and the fit, with
#' [bayesplot::pp_check()] on the result, before reading the estimates.
#'
#' The model is in `system.file("stan", "trait_gradient.stan", package = "comdynamite")`. It
#' is compiled by cmdstanr the first time it is used and kept in
#' `tools::R_user_dir("comdynamite", "cache")`. The pointwise log-likelihoods are in
#' `log_lik_p` (site means) and `log_lik_t` (trait values): cross-validation, for instance to
#' compare values of `ou_alpha`, needs one response at a time, `fit$loo(variables =
#' "log_lik_p")`, because the two responses are fitted to different sets of populations.
#'
#' @section Limits:
#' The phylogenetic part is Brownian motion or an Ornstein-Uhlenbeck process with a given
#' strength, shared by both responses; [trait_gradient_phylo()] estimates the strength. Site
#' means of species that share a site are almost the same number, but the model treats them as
#' independent observations, so the posterior of \eqn{\beta} is somewhat too narrow when
#' species share many sites. The site mean of a population also contains the trait values of
#' the species it lives with, so when relatives co-occur, the phylogenetic correlation
#' \eqn{\rho_u} of \eqn{\beta} and \eqn{\bar t} is inflated by their shared values, not only by
#' niche conservatism. The residuals of the two responses are taken as independent, which
#' ignores the slope of the trait on the gradient within species; [trait_gradient()] estimates
#' those slopes. Sites where a species is the only species have no site mean without it and
#' enter only its trait; a species that is alone in all its sites has no data on \eqn{\beta},
#' and its \eqn{\beta} and \eqn{\alpha} are returned as `NA`.
#'
#' Splitting the variance among species into a phylogenetic and a non-phylogenetic part needs
#' many species, and the intervals of \eqn{\lambda} stay wide. In 60 trials with 40--45
#' species, 30 sites, occurrence weights and the default settings, where the niches evolved
#' with Pagel's \eqn{\lambda} of 0, 0.5 or 1 (20 trials each):
#' \tabular{llll}{
#'   \tab \eqn{\lambda = 0} \tab \eqn{\lambda = 0.5} \tab \eqn{\lambda = 1} \cr
#'   median of \eqn{\lambda_\beta} (range of the 20 medians) \tab 0.07 (0.05--0.31) \tab
#'     0.24 (0.04--0.73) \tab 0.94 (0.38--0.99) \cr
#'   95\% interval contains 0.5 \tab 4 \tab 17 \tab 9 \cr
#'   lower bound above 0.05 / upper bound below 0.95 \tab 0 / 20 \tab 5 / 20 \tab 16 / 1 \cr
#'   passed every diagnostic \tab 12 \tab 11 \tab 8
#' }
#' A fit took 9--30 s with four chains in parallel. One of the 60 fits had a divergent
#' transition; the fits that failed the diagnostics did so by an effective sample size below
#' 400 or an R-hat between 1.01 and 1.03, for which longer chains (`iter_sampling`) are the
#' first remedy (not tested here).
#' Intermediate signal is estimated low (interpretation: the \eqn{\lambda} of \eqn{\beta} is
#' below that of the niches, because \eqn{\beta} also carries the error of the site means),
#' and its interval rarely excludes 0: a wide interval for \eqn{\lambda_\beta} means that the
#' data do not tell. An earlier version of this model fitted with brms, with independent
#' priors on the two standard deviations, had divergent transitions in most fits and pulled
#' intermediate signal to 0 (median 0.04, interval containing 0.5 in 12 of the same 20
#' trials). In trials with 18 species the chains of that version did not mix. Use
#' [trait_gradient_phylo()] to test the signal as well, and do not read a fit that fails
#' `diagnostics`.
#'
#' @param x Result of [trait_gradient()].
#' @param tree Phylogeny (`phylo`) with branch lengths; species not in it are left out.
#' @param ou_alpha `NULL` for Brownian motion, or the strength \eqn{a} of an
#'   Ornstein-Uhlenbeck process, in units of 1 / branch length.
#' @param lambda_prior The two shape parameters of the beta prior of both \eqn{\lambda}.
#' @param chains,parallel_chains,iter_warmup,iter_sampling,adapt_delta,seed Passed to the
#'   `$sample()` method of cmdstanr.
#' @param ... Further arguments to `$sample()`, such as `max_treedepth` or `refresh`.
#' @return An object of class `mc_trait_gradient_bayes`, a list with
#'   \describe{
#'     \item{signal}{the phylogenetic share \eqn{\lambda} of \eqn{\beta}, \eqn{\alpha} and the
#'       mean, and the phylogenetic correlations of \eqn{\beta} with the mean and with
#'       \eqn{\alpha}: `estimate` (posterior median), `lower`, `upper` (95% interval).}
#'     \item{variances}{the phylogenetic and species variances of the three quantities, in
#'       squared trait units.}
#'     \item{species}{posterior median and 95% interval of \eqn{\beta}, \eqn{\alpha} and the
#'       mean of each species.}
#'     \item{diagnostics}{largest R-hat, smallest bulk and tail effective sample sizes,
#'       number of divergent transitions, and `ok`.}
#'     \item{fit}{the `CmdStanMCMC` object; its draws live in temporary files of the session,
#'       so save them with `fit$save_object()` to keep them.}
#'     \item{data, A, ou_alpha, lambda_prior}{the data of the model, the phylogenetic
#'       correlation matrix and the settings.}
#'   }
#' @seealso [trait_gradient()], [trait_gradient_phylo()]
#' @family trait variation among individuals
#' @examplesIf interactive() && requireNamespace("cmdstanr", quietly = TRUE) && requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(40)
#' sim <- simulate_metacommunity(landscape_implicit(30), species_traits(tree, lambda = 1),
#'                               competition_matrix(tree, "none"), n_steps = 100,
#'                               seed_steps = 1:5)
#' ind <- sample_individuals(sim, n = 5, local_adaptation = 0.5, intraspecific_sd = 1)
#' tg <- trait_gradient(ind, abundance = "abundance", species_weights = "occurrence")
#' bayes <- trait_gradient_bayes(tg, tree, parallel_chains = 4)
#' bayes$diagnostics
#' bayes$signal
#' bayesplot::pp_check(bayes, resp = "sitemean")
#' @export
trait_gradient_bayes <- function(x, tree, ou_alpha = NULL, lambda_prior = c(1, 1), chains = 4,
                                 parallel_chains = 1, iter_warmup = 1000,
                                 iter_sampling = 1000, adapt_delta = 0.95, seed = NULL, ...) {
  if (!inherits(x, "mc_trait_gradient")) {
    stop("`x` must come from trait_gradient()", call. = FALSE)
  }
  for (pkg in c("cmdstanr", "posterior")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("trait_gradient_bayes() needs the package ", pkg, # nocov
           " (cmdstanr from https://stan-dev.r-universe.dev)", call. = FALSE) # nocov
    }
  }
  check_tree(tree)
  if (!is.numeric(lambda_prior) || length(lambda_prior) != 2 || anyNA(lambda_prior) ||
      any(lambda_prior <= 0)) {
    stop("`lambda_prior` must be two positive numbers", call. = FALSE)
  }
  model <- tg_bayes_data(x, tree, ou_alpha, lambda_prior)
  if (model$kish < 0.7) {
    warning("the weights leave an effective sample of ", round(100 * model$kish),
            "% of the populations, and the posterior treats the weighted data as if all ",
            "populations were full observations, so its intervals are too narrow; consider ",
            "species_weights = \"occurrence\" in trait_gradient()", call. = FALSE)
  }
  run <- run_stan(model$stan_data, chains = chains, parallel_chains = parallel_chains,
                  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                  adapt_delta = adapt_delta, seed = seed, ...)
  summary <- tryCatch(
    tg_bayes_summary(posterior::as_draws_df(run$draws), model$species, model$no_beta),
    error = function(e) {
      warning("could not summarise the posterior (", conditionMessage(e), "); the fit is in ",
              "`fit`", call. = FALSE)
      list(signal = NULL, variances = NULL, species = NULL)
    })
  diagnostics <- tg_diagnostics_from(run$draws, run$divergent)
  if (!diagnostics$ok) {
    warning("the model has not converged or has divergent transitions; see `diagnostics`",
            call. = FALSE)
  }
  structure(c(summary, list(diagnostics = diagnostics, fit = run$fit, data = model$data,
                            A = model$A, ou_alpha = ou_alpha, lambda_prior = lambda_prior)),
            class = "mc_trait_gradient_bayes")
}

# the variables read back from the sampler
tg_bayes_variables <- function() {
  c("lp__", "mu", "tau", "lambda", "rho_u", "rho_v", "sigma", "lambda_beta", "lambda_alpha",
    "lambda_mean", "cor_phylo_beta_mean", "cor_phylo_beta_alpha", "var_phylo_beta",
    "var_species_beta", "var_phylo_alpha", "var_species_alpha", "var_phylo_mean",
    "var_species_mean", "beta", "alpha", "mean_trait")
}

# compiled model, kept in the user cache
tg_stan_model <- function() {
  file <- system.file("stan", "trait_gradient.stan", package = "comdynamite") # nocov start
  cache <- tools::R_user_dir("comdynamite", "cache")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  local <- file.path(cache, "trait_gradient.stan")
  if (!file.exists(local) || !identical(readLines(local), readLines(file))) {
    file.copy(file, local, overwrite = TRUE)
  }
  cmdstanr::cmdstan_model(local, quiet = TRUE) # nocov end
}

# sampling, apart so that tests can replace it: the fit, its draws and divergences
run_stan <- function(data, ...) {
  fit <- tg_stan_model()$sample(data = data, ...) # nocov start
  list(fit = fit, draws = fit$draws(variables = tg_bayes_variables()),
       divergent = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)) # nocov end
}

# data of the model, phylogenetic correlation matrix and standardisation
tg_bayes_data <- function(x, tree, ou_alpha = NULL, lambda_prior = c(1, 1)) {
  p <- x$populations
  dropped <- setdiff(unique(p$species), tree$tip.label)
  if (length(dropped) > 0) {
    message("species not in `tree`, left out: ", paste(dropped, collapse = ", "))
  }
  p <- p[p$species %in% tree$tip.label, , drop = FALSE]
  species <- sort(unique(p$species))
  if (length(species) < 4) stop("fewer than 4 species of `x` are in `tree`", call. = FALSE)
  has_p <- !is.na(p$site_mean)
  center <- mean(p$trait[p$measured])
  scale <- stats::sd(p$trait[p$measured])
  if (!(scale > 0)) stop("the trait does not vary", call. = FALSE)
  A <- phylo_correlation(ape::keep.tip(tree, species), species, ou_alpha)
  if (!is.null(ou_alpha) && max(abs(A[upper.tri(A)])) < 0.01) {
    warning("with this `ou_alpha` the phylogenetic correlations are all below 0.01, so the ",
            "phylogenetic and species effects cannot be told apart", call. = FALSE)
  }
  L_A <- tryCatch(t(chol(A)), error = function(e) {
    stop("the phylogenetic correlation matrix is singular (tips at distance 0?)",
         call. = FALSE)
  })
  w <- p$weight / mean(p$weight)
  index <- match(p$species, species)
  stan_data <- list(S = length(species), L_A = L_A,
                    N_p = sum(has_p), y_p = as.array((p$site_mean[has_p] - center) / scale),
                    sp_p = as.array(index[has_p]), w_p = as.array(w[has_p]),
                    N_t = sum(p$measured),
                    y_t = as.array((p$trait[p$measured] - center) / scale),
                    sp_t = as.array(index[p$measured]), w_t = as.array(w[p$measured]),
                    lambda_prior = lambda_prior, center = center, scale = scale)
  list(stan_data = stan_data, data = p, A = A, species = species, center = center,
       scale = scale, kish = sum(w)^2 / sum(w^2) / length(w),
       no_beta = setdiff(species, unique(p$species[has_p])))
}

# correlation matrix of the tips under Brownian motion (ou_alpha NULL) or an OU process with
# a fixed root and strength ou_alpha
phylo_correlation <- function(tree, species, ou_alpha = NULL) {
  shared <- ape::vcv(tree)[species, species]
  if (is.null(ou_alpha)) return(stats::cov2cor(shared))
  if (!is.numeric(ou_alpha) || length(ou_alpha) != 1 || is.na(ou_alpha) || !(ou_alpha > 0)) {
    stop("`ou_alpha` must be NULL or one positive number", call. = FALSE)
  }
  d <- diag(shared)
  cov <- exp(-ou_alpha * (outer(d, d, "+") - 2 * shared)) * -expm1(-2 * ou_alpha * shared) /
    (2 * ou_alpha)
  stats::cov2cor(cov)
}

# signal, variances and species values from a data frame of posterior draws
tg_bayes_summary <- function(draws, species, no_beta = character()) {
  col <- function(name) {
    if (!name %in% names(draws)) stop("the posterior has no parameter ", name, call. = FALSE)
    draws[[name]]
  }
  interval <- function(v) stats::quantile(v, c(0.5, 0.025, 0.975), names = FALSE)
  table_of <- function(names, label) {
    out <- data.frame(names, t(vapply(names, function(n) interval(col(n)), numeric(3))),
                      row.names = NULL)
    names(out) <- c(label, "estimate", "lower", "upper")
    out
  }
  signal <- table_of(c("lambda_beta", "lambda_alpha", "lambda_mean", "cor_phylo_beta_mean",
                       "cor_phylo_beta_alpha"), "quantity")
  variances <- table_of(paste0("var_", c("phylo_beta", "species_beta", "phylo_alpha",
                                         "species_alpha", "phylo_mean", "species_mean")),
                        "component")
  variances$component <- sub("^var_", "", variances$component)
  per_species <- lapply(seq_along(species), function(i) {
    v <- c(interval(col(paste0("beta[", i, "]"))), interval(col(paste0("alpha[", i, "]"))),
           interval(col(paste0("mean_trait[", i, "]"))))
    # a species always alone in its sites has no site mean: its beta is only the prior
    if (species[i] %in% no_beta) v[1:6] <- NA_real_
    v
  })
  species_table <- data.frame(species = species, do.call(rbind, per_species), row.names = NULL)
  names(species_table)[-1] <- c("beta", "beta_lower", "beta_upper", "alpha", "alpha_lower",
                                "alpha_upper", "mean", "mean_lower", "mean_upper")
  list(signal = signal, variances = variances, species = species_table)
}

# convergence over the draws read back: largest R-hat, smallest bulk and tail ESS
tg_diagnostics_from <- function(draws, divergent) {
  s <- posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail")
  tg_diagnostics_table(max(s$rhat, na.rm = TRUE), min(s$ess_bulk, na.rm = TRUE),
                       min(s$ess_tail, na.rm = TRUE), divergent)
}

tg_diagnostics_table <- function(rhat, ess_bulk, ess_tail, divergent) {
  data.frame(max_rhat = rhat, min_ess_bulk = ess_bulk, min_ess_tail = ess_tail,
             divergent = as.integer(divergent),
             ok = rhat < 1.01 && ess_bulk > 400 && ess_tail > 400 && divergent == 0)
}

#' @export
print.mc_trait_gradient_bayes <- function(x, ...) {
  cat("Bayesian phylogenetic model of the trait gradient:", nrow(x$A), "species,",
      if (is.null(x$ou_alpha)) "Brownian motion\n" else
        paste0("Ornstein-Uhlenbeck (a = ", format(signif(x$ou_alpha, 3)), ")\n"))
  if (is.null(x$signal)) {
    cat("The posterior could not be summarised; see `fit`\n")
    return(invisible(x))
  }
  d <- x$diagnostics
  cat("Convergence: max R-hat ", format(round(d$max_rhat, 3)), ", min ESS bulk ",
      round(d$min_ess_bulk), ", tail ", round(d$min_ess_tail), ", divergent transitions ",
      d$divergent, if (d$ok) " (ok)" else " (CHECK BEFORE READING)", "\n", sep = "")
  s <- x$signal
  s[-1] <- lapply(s[-1], function(v) round(v, 3))
  print(s, row.names = FALSE)
  invisible(x)
}

#' Posterior predictive check of the trait-gradient model
#'
#' Overlays the density of the observed site means or trait values on the densities of
#' replicated data drawn from the posterior of [trait_gradient_bayes()].
#'
#' @param object Result of [trait_gradient_bayes()].
#' @param resp `"sitemean"` or `"trait"`.
#' @param ndraws Number of posterior draws to plot.
#' @param ... Further arguments to [bayesplot::ppc_dens_overlay()].
#' @return A ggplot object from bayesplot.
#' @exportS3Method bayesplot::pp_check
pp_check.mc_trait_gradient_bayes <- function(object, resp = c("sitemean", "trait"),
                                             ndraws = 50, ...) {
  resp <- match.arg(resp)
  p <- object$data
  if (resp == "sitemean") {
    y <- p$site_mean[!is.na(p$site_mean)]
    variable <- "y_rep_p"
  } else {
    y <- p$trait[p$measured]
    variable <- "y_rep_t"
  }
  yrep <- posterior::as_draws_matrix(object$fit$draws(variables = variable)) # nocov start
  rows <- sample.int(nrow(yrep), min(ndraws, nrow(yrep)))
  bayesplot::ppc_dens_overlay(y, unclass(yrep)[rows, , drop = FALSE], ...) # nocov end
}
