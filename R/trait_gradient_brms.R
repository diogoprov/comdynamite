#' Bayesian phylogenetic model of the trait gradient, with latent beta
#'
#' Fits with brms a bivariate phylogenetic mixed model in which the position of each species
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
#' (\eqn{A} the phylogenetic correlation matrix under Brownian motion), \eqn{(v^p, v^t)} are
#' species effects without phylogeny with covariance \eqn{\Sigma_v \otimes I}, and the residuals
#' are independent. So \eqn{\beta_i = \mu_p + u^p_i + v^p_i}, \eqn{\bar t_i = \mu_t + u^t_i +
#' v^t_i} and \eqn{\alpha_i = \bar t_i - \beta_i}, computed for each draw of the posterior.
#'
#' With `ou_alpha`, \eqn{A} is the correlation matrix of an Ornstein-Uhlenbeck process with
#' strength \eqn{a} = `ou_alpha` and a fixed root, as in [trait_gradient_phylo()]:
#' \eqn{C_{ij} = e^{-a(d_i + d_j - 2 s_{ij})} (1 - e^{-2 a s_{ij}}) / (2a)}, scaled to
#' correlations. \eqn{a} is taken as known, not estimated: take it from the OU fit of
#' [trait_gradient_phylo()], or fit several values and compare them with [brms::loo()]. The same
#' \eqn{a} applies to \eqn{\beta} and \eqn{\bar t}, because their phylogenetic correlation needs
#' one matrix. As \eqn{a} goes to 0 the matrix is that of Brownian motion.
#'
#' The phylogenetic share of the variance among species is
#' \eqn{\lambda = \sigma^2_u / (\sigma^2_u + \sigma^2_v)}, between 0 (no phylogenetic signal)
#' and 1 (Brownian motion). For \eqn{\alpha}, \eqn{\sigma^2_u(\alpha) = \sigma^2_u(t) +
#' \sigma^2_u(p) - 2 \rho_u \sigma_u(t) \sigma_u(p)}, and likewise for \eqn{\sigma^2_v(\alpha)}.
#' The phylogenetic correlation of \eqn{\beta} with \eqn{\bar t} is \eqn{\rho_u}, and with
#' \eqn{\alpha} it follows from the same terms.
#'
#' Both responses are centred on the mean trait value and divided by its standard deviation, so
#' that \eqn{\alpha} keeps its meaning, and the results are returned in trait units. Weights
#' (from `species_weights` in [trait_gradient()]) are divided by their mean, so that they change
#' the relative influence of populations without changing the amount of data.
#'
#' **Priors**, on the standardised scale: `normal(0, 1)` for the intercepts, `exponential(1)`
#' for the standard deviations of species effects and residuals, and `lkj(2)` for the
#' correlations; give `prior` to change them. Check convergence (`diagnostics`: R-hat below
#' 1.01, bulk and tail effective sample sizes above 400, no divergent transitions) and the fit,
#' with [bayesplot::pp_check()] on the result, before reading the estimates.
#'
#' @section Limits:
#' The phylogenetic part is Brownian motion or an Ornstein-Uhlenbeck process with a given
#' strength, shared by both responses; [trait_gradient_phylo()] estimates the strength. Site
#' means of species that share a site are almost the same
#' number, but the model treats them as independent observations, so the posterior of
#' \eqn{\beta} is somewhat too narrow when species share many sites. Sites where a species is
#' the only species have no site mean without it and enter only its trait.
#'
#' Splitting the variance among species into a phylogenetic and a non-phylogenetic part needs
#' many species. In trials with 18 species the posterior of \eqn{\lambda} had modes near 0 and
#' 1 and the chains did not mix (R-hat 1.04--1.27), even with `adapt_delta = 0.99` and
#' tighter priors; with 41 species they converged (R-hat 1.01, bulk ESS 800) and recovered
#' \eqn{\lambda_\beta} near 1 (interval 0.99--1) for niches that evolved by Brownian motion and
#' near 0 (0--0.02) for niches without phylogenetic signal. Check `diagnostics` always.
#'
#' @param x Result of [trait_gradient()].
#' @param tree Phylogeny (`phylo`) with branch lengths; species not in it are left out.
#' @param ou_alpha `NULL` for Brownian motion, or the strength \eqn{a} of an
#'   Ornstein-Uhlenbeck process, in units of 1 / branch length.
#' @param prior Priors as a `brmsprior` object (from [brms::prior()]), on the standardised
#'   scale, or `NULL` for the defaults in Details.
#' @param chains,iter,warmup,cores,seed Passed to [brms::brm()].
#' @param ... Further arguments to [brms::brm()], such as `backend = "cmdstanr"`; `control`
#'   defaults to `list(adapt_delta = 0.95)`.
#' @return An object of class `mc_trait_gradient_brms`, a list with
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
#'     \item{fit}{the `brmsfit`.}
#'     \item{data, center, scale}{the data of the model and the standardisation.}
#'   }
#' @seealso [trait_gradient()], [trait_gradient_phylo()]
#' @family trait variation among individuals
#' @examplesIf interactive() && requireNamespace("brms", quietly = TRUE) && requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(20)
#' sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree, lambda = 1),
#'                               competition_matrix(tree, "none"), n_steps = 100,
#'                               seed_steps = 1:5)
#' ind <- sample_individuals(sim, n = 5, local_adaptation = 0.5, intraspecific_sd = 1)
#' tg <- trait_gradient(ind, abundance = "abundance")
#' bayes <- trait_gradient_brms(tg, tree, chains = 4, cores = 4)
#' bayes$diagnostics
#' bayes$signal
#' bayesplot::pp_check(bayes, resp = "sitemean")
#' @export
trait_gradient_brms <- function(x, tree, ou_alpha = NULL, prior = NULL, chains = 4,
                                iter = 2000,
                                warmup = floor(iter / 2), cores = 1, seed = NA, ...) {
  if (!inherits(x, "mc_trait_gradient")) {
    stop("`x` must come from trait_gradient()", call. = FALSE)
  }
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("trait_gradient_brms() needs the package brms: install.packages(\"brms\")", # nocov
         call. = FALSE) # nocov
  }
  check_tree(tree)
  model <- tg_brms_data(x, tree, ou_alpha)
  prior <- prior %||% tg_brms_priors()
  args <- list(...)
  if (is.null(args$control)) args$control <- list(adapt_delta = 0.95)
  fit <- do.call(run_brm, c(list(model$formula, data = model$data, data2 = list(A = model$A),
                                 prior = prior, chains = chains, iter = iter, warmup = warmup,
                                 cores = cores, seed = seed), args))
  draws <- as.data.frame(posterior::as_draws_df(fit))
  summary <- tg_brms_summary(draws, model$species, model$center, model$scale)
  diagnostics <- tg_brms_diagnostics(fit)
  if (!diagnostics$ok) {
    warning("the model has not converged or has divergent transitions; see `diagnostics`",
            call. = FALSE)
  }
  structure(c(summary, list(diagnostics = diagnostics, fit = fit, data = model$data,
                            center = model$center, scale = model$scale, A = model$A,
                            ou_alpha = ou_alpha)),
            class = "mc_trait_gradient_brms")
}

# brms::brm(), apart so that tests can replace it
run_brm <- function(...) brms::brm(...) # nocov

# data, phylogenetic correlation matrix and formula of the bivariate model
tg_brms_data <- function(x, tree, ou_alpha = NULL) {
  p <- x$populations
  keep <- p$species %in% tree$tip.label
  dropped <- setdiff(unique(p$species), tree$tip.label)
  if (length(dropped) > 0) {
    message("species not in `tree`, left out: ", paste(dropped, collapse = ", "))
  }
  p <- p[keep, , drop = FALSE]
  species <- sort(unique(p$species))
  if (length(species) < 4) stop("fewer than 4 species of `x` are in `tree`", call. = FALSE)
  has_p <- !is.na(p$site_mean)
  center <- mean(p$trait[p$measured])
  scale <- stats::sd(p$trait[p$measured])
  if (!(scale > 0)) stop("the trait does not vary", call. = FALSE)
  # rows outside a subset are not used, but brms needs a number
  data <- data.frame(sitemean = ifelse(has_p, (p$site_mean - center) / scale, 0),
                     trait = ifelse(p$measured, (p$trait - center) / scale, 0),
                     w = p$weight / mean(p$weight), has_p = has_p, measured = p$measured,
                     phylo = p$species, species = p$species, site = p$site,
                     stringsAsFactors = FALSE)
  A <- phylo_correlation(ape::keep.tip(tree, species), species, ou_alpha)
  formula <- brms::bf(sitemean | weights(w) + subset(has_p) ~
                        1 + (1 | p | gr(phylo, cov = A)) + (1 | q | species)) +
    brms::bf(trait | weights(w) + subset(measured) ~
               1 + (1 | p | gr(phylo, cov = A)) + (1 | q | species)) +
    brms::set_rescor(FALSE)
  list(data = data, A = A, formula = formula, species = species, center = center,
       scale = scale)
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

# default priors on the standardised scale
tg_brms_priors <- function() {
  c(brms::prior("normal(0, 1)", class = "Intercept", resp = "sitemean"),
    brms::prior("normal(0, 1)", class = "Intercept", resp = "trait"),
    brms::prior("exponential(1)", class = "sd", resp = "sitemean"),
    brms::prior("exponential(1)", class = "sd", resp = "trait"),
    brms::prior("exponential(1)", class = "sigma", resp = "sitemean"),
    brms::prior("exponential(1)", class = "sigma", resp = "trait"),
    brms::prior("lkj(2)", class = "cor"))
}

# signal, variances and species values from a data frame of posterior draws
tg_brms_summary <- function(draws, species, center, scale) {
  col <- function(name) {
    if (!name %in% names(draws)) stop("the posterior has no parameter ", name, call. = FALSE)
    draws[[name]]
  }
  sp_b <- col("sd_phylo__sitemean_Intercept") * scale
  sp_t <- col("sd_phylo__trait_Intercept") * scale
  ss_b <- col("sd_species__sitemean_Intercept") * scale
  ss_t <- col("sd_species__trait_Intercept") * scale
  r_p <- col("cor_phylo__sitemean_Intercept__trait_Intercept")
  r_s <- col("cor_species__sitemean_Intercept__trait_Intercept")
  vp_beta <- sp_b^2
  vs_beta <- ss_b^2
  vp_mean <- sp_t^2
  vs_mean <- ss_t^2
  vp_alpha <- sp_t^2 + sp_b^2 - 2 * r_p * sp_t * sp_b
  vs_alpha <- ss_t^2 + ss_b^2 - 2 * r_s * ss_t * ss_b
  share <- function(u, v) u / (u + v)
  quantities <- list(
    lambda_beta = share(vp_beta, vs_beta),
    lambda_alpha = share(vp_alpha, vs_alpha),
    lambda_mean = share(vp_mean, vs_mean),
    cor_phylo_beta_mean = r_p,
    cor_phylo_beta_alpha = (r_p * sp_t * sp_b - sp_b^2) / sqrt(vp_beta * vp_alpha))
  interval <- function(v) stats::quantile(v, c(0.5, 0.025, 0.975), names = FALSE)
  table_of <- function(l, label) {
    out <- data.frame(names(l), t(vapply(l, interval, numeric(3))), row.names = NULL)
    names(out) <- c(label, "estimate", "lower", "upper")
    out
  }
  signal <- table_of(quantities, "quantity")
  variances <- table_of(list(phylo_beta = vp_beta, species_beta = vs_beta,
                             phylo_alpha = vp_alpha, species_alpha = vs_alpha,
                             phylo_mean = vp_mean, species_mean = vs_mean), "component")
  effect <- function(resp, sp) {
    col(paste0("b_", resp, "_Intercept")) + col(paste0("r_phylo__", resp, "[", sp, ",Intercept]")) +
      col(paste0("r_species__", resp, "[", sp, ",Intercept]"))
  }
  per_species <- lapply(species, function(sp) {
    beta <- center + scale * effect("sitemean", sp)
    mean <- center + scale * effect("trait", sp)
    alpha <- mean - beta
    v <- c(interval(beta), interval(alpha), interval(mean))
    names(v) <- c("beta", "beta_lower", "beta_upper", "alpha", "alpha_lower", "alpha_upper",
                  "mean", "mean_lower", "mean_upper")
    v
  })
  species_table <- data.frame(species = species, do.call(rbind, per_species), row.names = NULL)
  list(signal = signal, variances = variances, species = species_table)
}

# convergence: largest R-hat, smallest bulk and tail ESS of the parameters, divergences
tg_brms_diagnostics <- function(fit) {
  tg_diagnostics_from(posterior::as_draws_array(fit), brms::nuts_params(fit)) # nocov
}

tg_diagnostics_from <- function(draws, np) {
  keep <- grep("^(b_|sd_|cor_|sigma_)", posterior::variables(draws), value = TRUE)
  s <- posterior::summarise_draws(posterior::subset_draws(draws, variable = keep),
                                  "rhat", "ess_bulk", "ess_tail")
  divergent <- sum(np$Value[np$Parameter == "divergent__"])
  tg_diagnostics_table(max(s$rhat, na.rm = TRUE), min(s$ess_bulk, na.rm = TRUE),
                       min(s$ess_tail, na.rm = TRUE), divergent)
}

tg_diagnostics_table <- function(rhat, ess_bulk, ess_tail, divergent) {
  data.frame(max_rhat = rhat, min_ess_bulk = ess_bulk, min_ess_tail = ess_tail,
             divergent = as.integer(divergent),
             ok = rhat < 1.01 && ess_bulk > 400 && ess_tail > 400 && divergent == 0)
}

#' @export
print.mc_trait_gradient_brms <- function(x, ...) {
  cat("Bayesian phylogenetic model of the trait gradient:", nrow(x$species), "species,",
      if (is.null(x$ou_alpha)) "Brownian motion\n" else
        paste0("Ornstein-Uhlenbeck (a = ", format(signif(x$ou_alpha, 3)), ")\n"))
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
#' @param object Result of [trait_gradient_brms()].
#' @param resp `"sitemean"` or `"trait"`.
#' @param ... Further arguments to [brms::pp_check()].
#' @return A ggplot object from bayesplot.
#' @exportS3Method bayesplot::pp_check
pp_check.mc_trait_gradient_brms <- function(object, resp = c("sitemean", "trait"), ...) {
  resp <- match.arg(resp)
  brms::pp_check(object$fit, resp = resp, ...) # nocov
}
