# Accessor methods (coef, confint, vcov, nobs, fitted, residuals, formula) for the regression
# results of the package

#' Coefficients, intervals and fitted values of the test of drift
#'
#' Accessors for the result of [drift_test()], which regresses the logarithm of the variance
#' among group means on the logarithm of the eigenvalues of `W`, one point per principal
#' component.
#'
#' @param object Result of [drift_test()].
#' @param parm Names of the coefficients (`"intercept"`, `"slope"`); all by default.
#' @param level Confidence level, one number between 0 and 1, both excluded.
#' @param ... Not used.
#' @return `coef()`: a named vector with `intercept` and `slope`. `confint()`: a matrix with
#'   one row per coefficient and the lower and upper limits, from the \eqn{t} distribution
#'   with `n_traits - 2` degrees of freedom. `vcov()`: their covariance matrix. `nobs()`: the
#'   number of points of the regression, the number of traits. `fitted()` and `residuals()`:
#'   one value per principal component, on the scale of `log(between)`.
#' @seealso [drift_test()].
#' @family trait variation among individuals
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(300 * 4), 300, 4) %*% diag(c(3, 2, 1.5, 1))
#' group <- rep(1:30, each = 10)
#' x <- x + (matrix(rnorm(30 * 4), 30, 4) %*% diag(c(3, 2, 1.5, 1)))[group, ]
#' res <- drift_test(x, group)
#' coef(res)
#' confint(res)
#' nobs(res)
#' @exportS3Method stats::coef
coef.drift_test <- function(object, ...) {
  c(intercept = object$intercept, slope = object$slope)
}

#' @rdname coef.drift_test
#' @exportS3Method stats::confint
confint.drift_test <- function(object, parm = NULL, level = 0.95, ...) {
  check_number(level, "level", 0, 1, open = TRUE)
  out <- stats::confint(object$fit, level = level)
  rownames(out) <- c("intercept", "slope")
  if (!is.null(parm)) out <- out[pick_parm(parm, rownames(out)), , drop = FALSE]
  out
}

#' @rdname coef.drift_test
#' @exportS3Method stats::vcov
vcov.drift_test <- function(object, ...) {
  out <- stats::vcov(object$fit)
  dimnames(out) <- list(c("intercept", "slope"), c("intercept", "slope"))
  out
}

#' @rdname coef.drift_test
#' @exportS3Method stats::nobs
nobs.drift_test <- function(object, ...) object$n_traits

#' @rdname coef.drift_test
#' @exportS3Method stats::fitted
fitted.drift_test <- function(object, ...) unname(stats::fitted(object$fit))

#' @rdname coef.drift_test
#' @exportS3Method stats::residuals
residuals.drift_test <- function(object, ...) unname(stats::residuals(object$fit))

#' Coefficients, intervals and fitted values of population shifts
#'
#' Accessors for the result of [population_shifts()].
#'
#' @param object Result of [population_shifts()].
#' @param parm Names of the coefficients (the slopes `local_adaptation...` and
#'   `displacement`); all by default.
#' @param level Confidence level. The intervals are those computed by [population_shifts()]
#'   with its `level`; another level stops with an error, since the bootstrap must be run
#'   again.
#' @param ... Not used.
#' @return `coef()`: a named vector with the slopes on each environmental variable and the
#'   strength of displacement `c`. `confint()`: a matrix with one row per coefficient and the
#'   lower and upper limits of the basic bootstrap intervals (needs `nboot > 0`). `vcov()`:
#'   the covariance of the coefficients among the bootstrap simulations from the fitted model
#'   (needs `nboot > 0`). `nobs()`: the number of populations that entered as responses.
#'   `fitted()` and `residuals()`: one value per half-sample mean of a population, centred
#'   within species, in the order of the attribute `model` of `object`, at the chosen width.
#' @seealso [population_shifts()].
#' @family trait variation among individuals
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(20)
#' sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
#'                               competition_matrix(tree, "none"), n_steps = 100,
#'                               seed_steps = 1:5)
#' ind <- sample_individuals(sim, n = 20, local_adaptation = 0.5, displacement = 2)
#' fit <- population_shifts(ind, width = 0.5, nboot = 19)
#' coef(fit)
#' confint(fit)
#' nobs(fit)
#' @exportS3Method stats::coef
coef.population_shifts <- function(object, ...) {
  names <- shifts_coef_names(object)
  stats::setNames(as.numeric(unlist(as.data.frame(object)[1, names])), names)
}

#' @rdname coef.population_shifts
#' @exportS3Method stats::confint
confint.population_shifts <- function(object, parm = NULL, level = attr(object, "level"),
                                      ...) {
  stored <- attr(object, "level")
  if (is.null(stored)) stop("intervals need `nboot > 0` in population_shifts()", call. = FALSE)
  check_number(level, "level", 0, 1, open = TRUE)
  if (abs(level - stored) > 1e-12) {
    stop("the intervals were computed with `level = ", stored,
         "`; run population_shifts() again with the level wanted", call. = FALSE)
  }
  names <- shifts_coef_names(object)
  d <- as.data.frame(object)
  out <- cbind(vapply(names, function(v) d[[paste0(v, "_low")]], numeric(1)),
               vapply(names, function(v) d[[paste0(v, "_high")]], numeric(1)))
  dimnames(out) <- list(names, paste(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                                            trim = TRUE), "%"))
  if (!is.null(parm)) out <- out[pick_parm(parm, rownames(out)), , drop = FALSE]
  out
}

#' @rdname coef.population_shifts
#' @exportS3Method stats::vcov
vcov.population_shifts <- function(object, ...) {
  draws <- attr(object, "bootstrap")
  if (is.null(draws)) stop("the covariance needs `nboot > 0` in population_shifts()",
                           call. = FALSE)
  names <- shifts_coef_names(object)
  stats::cov(draws[draws$set == "fitted", names, drop = FALSE], use = "complete.obs")
}

#' @rdname coef.population_shifts
#' @exportS3Method stats::nobs
nobs.population_shifts <- function(object, ...) as.data.frame(object)$n_populations[1]

#' @rdname coef.population_shifts
#' @exportS3Method stats::fitted
fitted.population_shifts <- function(object, ...) attr(object, "model")$fitted

#' @rdname coef.population_shifts
#' @exportS3Method stats::residuals
residuals.population_shifts <- function(object, ...) attr(object, "model")$residual

#' the slope columns and `displacement`, in the order of the table
#' @noRd
shifts_coef_names <- function(object) {
  d <- as.data.frame(object)
  slopes <- grep("^local_adaptation", names(d), value = TRUE)
  c(slopes[!grepl("_(low|high)$", slopes)], "displacement")
}

#' rows of a coefficient table chosen by name or position
#' @noRd
pick_parm <- function(parm, available) {
  if (is.numeric(parm)) parm <- available[parm]
  if (!is.character(parm) || anyNA(parm) || !all(parm %in% available)) {
    stop("`parm` must name coefficients among ", toString(available),
         call. = FALSE)
  }
  parm
}

#' Model, fitted values and residuals of the decomposition of community trait means
#'
#' Accessors and print method for the result of [trait_flex()].
#'
#' @param x,object Result of [trait_flex()].
#' @param digits Number of significant digits printed. One whole number, 1 or more.
#' @param ... Not used.
#' @return `formula()`: the one-sided formula of the environmental model. `fitted()` and
#'   `residuals()`: the column of the same name of `object$model`, one value per site, trait
#'   and component (`turnover`, `intraspecific`, `total`), in the units of the trait; without
#'   a model the fitted value is the mean over sites. `nobs()`: the number of sites with
#'   values of all variables of the model. `print()` returns `x` invisibly.
#' @seealso [trait_flex()].
#' @family trait variation among individuals
#' @examples
#' set.seed(1)
#' ind <- data.frame(site = rep(paste0("s", 1:6), each = 20),
#'                   species = sample(c("a", "b", "c"), 120, replace = TRUE))
#' env <- data.frame(site = paste0("s", 1:6), moisture = 1:6)
#' ind$trait <- c(a = 10, b = 20, c = 30)[ind$species] +
#'   env$moisture[match(ind$site, env$site)] + rnorm(120)
#' res <- trait_flex(ind, formula = ~ moisture, env = env)
#' res
#' formula(res)
#' head(residuals(res))
#' @export
print.trait_flex <- function(x, digits = 3, ...) {
  check_count(digits, "digits")
  dec <- x$decomposition
  traits <- unique(dec$trait)
  cat("Decomposition of community trait means (Lep\u0161 et al. 2011):",
      length(unique(x$sites$site)), "sites,", length(traits),
      if (length(traits) == 1) "trait\n" else "traits\n")
  cat("Model:", deparse(x$formula), "\n")
  shown <- dec
  for (v in c("explained", "residual", "total", "p")) {
    shown[[v]] <- signif(shown[[v]], digits)
  }
  print(shown, row.names = FALSE)
  invisible(x)
}

#' @rdname print.trait_flex
#' @exportS3Method stats::formula
formula.trait_flex <- function(x, ...) x$formula

#' @rdname print.trait_flex
#' @exportS3Method stats::fitted
fitted.trait_flex <- function(object, ...) object$model$fitted

#' @rdname print.trait_flex
#' @exportS3Method stats::residuals
residuals.trait_flex <- function(object, ...) object$model$residual

#' @rdname print.trait_flex
#' @exportS3Method stats::nobs
nobs.trait_flex <- function(object, ...) {
  m <- object$model
  first <- m$trait == m$trait[1] & m$component == "total"
  sum(!is.na(m$fitted[first]))
}

#' Coefficients and numbers of observations of the trait-gradient analysis
#'
#' Accessors for the results of [trait_gradient()] and [trait_gradient_phylo()].
#'
#' @param object Result of [trait_gradient()] or [trait_gradient_phylo()].
#' @param ... Not used.
#' @return For [trait_gradient()], `coef()` gives a matrix with one row per species and the
#'   columns `beta`, `alpha`, `intercept` and `slope` (the last two `NA` for species with too
#'   few sites), and `nobs()` the number of populations (site by species). For
#'   [trait_gradient_phylo()], `coef()` gives a matrix with one row per response and model
#'   (named `response:model`) and the columns `mu`, `sigma2` and `a`, and `nobs()` the number
#'   of species in the fits.
#' @seealso [trait_gradient()], [trait_gradient_phylo()].
#' @family trait variation among individuals
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(20)
#' sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
#'                               competition_matrix(tree, "none"), n_steps = 100,
#'                               seed_steps = 1:5)
#' tg <- trait_gradient(sample_individuals(sim, n = 3, intraspecific_sd = 1))
#' head(coef(tg))
#' nobs(tg)
#' @exportS3Method stats::coef
coef.mc_trait_gradient <- function(object, ...) {
  s <- object$species
  out <- as.matrix(s[c("beta", "alpha", "intercept", "slope")])
  rownames(out) <- s$species
  out
}

#' @rdname coef.mc_trait_gradient
#' @exportS3Method stats::nobs
nobs.mc_trait_gradient <- function(object, ...) nrow(object$populations)

#' @rdname coef.mc_trait_gradient
#' @exportS3Method stats::coef
coef.mc_trait_gradient_phylo <- function(object, ...) {
  m <- object$models
  out <- as.matrix(m[c("mu", "sigma2", "a")])
  rownames(out) <- paste(m$response, m$model, sep = ":")
  out
}

#' @rdname coef.mc_trait_gradient
#' @exportS3Method stats::nobs
nobs.mc_trait_gradient_phylo <- function(object, ...) nrow(object$data)

#' Numbers of observations of the node test and of ABC
#'
#' @param object Result of [node_test()] or [abc_fit()].
#' @param ... Not used.
#' @return For [node_test()], the number of sites at the root, which holds all the sites
#'   analysed (each node uses the sites where its clade occurs, given in the column `n_sites`
#'   of `nodes`). For [abc_fit()], the number of complete simulations in the training table.
#' @seealso [node_test()], [abc_fit()].
#' @family inference
#' @exportS3Method stats::nobs
nobs.mc_node_test <- function(object, ...) as.integer(max(object$nodes$n_sites, na.rm = TRUE))

#' @rdname nobs.mc_node_test
#' @exportS3Method stats::nobs
nobs.mc_abc_fit <- function(object, ...) nrow(object$training)
