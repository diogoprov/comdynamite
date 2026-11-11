#' Reference table for approximate Bayesian computation
#'
#' Draws parameters from uniform priors, simulates one data set per draw and summarises it, as
#' the first step of inference by approximate Bayesian computation (ABC). What `setup` returns is
#' drawn once and shared by every simulation, so it should hold what is fixed by the data, such
#' as the observed tree and landscape. A parameter that spans orders of magnitude, such as
#' dispersal, is better drawn on the log scale and transformed back inside `simulate`.
#'
#' @param prior Named list with one numeric vector of length 2 per parameter, the lower and upper
#'   limits of its uniform prior.
#' @param simulate Function of `scenario`, a list with one value per parameter, and `shared`,
#'   what `setup` returned, which gives a simulation. As in [run_scenarios()].
#' @param n Number of simulations.
#' @param setup Optional function without arguments that gives what all simulations share.
#' @param summarise Function of `sim` and `shared` that gives a named numeric vector. The default
#'   is [summary_stats()], with `shared$tree` if there is one.
#' @param seed Integer that sets all the seeds.
#' @return An `mc_abc_table`: a data frame with one row per simulation, the parameters and the
#'   statistics, with the attributes `prior`, `parameters` and `statistics`.
#' @family inference
#' @seealso [abc_fit()] to estimate the parameters, [abc_recovery()] to check whether they can be
#'   recovered.
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' shared_pool <- function() {
#'   list(tree = treeshapesim::sim_tree_shape(10), landscape = landscape_implicit(8))
#' }
#' one_run <- function(scenario, shared) {
#'   simulate_metacommunity(shared$landscape, species_traits(shared$tree, scenario$lambda),
#'                          competition_matrix(shared$tree, "none"),
#'                          dispersal = 10^scenario$log_dispersal, n_steps = 30,
#'                          seed_steps = 1:3)
#' }
#' prior <- list(log_dispersal = c(-3, -0.3), lambda = c(0, 1))
#' table <- abc_table(prior, one_run, n = 20, setup = shared_pool)
#' head(table[, 1:5])
#' @export
abc_table <- function(prior, simulate, n, setup = NULL,
                      summarise = function(sim, shared) summary_stats(sim, shared$tree),
                      seed = 1) {
  check_prior(prior)
  stopifnot(is.function(simulate), length(n) == 1, n >= 1, length(seed) == 1)
  old <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, envir = globalenv()))
  cell_seed(seed)
  parameters <- as.data.frame(lapply(prior, function(limits) stats::runif(n, limits[1], limits[2])))
  table <- run_scenarios(parameters, simulate, setup, summarise, n_rep = 1, seed = seed + 1)
  table$replicate <- NULL
  structure(table, class = c("mc_abc_table", "data.frame"), prior = prior,
            parameters = names(prior), statistics = setdiff(names(table), names(prior)))
}

check_prior <- function(prior) {
  ok <- is.list(prior) && length(prior) >= 1 && !is.null(names(prior)) &&
    all(nzchar(names(prior))) && !anyDuplicated(names(prior)) &&
    all(vapply(prior, function(p) is.numeric(p) && length(p) == 2 && all(is.finite(p)) &&
                 p[1] < p[2], logical(1)))
  if (!ok) {
    stop("`prior` must be a named list with one increasing pair of finite limits per parameter.",
         call. = FALSE)
  }
  invisible(prior)
}

#' Estimate parameters by ABC with random forests
#'
#' Fits one regression random forest per parameter, with the statistics as predictors, following
#' Raynal et al. (2019) and the package abcrf. Random forests need no tolerance and no choice of
#' statistics, cope with many correlated statistics and rank their importance. [predict()] gives
#' the posterior median and a credible interval from quantile regression forests (Meinshausen
#' 2006). Statistics with more than 5% of missing or infinite values and constant statistics are
#' dropped, and then simulations with any missing statistic.
#'
#' Check with [abc_recovery()] that a parameter can be recovered from simulated data before
#' trusting its estimate from real data. A good estimate does not validate the model either.
#'
#' @param table An `mc_abc_table` from [abc_table()], or a data frame with the parameters and
#'   the statistics.
#' @param parameters Names of the parameters to estimate. By default those of `table`.
#' @param statistics Names of the statistics to use. By default all the others.
#' @param ntree Number of trees per forest.
#' @param ncores Number of threads for the forests.
#' @return An `mc_abc_fit`, a list with `forests` (one per parameter), `training`, `parameters`,
#'   `statistics` and `prior`.
#' @references
#' Meinshausen, N. (2006). Quantile regression forests. *Journal of Machine Learning Research*,
#' 7, 983–999.
#'
#' Raynal, L., Marin, J.-M., Pudlo, P., Ribatet, M., Robert, C. P., & Estoup, A. (2019). ABC
#' random forests for Bayesian parameter inference. *Bioinformatics*, 35, 1720–1728.
#' \doi{10.1093/bioinformatics/bty867}
#' @family inference
#' @examplesIf requireNamespace("abcrf", quietly = TRUE)
#' set.seed(1)
#' table <- data.frame(a = runif(300), b = runif(300))
#' table$s1 <- table$a + rnorm(300, 0, 0.05)
#' table$s2 <- table$b * table$a + rnorm(300, 0, 0.05)
#' fit <- abc_fit(table, parameters = c("a", "b"), ntree = 100)
#' predict(fit, c(s1 = 0.5, s2 = 0.2))
#' @export
abc_fit <- function(table, parameters = NULL, statistics = NULL, ntree = 500, ncores = 1) {
  if (!requireNamespace("abcrf", quietly = TRUE)) {
    stop("abc_fit() needs the package abcrf.", call. = FALSE) # nocov
  }
  stopifnot(is.data.frame(table))
  parameters <- parameters %||% attr(table, "parameters")
  if (is.null(parameters)) stop("Give `parameters` for a plain data frame.", call. = FALSE)
  statistics <- statistics %||% attr(table, "statistics") %||% setdiff(names(table), parameters)
  missing <- setdiff(c(parameters, statistics), names(table))
  if (length(missing)) {
    stop("Not in `table`: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  statistics <- setdiff(statistics, parameters)
  values <- as.matrix(as.data.frame(table)[, statistics, drop = FALSE])
  usable <- colMeans(is.finite(values)) > 0.95 &
    apply(values, 2, function(v) isTRUE(stats::sd(v[is.finite(v)]) > 0))
  if (!any(usable)) stop("No statistic varies among the simulations.", call. = FALSE)
  statistics <- statistics[usable]
  training <- as.data.frame(table)[, c(parameters, statistics), drop = FALSE]
  training <- training[stats::complete.cases(training) &
                         rowSums(!is.finite(as.matrix(training))) == 0, , drop = FALSE]
  if (nrow(training) < 10) stop("Fewer than 10 complete simulations.", call. = FALSE)
  forests <- lapply(stats::setNames(parameters, parameters), function(p) {
    data <- data.frame(.parameter = training[[p]], training[statistics])
    abcrf::regAbcrf(.parameter ~ ., data, ntree = ntree, paral = ncores > 1, ncores = ncores)
  })
  structure(list(forests = forests, training = training, parameters = parameters,
                 statistics = statistics, prior = attr(table, "prior")),
            class = "mc_abc_fit")
}

#' @describeIn abc_fit Posterior median and credible interval for each observation and parameter.
#' @param object An `mc_abc_fit`.
#' @param newdata Named numeric vector of observed statistics, as from [summary_stats()], or a
#'   data frame with one observation per row.
#' @param level Probability of the credible interval.
#' @param ... Not used.
#' @return For `predict()`, a data frame with `observation`, `parameter`, `estimate` (posterior
#'   median), `low`, `high` and `mean`.
#' @export
predict.mc_abc_fit <- function(object, newdata, level = 0.9, ncores = 1, ...) {
  stopifnot(length(level) == 1, level > 0, level < 1)
  if (!is.data.frame(newdata)) newdata <- as.data.frame(as.list(newdata))
  missing <- setdiff(object$statistics, names(newdata))
  if (length(missing)) {
    stop("Statistics missing from `newdata`: ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }
  observed <- newdata[, object$statistics, drop = FALSE]
  if (any(!is.finite(as.matrix(observed)))) {
    stop("`newdata` has missing or infinite statistics.", call. = FALSE)
  }
  probs <- c((1 - level) / 2, (1 + level) / 2)
  out <- lapply(object$parameters, function(p) {
    training <- data.frame(.parameter = object$training[[p]], object$training[object$statistics])
    post <- stats::predict(object$forests[[p]], observed, training, quantiles = probs,
                           paral = ncores > 1, ncores = ncores)
    data.frame(observation = seq_len(nrow(observed)), parameter = p,
               estimate = as.vector(post$med), low = post$quantiles[, 1],
               high = post$quantiles[, 2], mean = post$expectation)
  })
  out <- do.call(rbind, out)
  out[order(out$observation), , drop = FALSE] |> `rownames<-`(NULL)
}

#' @describeIn abc_fit Importance of each statistic for each parameter, as the decrease in
#'   impurity scaled to sum to one.
#' @param x An `mc_abc_fit`.
#' @export
abc_importance <- function(x) {
  stopifnot(inherits(x, "mc_abc_fit"))
  out <- lapply(x$parameters, function(p) {
    imp <- x$forests[[p]]$model.rf$variable.importance
    data.frame(parameter = p, statistic = names(imp), importance = imp / sum(imp))
  })
  out <- do.call(rbind, out)
  out <- out[order(out$parameter, -out$importance), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @export
print.mc_abc_fit <- function(x, ...) {
  cat("ABC random forests for", length(x$parameters), "parameters (",
      paste(x$parameters, collapse = ", "), ")\n")
  cat(nrow(x$training), "simulations and", length(x$statistics), "statistics\n")
  invisible(x)
}

#' Can the parameters be recovered by ABC?
#'
#' Leaves `n_test` simulations out of the reference table as pseudo-observations with known
#' parameters, fits [abc_fit()] on the rest and compares the estimates with the truth. A
#' parameter that is not recovered from simulated data should not be inferred from real data
#' (Pontarp et al. 2019).
#'
#' @inheritParams abc_fit
#' @param n_test Number of pseudo-observations.
#' @param level Probability of the credible interval.
#' @param seed Integer for the choice of the pseudo-observations.
#' @return A data frame with one row per parameter: `n`, `r2` (1 minus the squared error over the
#'   variance of the truth), `bias` and `width` (mean width of the interval) as fractions of the
#'   prior range, and `coverage`, the fraction of intervals that contain the truth, with its
#'   standard error. The attribute `estimates` holds every estimate. A coverage near `level`
#'   says the intervals are calibrated; R2 says how much the data narrow the prior.
#' @references
#' Pontarp, M., Brännström, Å., & Petchey, O. L. (2019). Inferring community assembly processes
#' from macroscopic patterns using dynamic eco-evolutionary models and Approximate Bayesian
#' Computation (ABC). *Methods in Ecology and Evolution*, 10, 450–460.
#' \doi{10.1111/2041-210X.13129}
#' @family inference
#' @examplesIf requireNamespace("abcrf", quietly = TRUE)
#' set.seed(1)
#' table <- data.frame(a = runif(400), b = runif(400))
#' table$s1 <- table$a + rnorm(400, 0, 0.05)
#' table$s2 <- rnorm(400)
#' abc_recovery(table, parameters = c("a", "b"), n_test = 50, ntree = 100)
#' @export
abc_recovery <- function(table, parameters = NULL, statistics = NULL, n_test = 200,
                         level = 0.9, ntree = 500, ncores = 1, seed = 1) {
  stopifnot(is.data.frame(table), length(n_test) == 1, n_test >= 2, n_test < nrow(table) - 10)
  parameters <- parameters %||% attr(table, "parameters")
  if (is.null(parameters)) stop("Give `parameters` for a plain data frame.", call. = FALSE)
  old <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, envir = globalenv()))
  cell_seed(seed)
  test <- sort(sample.int(nrow(table), n_test))
  train <- table[-test, , drop = FALSE]
  attributes(train)[c("prior", "parameters", "statistics")] <-
    attributes(table)[c("prior", "parameters", "statistics")]
  fit <- abc_fit(train, parameters, statistics, ntree = ntree, ncores = ncores)
  pseudo <- as.data.frame(table)[test, , drop = FALSE]
  complete <- rowSums(!is.finite(as.matrix(pseudo[, fit$statistics, drop = FALSE]))) == 0
  pseudo <- pseudo[complete, , drop = FALSE]
  estimates <- stats::predict(fit, pseudo, level = level, ncores = ncores)
  estimates$truth <- unlist(lapply(seq_len(nrow(pseudo)), function(i) pseudo[i, parameters]))
  prior <- attr(table, "prior")
  out <- do.call(rbind, lapply(parameters, function(p) {
    d <- estimates[estimates$parameter == p, ]
    range <- if (!is.null(prior[[p]])) diff(prior[[p]]) else diff(range(table[[p]]))
    covered <- d$truth >= d$low & d$truth <= d$high
    data.frame(parameter = p, n = nrow(d),
               r2 = 1 - sum((d$truth - d$estimate)^2) / sum((d$truth - mean(d$truth))^2),
               bias = mean(d$estimate - d$truth) / range,
               coverage = mean(covered),
               coverage_se = sqrt(mean(covered) * (1 - mean(covered)) / nrow(d)),
               width = mean(d$high - d$low) / range)
  }))
  attr(out, "estimates") <- estimates
  out
}
