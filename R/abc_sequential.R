#' Sequential ABC with random forests
#'
#' Narrows a wide prior towards the region that fits one observation, by rounds of simulation
#' and ABC with random forests: the ABC sequential Monte Carlo with random forests (ABC-SMC-RF)
#' of Dinh et al. (2025). The first round is a reference table drawn from the prior, as in
#' [abc_table()], and one forest per parameter as in [abc_fit()]. Each later round draws new
#' parameters from the weighted parameters of the round before, perturbs them with a normal
#' kernel of twice their weighted variance (Beaumont et al. 2009), draws again those outside the
#' prior, simulates and fits new forests on the new simulations only. The variance of the kernel
#' is that of the parameters of the round before weighted by their weights. The weights of a
#' round are those the forests give to its simulations for the observation, times the prior and
#' divided by the density of the proposal they were drawn from, so that the last round estimates
#' the posterior under the prior and not under the proposal.
#'
#' The gain is in precision for the same number of simulations when the prior is wide and the
#' posterior narrow, since later rounds simulate where the posterior is. The result belongs to
#' one observation: a new observation needs new rounds, [abc_recovery()] and [abc_calibrate()]
#' do not apply, and the coverage of the intervals can only be measured by repeating the whole
#' procedure on pseudo-observations (`data-raw/abc_sequential.R`). As in [abc_fit()], each
#' parameter has its own forest and weights, so the proposal draws each parameter on its own
#' and the correlation among parameters is not kept. The number of rounds is fixed: compare the
#' posteriors of the last rounds, with [predict()][abc_sequential()] and its `round`, to see
#' whether more would help.
#'
#' The rounds can make the intervals too narrow. The posterior of a forest rests on few
#' simulations, its effective number (`ess` of [predict()][abc_sequential()]), and its error of
#' location follows the spread of the reference table rather than the width of the posterior.
#' As the rounds concentrate the table, this error grows relative to the interval and builds up
#' from round to round. In a toy model with a known posterior and one informative statistic
#' (`data-raw/abc_sequential.R`), about 14 effective simulations per round brought the coverage
#' of 90% intervals from 0.99 in the first round to 0.90 in the third and 0.85 in the fourth,
#' while the same rounds with exact weights kept it at 0.92; wider kernels restored the coverage
#' only by giving back the narrower intervals (`data-raw/abc_sequential_diagnostics.R`). In the
#' community model of `data-raw/abc_recovery.R`, with 117 to 161 effective simulations, the
#' coverage held. A warning is given when the last round has fewer than 50 effective
#' simulations for some parameter; there, prefer fewer rounds and compare them.
#'
#' @param prior Named list with one numeric vector of length 2 per parameter, the lower and upper
#'   limits of its uniform prior, as in [abc_table()].
#' @param simulate,setup,summarise Functions to simulate and summarise, as in [abc_table()].
#'   `setup` is called once, with the seed that [abc_table()] would use, and what it returns is
#'   shared by all rounds.
#' @param observed Named numeric vector of observed statistics, as from [summary_stats()], or a
#'   data frame with one row.
#' @param n Number of simulations per round. One whole number, 10 or more.
#' @param rounds Number of rounds, the first from the prior. One whole number, 1 or more.
#' @param table Optional reference table from the prior, as from [abc_table()] with the same
#'   `prior`, `simulate`, `setup` and `summarise`, used as the first round instead of new
#'   simulations. Its rows need not number `n`.
#' @param ntree,min_node_size,ncores Passed to [abc_fit()].
#' @param seed Integer that sets all the seeds. One finite whole number.
#' @return An `mc_abc_sequential`: a list with `rounds`, one element per round with its `table`
#'   (an `mc_abc_table`), `fit` (an `mc_abc_fit`) and `weights` (a data frame with the
#'   parameters of the simulations kept by the fit and one column of normalised weights per
#'   parameter, named `w_<parameter>`); and `parameters`, `statistics`, `prior` and `observed`.
#' @family inference
#' @references
#' Beaumont, M. A., Cornuet, J.-M., Marin, J.-M. & Robert, C. P. (2009). Adaptive approximate
#' Bayesian computation. *Biometrika* 96: 983-990. \doi{10.1093/biomet/asp052}
#'
#' Dinh, K. N., Liu, C., Xiang, Z., Liu, Z. & Tavaré, S. (2025). Approximate Bayesian
#' computation sequential Monte Carlo via random forests. *Statistics and Computing* 35: 219.
#' \doi{10.1007/s11222-025-10748-x}
#' @examplesIf requireNamespace("abcrf", quietly = TRUE)
#' # one statistic, the parameter plus noise, and a prior much wider than the posterior; with
#' # so few simulations the forests rest on a handful of them, hence the warning
#' toy <- function(scenario, shared) scenario$a + rnorm(1, 0, 0.1)
#' prior <- list(a = c(-10, 10))
#' seq_fit <- abc_sequential(prior, toy, observed = c(s = 2), n = 200, rounds = 3,
#'                           summarise = function(sim, shared) c(s = sim), ntree = 100)
#' seq_fit
#' predict(seq_fit, round = 1)
#' @export
abc_sequential <- function(prior, simulate, observed, n, rounds = 4, setup = NULL,
                           summarise = function(sim, shared) summary_stats(sim, shared$tree),
                           table = NULL, ntree = 500, min_node_size = 5, ncores = 1, seed = 1) {
  cite_method("abc_sequential")
  if (!requireNamespace("abcrf", quietly = TRUE)) {
    stop("abc_sequential() needs the package abcrf.", call. = FALSE) # nocov
  }
  check_prior(prior)
  check_function_arg(simulate, "simulate")
  check_function_arg(setup, "setup", null_ok = TRUE)
  check_function_arg(summarise, "summarise")
  check_count(n, "n")
  if (n < 10) stop("`n` must be 10 or more.", call. = FALSE)
  check_count(rounds, "rounds")
  check_seed(seed)
  observed <- sequential_observed(observed)
  if (!is.null(table)) {
    check_abc_table(table)
    missing <- setdiff(names(prior), names(table))
    if (length(missing)) stop("Not in `table`: ", toString(missing), ".", call. = FALSE)
  }
  old <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  cell_seed(seed + 1)
  shared <- if (is.null(setup)) list() else setup()
  shared_setup <- function() shared
  if (is.null(table)) {
    table <- abc_table(prior, simulate, n, setup = shared_setup, summarise = summarise,
                       seed = seed)
  }
  parameters <- names(prior)
  statistics <- attr(table, "statistics") %||% setdiff(names(table), parameters)
  out <- vector("list", rounds)
  previous <- NULL
  for (r in seq_len(rounds)) {
    if (r > 1) {
      cell_seed(seed + 2 + r)
      proposal <- sequential_proposal(previous, prior, n)
      table <- run_scenarios(proposal, simulate, shared_setup, summarise, n_rep = 1,
                             seed = seed + r * 1e7)
      table$replicate <- NULL
      table <- structure(table, class = c("mc_abc_table", "data.frame"), prior = prior,
                         parameters = parameters, statistics = statistics,
                         seed = seed + r * 1e7)
    }
    fit <- abc_fit(table, parameters = parameters, statistics = statistics, ntree = ntree,
                   ncores = ncores, min_node_size = min_node_size)
    obs <- observed_statistics(observed, fit$statistics)
    weights <- fit$training[parameters]
    for (p in parameters) {
      w <- parameter_weights(fit, p, obs)[, 1]
      if (!is.null(previous)) {
        w <- w / sequential_density(weights[[p]], previous[[p]], previous[[paste0("w_", p)]])
      }
      weights[[paste0("w_", p)]] <- w / sum(w)
    }
    out[[r]] <- list(table = table, fit = fit, weights = weights)
    previous <- weights
  }
  ess <- vapply(parameters, function(p) 1 / sum(weights[[paste0("w_", p)]]^2), numeric(1))
  if (rounds > 1 && any(ess < 50)) {
    warning("abc_sequential(): Fewer than 50 effective simulations in the last round for ",
            toString(parameters[ess < 50]), "; the intervals may be too narrow. ",
            "Compare the rounds with predict(round = ).", call. = FALSE)
  }
  structure(list(rounds = out, parameters = parameters, statistics = fit$statistics,
                 prior = prior, observed = observed),
            class = "mc_abc_sequential")
}

#' the observation as a one-row data frame of finite numbers
#' @noRd
sequential_observed <- function(observed) {
  if (!is.data.frame(observed)) {
    if (!is.numeric(observed) || is.null(names(observed))) {
      stop("`observed` must be a named numeric vector or a data frame with one row.",
           call. = FALSE)
    }
    observed <- as.data.frame(as.list(observed))
  }
  if (nrow(observed) != 1) stop("`observed` must hold one observation.", call. = FALSE)
  observed
}

#' the standard deviation of the normal kernel: the square root of twice the weighted variance
#' @noRd
kernel_sd <- function(x, w) {
  m <- sum(w * x)
  sqrt(2 * sum(w * (x - m)^2))
}

#' `n` parameter sets drawn from the weighted parameters of the round before, each parameter on
#' its own, perturbed and drawn again until inside the prior
#' @noRd
sequential_proposal <- function(previous, prior, n) {
  as.data.frame(lapply(stats::setNames(names(prior), names(prior)), function(p) {
    x <- previous[[p]]
    w <- previous[[paste0("w_", p)]]
    sd <- kernel_sd(x, w)
    limits <- prior[[p]]
    out <- rep(NA_real_, n)
    todo <- seq_len(n)
    for (attempt in 1:1000) {
      draw <- x[sample.int(length(x), length(todo), replace = TRUE, prob = w)] +
        stats::rnorm(length(todo), 0, sd)
      inside <- draw >= limits[1] & draw <= limits[2]
      out[todo[inside]] <- draw[inside]
      todo <- todo[!inside]
      if (!length(todo)) break
    }
    if (length(todo)) {
      stop("abc_sequential(): Could not draw ", p, " inside its prior.", call. = FALSE) # nocov
    }
    out
  }))
}

#' density of the proposal, a mixture of normal kernels around the weighted parameters `x`,
#' at `theta`; the prior is uniform and constant inside its limits, and the truncation to the
#' prior changes the density by a constant, so both drop out when the weights are normalised
#' @noRd
sequential_density <- function(theta, x, w) {
  sd <- kernel_sd(x, w)
  keep <- w > 0
  x <- x[keep]
  w <- w[keep]
  vapply(theta, function(t) sum(w * stats::dnorm(t, x, sd)), numeric(1))
}

#' quantiles of `x` with weights `w`, the smallest value whose cumulative weight reaches each
#' probability
#' @noRd
weighted_quantile <- function(x, w, probs) {
  o <- order(x)
  x <- x[o]
  cw <- cumsum(w[o]) / sum(w)
  vapply(probs, function(p) x[min(which(cw >= p - 1e-12))], numeric(1))
}

#' Estimates from sequential ABC with random forests
#'
#' @param object An `mc_abc_sequential` from [abc_sequential()].
#' @param level Probability of the credible interval. One number greater than 0 and less
#'   than 1.
#' @param round Round whose weighted parameters give the estimates. By default the last.
#' @param ... Not used.
#' @return A data frame with one row per parameter: `round`, `parameter`, `estimate` (weighted
#'   median), `low` and `high` (the limits of the central credible interval), `mean` (weighted
#'   mean) and `ess`, the effective number of simulations behind the weights,
#'   1 / sum(weights^2).
#' @rdname abc_sequential
#' @export
predict.mc_abc_sequential <- function(object, level = 0.9, round = length(object$rounds), ...) {
  check_number(level, "level", 0, 1, open = TRUE)
  check_count(round, "round")
  if (round > length(object$rounds)) {
    stop("`round` must be ", length(object$rounds), " or less.", call. = FALSE)
  }
  weights <- object$rounds[[round]]$weights
  probs <- c(0.5, (1 - level) / 2, (1 + level) / 2)
  out <- lapply(object$parameters, function(p) {
    x <- weights[[p]]
    w <- weights[[paste0("w_", p)]]
    q <- weighted_quantile(x, w, probs)
    data.frame(round = round, parameter = p, estimate = q[1], low = q[2], high = q[3],
               mean = sum(w * x), ess = 1 / sum(w^2))
  })
  do.call(rbind, out)
}

#' @export
print.mc_abc_sequential <- function(x, ...) {
  sizes <- vapply(x$rounds, function(r) nrow(r$weights), integer(1))
  cat("Sequential ABC with random forests:", length(x$rounds), "rounds,",
      sum(sizes), "simulations used\n")
  cat("Parameters:", toString(x$parameters), "\n")
  cat("Statistics:", length(x$statistics), "\n")
  est <- stats::predict(x)
  print(est[c("parameter", "estimate", "low", "high", "ess")], row.names = FALSE,
        digits = 3)
  ess <- do.call(rbind, lapply(seq_along(x$rounds), function(r) stats::predict(x, round = r)$ess))
  dimnames(ess) <- list(paste("round", seq_along(x$rounds)), x$parameters)
  cat("Effective simulations per round:\n")
  print(round(ess, 1))
  invisible(x)
}

#' @param x An `mc_abc_sequential` from [abc_sequential()].
#' @param bins Number of bins of the histograms. One whole number, 1 or more.
#' @rdname abc_sequential
#' @export
plot.mc_abc_sequential <- function(x, bins = 30, ...) {
  need_ggplot("plot() of abc_sequential()")
  check_count(bins, "bins")
  .data <- ggplot2::.data
  d <- do.call(rbind, lapply(seq_along(x$rounds), function(r) {
    weights <- x$rounds[[r]]$weights
    do.call(rbind, lapply(x$parameters, function(p) {
      data.frame(parameter = p, round = r, value = weights[[p]],
                 weight = weights[[paste0("w_", p)]])
    }))
  }))
  d$parameter <- factor(d$parameter, x$parameters)
  d$round <- factor(d$round)
  limits <- data.frame(parameter = factor(names(x$prior), x$parameters),
                       value = c(vapply(x$prior, `[`, numeric(1), 1),
                                 vapply(x$prior, `[`, numeric(1), 2)),
                       weight = 1, round = factor(1, levels(d$round)))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$value, weight = .data$weight,
                                  fill = .data$round)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = bins,
                            position = "identity", alpha = 0.5) +
    ggplot2::geom_blank(data = limits) +
    ggplot2::scale_fill_viridis_d(end = 0.9) +
    ggplot2::facet_wrap(~parameter, scales = "free") +
    ggplot2::labs(x = "Value", y = "Posterior density", fill = "Round") +
    ggplot2::theme_minimal()
}
