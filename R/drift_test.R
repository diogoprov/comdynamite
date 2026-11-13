#' Test whether divergence among groups is proportional to variation within them
#'
#' Tests the prediction of neutral evolution that the variance among group means is
#' proportional to the variance within groups (Lande 1979): under genetic drift alone the
#' covariance matrix of the means, \eqn{B}{B}, is expected to be \eqn{(t / N_e) G}{(t / Ne)
#' G}, where \eqn{G}{G} is the additive genetic covariance matrix, \eqn{t}{t} the number of
#' generations and \eqn{N_e}{Ne} the effective population size. Following Ackermann & Cheverud
#' (2004), the group means are projected on the principal components of the within-group
#' matrix \eqn{W}{W}, and the logarithm of the variance among groups on each component is
#' regressed on the logarithm of its eigenvalue. Under drift the slope is 1; a slope above 1
#' means more divergence along the main axes of variation than drift predicts (for example,
#' selection along them), and a slope below 1 means more divergence along the minor axes.
#'
#' By default \eqn{W}{W} is the pooled phenotypic covariance within groups, which stands in for
#' \eqn{G}{G} when it is unknown. This replacement assumes that \eqn{G}{G} and \eqn{P}{P} are
#' proportional; Prôa et al. (2013) found acceptable type I error rates with it only when the
#' matrices are similar and heritabilities are high, or when \eqn{t / N_e}{t / Ne} is small. A
#' known \eqn{G}{G} (for example, from [genetic_matrices()] or the `genetic_var` of a simulation)
#' can be given in `W` instead. Individuals should be grouped by population when `W` is
#' estimated from them: a `W` pooled within species whose populations diverged (for example by
#' local adaptation along a gradient) carries that divergence and can hide selection among
#' species.
#'
#' The variance of the observed means includes their sampling error, \eqn{P / n}{P / n} for groups
#' of \eqn{n}{n} individuals. Unless \eqn{P}{P} is proportional to `W`, this error flattens the
#' slope (most with small groups and little divergence), so by default (`correct = TRUE`) the
#' expected sampling variance on each component, \eqn{v_k^\top P v_k}{v_k' P v_k} times the mean of
#' \eqn{1 / n}{1 / n} over groups, with \eqn{P}{P} the pooled covariance within groups, is
#' subtracted from the variance among means before the regression. The regression treats each component as an independent
#' observation, so with few traits the test has little power, and the slope has \eqn{q - 2}{q -
#' 2} degrees of freedom for \eqn{q}{q} traits.
#'
#' @param x A data frame or a numeric matrix with one row per individual (for example, from
#'   [sample_individuals()]).
#' @param group The groups (populations, sites or species): a vector with one value per row of
#'   `x`, or the name of a column of `x`.
#' @param traits With a data frame, the names of the trait columns. `NULL` (the default) takes
#'   every numeric column except `group` (or the column equal to it); with the output of
#'   [sample_individuals()], which also has numeric columns of abundance, environment and means,
#'   name the traits.
#' @param W `NULL` (the default) for the pooled covariance within groups, or a covariance matrix
#'   with one row and column per trait (for example, a known \eqn{G}{G}).
#' @param correct Whether to subtract the sampling variance of the group means (needs groups
#'   with more than one individual).
#' @param level The confidence level of the interval for the slope.
#' @return An object of class `drift_test`: a list with `slope`, its standard error `se`, the
#'   interval `conf_int`, the `statistic` (\eqn{t}{t}) and `p_value` of the test of slope 1, the
#'   `intercept` (the logarithm of the proportionality constant, an estimate of \eqn{\log(t /
#'   N_e)}{log(t / Ne)} when `W` is \eqn{G}{G} and `correct = TRUE`), `df`, `n_groups`, `n_traits`, the matrix `W`
#'   and a data frame `components` with the `eigenvalue`, the expected `sampling` variance of
#'   the means (0 without correction) and the `between` variance (corrected) of each principal
#'   component. It has `print()` and `plot()` methods.
#' @references
#' Ackermann, R. R. & Cheverud, J. M. (2004). Detecting genetic drift versus selection in human
#' evolution. *Proceedings of the National Academy of Sciences*, 101, 17946--17951.
#' \doi{10.1073/pnas.0405919102}
#'
#' Lande, R. (1979). Quantitative genetic analysis of multivariate evolution, applied to brain:body
#' size allometry. *Evolution*, 33, 402--416. \doi{10.1111/j.1558-5646.1979.tb04694.x}
#'
#' Prôa, M., O'Higgins, P. & Monteiro, L. R. (2013). Type I error rates for testing genetic drift
#' with phenotypic covariance matrices: a simulation study. *Evolution*, 67, 185--195.
#' \doi{10.1111/j.1558-5646.2012.01746.x}
#' @family trait variation among individuals
#' @examples
#' set.seed(1)
#' G <- genetic_matrices(1, 4, correlation = 0.5)$genetic_var[[1]]
#' # 30 populations that diverged by drift: means drawn from N(0, G t / Ne)
#' means <- matrix(rnorm(30 * 4), 30) %*% chol(G * 0.5)
#' ind <- means[rep(1:30, each = 20), ] + matrix(rnorm(600 * 4), 600) %*% chol(G)
#' test <- drift_test(ind, group = rep(1:30, each = 20))
#' test
#' plot(test)
#' @export
drift_test <- function(x, group, traits = NULL, W = NULL, correct = TRUE, level = 0.95) {
  if (is.data.frame(x)) {
    if (is.character(group) && length(group) == 1 && group %in% names(x)) {
      gname <- group
      group <- x[[group]]
    } else {
      gname <- NULL
    }
    if (is.null(traits)) {
      numeric_cols <- names(x)[vapply(x, is.numeric, logical(1))]
      same <- vapply(numeric_cols, function(v) {
        length(group) == nrow(x) && isTRUE(all.equal(x[[v]], group, check.attributes = FALSE))
      }, logical(1))
      traits <- setdiff(numeric_cols[!same], gname)
    }
    if (!is.character(traits) || !all(traits %in% names(x))) {
      stop("`traits` must name columns of `x`", call. = FALSE)
    }
    x <- as.matrix(x[traits])
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a data frame or a numeric matrix", call. = FALSE)
  }
  if (is.null(colnames(x))) colnames(x) <- paste0("trait", seq_len(ncol(x)))
  q <- ncol(x)
  if (q < 3) stop("the test needs at least three traits", call. = FALSE)
  if (length(group) != nrow(x)) {
    stop("`group` must have one value per row of `x` or name a column", call. = FALSE)
  }
  if (!is.logical(correct) || length(correct) != 1 || is.na(correct)) {
    stop("`correct` must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1 || is.na(level) || level <= 0 || level >= 1) {
    stop("`level` must be a number between 0 and 1", call. = FALSE)
  }
  keep <- stats::complete.cases(x) & !is.na(group)
  x <- x[keep, , drop = FALSE]
  group <- factor(group[keep])
  sizes <- as.vector(table(group))
  means <- rowsum(x, group) / sizes
  if (nrow(means) < 3) stop("the test needs at least three groups", call. = FALSE)
  centred <- x - means[as.integer(group), , drop = FALSE]
  df_within <- nrow(x) - nlevels(group)
  P <- if (df_within > 0) crossprod(centred) / df_within else NULL
  if (correct && is.null(P)) {
    stop("the sampling correction needs groups with more than one individual; ",
         "use `correct = FALSE`", call. = FALSE)
  }
  if (is.null(W)) {
    if (df_within < q) {
      stop("too few individuals within groups to estimate `W`; give `W`", call. = FALSE)
    }
    W <- P
  } else {
    W <- as.matrix(if (is.data.frame(W)) data.matrix(W) else W)
    if (!is.numeric(W) || !identical(dim(W), c(q, q)) || !isSymmetric(unname(W))) {
      stop("`W` must be a symmetric matrix with one row and column per trait", call. = FALSE)
    }
    if (!is.null(rownames(W)) && !is.null(colnames(W))) {
      if (!setequal(rownames(W), colnames(x)) || !setequal(colnames(W), colnames(x))) {
        stop("the names of `W` must be those of the traits", call. = FALSE)
      }
      W <- W[colnames(x), colnames(x)]
    }
  }
  dimnames(W) <- list(colnames(x), colnames(x))
  eig <- eigen(W, symmetric = TRUE)
  if (any(eig$values <= sqrt(.Machine$double.eps) * max(abs(eig$values)))) {
    stop("`W` must be positive definite", call. = FALSE)
  }
  scores <- means %*% eig$vectors
  between <- apply(scores, 2, stats::var)
  if (any(between <= 0)) stop("the group means do not vary on every component", call. = FALSE)
  sampling <- if (correct) {
    colSums(eig$vectors * (P %*% eig$vectors)) * mean(1 / sizes)
  } else {
    rep(0, q)
  }
  between <- between - sampling
  if (any(between <= 0)) {
    stop("on some components the variance among group means is no larger than their sampling ",
         "error; measure more individuals per group or use `correct = FALSE`", call. = FALSE)
  }
  components <- data.frame(component = seq_len(q), eigenvalue = eig$values, sampling = sampling,
                           between = between)
  fit <- stats::lm(log(between) ~ log(eigenvalue), data = components)
  coef <- summary(fit)$coefficients
  slope <- coef[2, 1]
  se <- coef[2, 2]
  df <- q - 2
  statistic <- if (se > 0) (slope - 1) / se else if (slope == 1) 0 else Inf * sign(slope - 1)
  crit <- stats::qt(1 - (1 - level) / 2, df)
  structure(list(slope = slope, se = se, conf_int = slope + c(-1, 1) * crit * se,
                 level = level, statistic = statistic,
                 p_value = 2 * stats::pt(-abs(statistic), df), intercept = coef[1, 1], df = df,
                 n_groups = nrow(means), n_traits = q, W = W, components = components),
            class = "drift_test")
}

#' @export
print.drift_test <- function(x, ...) {
  cat("Test of proportional divergence (Ackermann & Cheverud 2004)\n")
  cat(sprintf("%d groups, %d traits\n", x$n_groups, x$n_traits))
  cat(sprintf("slope of log(between variance) on log(eigenvalue): %.3f (SE %.3f)\n",
              x$slope, x$se))
  cat(sprintf("%g%% interval: %.3f to %.3f\n", 100 * x$level, x$conf_int[1], x$conf_int[2]))
  cat(sprintf("test of slope = 1 (drift): t = %.3f, df = %d, p = %.3g\n", x$statistic, x$df,
              x$p_value))
  invisible(x)
}

#' @export
plot.drift_test <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot() for drift_test() needs the package ggplot2", call. = FALSE) # nocov
  }
  d <- x$components
  .data <- ggplot2::.data
  mid <- c(mean(log(d$eigenvalue)), mean(log(d$between)))
  ggplot2::ggplot(d, ggplot2::aes(log(.data$eigenvalue), log(.data$between))) +
    ggplot2::geom_abline(intercept = mid[2] - mid[1], slope = 1, linetype = "dashed",
                         colour = "grey50") +
    ggplot2::geom_abline(intercept = x$intercept, slope = x$slope, colour = "#0072B2") +
    ggplot2::geom_point(colour = "#0072B2", size = 2) +
    ggplot2::labs(x = "log(eigenvalue of W)", y = "log(variance among groups)",
                  subtitle = sprintf("slope %.2f (dashed line: slope 1, drift)", x$slope)) +
    ggplot2::theme_minimal()
}
