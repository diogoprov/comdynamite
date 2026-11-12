#' Folds of sites for cross-validation
#'
#' Splits the sites into `k` folds, either as spatial blocks, clusters of neighbouring sites
#' found by k-means on the coordinates, or at random. With spatial autocorrelation, sites left
#' out at random have close neighbours among the sites kept, so random folds can make a model
#' look better at predicting new sites than it is; spatial blocks hold out whole regions
#' (Roberts et al. 2017).
#'
#' @param coords Matrix or data frame with one row per site and the coordinates in columns, such
#'   as `landscape$coords` or the `coords` of an [mc_data()].
#' @param k Number of folds.
#' @param method `"blocks"` for spatial blocks or `"random"`.
#' @param seed Integer that sets the seed.
#' @return Integer vector with the fold of each site.
#' @references
#' Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J., Guillera-Arroita, G., Hauenstein,
#' S., Lahoz-Monfort, J. J., Schröder, B., Thuiller, W., Warton, D. I., Wintle, B. A., Hartig, F.,
#' & Dormann, C. F. (2017). Cross-validation strategies for data with temporal, spatial,
#' hierarchical, or phylogenetic structure. *Ecography*, 40, 913–929. \doi{10.1111/ecog.02881}
#' @family inference
#' @examples
#' grid <- landscape_grid(side = 6)
#' table(spatial_folds(grid$coords, k = 4))
#' @export
spatial_folds <- function(coords, k = 5, method = c("blocks", "random"), seed = 1) {
  method <- match.arg(method)
  coords <- as.matrix(coords)
  stopifnot(is.numeric(coords), nrow(coords) >= 2, length(k) == 1, k >= 2, k <= nrow(coords))
  old <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, envir = globalenv()))
  cell_seed(seed)
  if (method == "random") {
    return(sample(rep_len(seq_len(k), nrow(coords))))
  }
  fit <- stats::kmeans(scale(coords, scale = FALSE), centers = k, nstart = 20, iter.max = 100)
  # number the blocks from the first coordinate, so the folds do not depend on the labels of
  # k-means
  order_blocks <- order(fit$centers[, 1], fit$centers[, ncol(fit$centers)])
  match(fit$cluster, order_blocks)
}

#' Statistics by fold of sites
#'
#' Turns a function that summarises a data set into one that summarises, for every fold, the
#' sites kept (`train`) and the sites left out (`test`). Use it as `summarise` in
#' [abc_table()], so that one simulation gives the statistics of every fold, and apply it to the
#' observed data before [abc_cv()]. The statistics are named `f1_train_alpha`,
#' `f1_test_alpha` and so on.
#'
#' @param folds Integer vector with the fold of each site, as from [spatial_folds()].
#' @param summarise Function of an [mc_data()] object and `shared` that gives a named numeric
#'   vector. The default is [summary_stats()].
#' @return A function of `sim`, a simulation or an [mc_data()] object, and `shared`.
#' @family inference
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' tree <- treeshapesim::sim_tree_shape(8)
#' grid <- landscape_grid(side = 4)
#' sim <- simulate_metacommunity(grid, species_traits(tree, 0.5), competition_matrix(tree, "none"),
#'                               dispersal = 0.05, n_steps = 40)
#' by_fold <- fold_summarise(spatial_folds(grid$coords, k = 2),
#'                           function(x, shared) summary_stats(x, n_null = 9)[1:3])
#' by_fold(sim, list())
#' @export
fold_summarise <- function(folds, summarise = function(x, shared) summary_stats(x)) {
  stopifnot(length(folds) >= 2, is.function(summarise))
  folds <- as.integer(folds)
  ids <- sort(unique(folds))
  if (length(ids) < 2) stop("`folds` needs at least two folds.", call. = FALSE)
  function(sim, shared) {
    data <- if (inherits(sim, "mc_data")) sim else as_mc_data(sim)
    if (nrow(data$comm) != length(folds)) {
      stop("`folds` needs one value per site.", call. = FALSE)
    }
    unlist(lapply(ids, function(f) {
      train <- summarise(data[folds != f, ], shared)
      test <- summarise(data[folds == f, ], shared)
      c(stats::setNames(train, paste0("f", f, "_train_", names(train))),
        stats::setNames(test, paste0("f", f, "_test_", names(test))))
    }))
  }
}

#' Cross-validation of an ABC model across sites
#'
#' For each fold, fits [abc_fit()] on the statistics of the sites kept, draws from the
#' approximate posterior given the observed statistics of those sites, and predicts the
#' statistics of the sites left out from the same rows of the reference table, as in
#' [abc_ppc()]. A model that fits the sites it saw but not the others gives observed test
#' statistics outside the predictive intervals. Use spatial blocks from [spatial_folds()]
#' so that the test sites are not neighbours of the training sites.
#'
#' @param table An `mc_abc_table` from [abc_table()] with `summarise = fold_summarise(...)`.
#' @param observed Named numeric vector with the same statistics for the observed data, from
#'   the same [fold_summarise()] function, or a data frame with one data set per row. The
#'   forests are fitted once per fold for all rows.
#' @param parameters Names of the parameters. By default those of `table`.
#' @param statistics Names of the statistics without the prefix, such as `"alpha"`. By default
#'   all.
#' @param level Probability of the predictive interval.
#' @param n Number of draws from the posterior of each fold.
#' @param ntree Number of trees per forest.
#' @param ncores Number of threads for the forests.
#' @param seed Integer that sets the seed.
#' @return A data frame with one row per observed data set, fold and statistic: `observed`, the median and the
#'   interval of the predicted statistic, `inside` (whether the observed value is in the
#'   interval), `p` (two-sided posterior predictive p-value) and `z`, the observed value minus
#'   the predicted median over the predictive standard deviation. `summary()` averages over
#'   folds and data sets (`n` rows per statistic): the coverage, the mean of `z`, the root mean square of `z` and the fraction of
#'   p-values below 0.05.
#' @family inference
#' @examplesIf requireNamespace("abcrf", quietly = TRUE)
#' set.seed(1)
#' table <- data.frame(a = runif(400))
#' for (f in 1:2) {
#'   table[[paste0("f", f, "_train_s")]] <- table$a + rnorm(400, 0, 0.05)
#'   table[[paste0("f", f, "_test_s")]] <- table$a + rnorm(400, 0, 0.1)
#' }
#' observed <- c(f1_train_s = 0.5, f1_test_s = 0.55, f2_train_s = 0.5, f2_test_s = 0.9)
#' cv <- abc_cv(table, observed, parameters = "a", ntree = 100)
#' cv
#' summary(cv)
#' @export
abc_cv <- function(table, observed, parameters = NULL, statistics = NULL, level = 0.9,
                   n = 500, ntree = 500, ncores = 1, seed = 1) {
  stopifnot(is.data.frame(table), length(level) == 1, level > 0, level < 1, n >= 2)
  parameters <- parameters %||% attr(table, "parameters")
  if (is.null(parameters)) stop("Give `parameters` for a plain data frame.", call. = FALSE)
  pattern <- "^f([0-9]+)_(train|test)_(.+)$"
  columns <- grep(pattern, names(table), value = TRUE)
  if (!length(columns)) {
    stop("`table` has no statistics by fold; see fold_summarise().", call. = FALSE)
  }
  info <- data.frame(column = columns, fold = as.integer(sub(pattern, "\\1", columns)),
                     part = sub(pattern, "\\2", columns), statistic = sub(pattern, "\\3", columns))
  if (!is.null(statistics)) info <- info[info$statistic %in% statistics, , drop = FALSE]
  if (!nrow(info)) stop("None of `statistics` is in `table`.", call. = FALSE)
  if (!is.data.frame(observed)) observed <- as.data.frame(as.list(observed))
  missing <- setdiff(info$column, names(observed))
  if (length(missing)) {
    stop("Statistics missing from `observed`: ", paste(utils::head(missing, 5), collapse = ", "),
         if (length(missing) > 5) ", ...", ".", call. = FALSE)
  }
  probs <- c((1 - level) / 2, (1 + level) / 2)
  old <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, envir = globalenv()))
  out <- lapply(sort(unique(info$fold)), function(f) {
    train <- info[info$fold == f & info$part == "train", ]
    test <- info[info$fold == f & info$part == "test", ]
    finite <- colSums(!is.finite(as.matrix(observed[, train$column, drop = FALSE]))) == 0
    train <- train[finite, , drop = FALSE]
    fold_table <- as.data.frame(table)[, c(parameters, train$column), drop = FALSE]
    fit <- abc_fit(fold_table, parameters, train$column, ntree = ntree, ncores = ncores)
    weights <- forest_weights(fit, observed[, fit$statistics, drop = FALSE])
    rows <- match(rownames(fit$training), rownames(table))
    do.call(rbind, lapply(seq_len(nrow(observed)), function(o) {
      cell_seed(seed + 1000 * f + o)
      drawn <- rows[sample.int(length(rows), n, replace = TRUE, prob = weights[, o])]
      do.call(rbind, lapply(seq_len(nrow(test)), function(i) {
        v <- table[[test$column[i]]][drawn]
        v <- v[is.finite(v)]
        obs <- observed[[test$column[i]]][o]
        if (!length(v) || !is.finite(obs)) return(NULL)
        interval <- unname(stats::quantile(v, probs))
        spread <- stats::sd(v)
        data.frame(observation = o, fold = f, statistic = test$statistic[i], observed = obs,
                   median = stats::median(v), low = interval[1], high = interval[2],
                   inside = obs >= interval[1] & obs <= interval[2],
                   p = min(1, 2 * min(mean(v >= obs), mean(v <= obs))),
                   z = if (isTRUE(spread > 0)) (obs - stats::median(v)) / spread else NA_real_)
      }))
    }))
  })
  out <- do.call(rbind, out)
  out <- out[order(out$observation, out$fold), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("mc_abc_cv", "data.frame")
  out
}

#' @export
summary.mc_abc_cv <- function(object, ...) {
  d <- as.data.frame(unclass(object))
  out <- do.call(rbind, lapply(split(d, d$statistic), function(s) {
    data.frame(statistic = s$statistic[1], n = nrow(s), coverage = mean(s$inside),
               z_mean = mean(s$z, na.rm = TRUE), z_rms = sqrt(mean(s$z^2, na.rm = TRUE)),
               rejected = mean(s$p < 0.05))
  }))
  rownames(out) <- NULL
  out
}

# mean over parameters of the weights that the forests give to the rows of the reference table
forest_weights <- function(fit, observed) {
  Reduce(`+`, lapply(fit$parameters, function(p) {
    training <- data.frame(.parameter = fit$training[[p]], fit$training[fit$statistics])
    w <- as.matrix(stats::predict(fit$forests[[p]], observed, training,
                                  rf.weights = TRUE)$weights)
    sweep(w, 2, colSums(w), "/")
  })) / length(fit$parameters)
}
