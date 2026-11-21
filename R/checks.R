# Checks shared by the analysis functions, with messages that say what to do.

# numeric community without missing or negative values
check_comm <- function(comm) {
  check_numeric_columns(comm, "`comm`", "Move names of sites to row names.")
  comm <- as.matrix(comm) * 1
  if (anyNA(comm)) stop("`comm` has missing values", call. = FALSE)
  if (any(is.infinite(comm))) stop("`comm` has infinite values", call. = FALSE)
  if (any(comm < 0)) stop("`comm` has negative values", call. = FALSE)
  comm
}

# numeric predictors in the order of the sites of `comm`, each with variance and, when
# `collinear` is FALSE, none a combination of the others
# stops unless `tree` is a phylo object (with branch lengths, if asked); `name` is the argument
check_tree <- function(tree, branch_lengths = TRUE, name = "`tree`") {
  if (!inherits(tree, "phylo")) {
    stop(name, " must be a `phylo` object, such as one from treeshapesim::sim_tree_shape()",
         call. = FALSE)
  }
  if (branch_lengths && is.null(tree$edge.length)) {
    stop(name, " needs branch lengths", call. = FALSE)
  }
  invisible(tree)
}

# stops unless `x` is one finite whole number between `min` and `max`; NULL passes when
# `null_ok` is TRUE
check_count <- function(x, name, min = 1, max = Inf, null_ok = FALSE) {
  if (null_ok && is.null(x)) return(invisible(x))
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x != round(x) || x < min ||
        x > max) {
    range <- if (is.finite(max)) paste0(" between ", min, " and ", max) else
      paste0(", ", min, " or more")
    stop("`", name, "` must be one whole number", range, call. = FALSE)
  }
  invisible(x)
}

# stops unless `x` is one number (not NA) between `lower` and `upper`, which are excluded when
# `open` is TRUE (or a pair of logicals for each end); infinite values pass only when
# `finite` is FALSE, and NULL only when `null_ok` is TRUE
check_number <- function(x, name, lower = -Inf, upper = Inf, open = FALSE, finite = TRUE,
                         null_ok = FALSE) {
  if (null_ok && is.null(x)) return(invisible(x))
  open <- rep_len(open, 2)
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && (!finite || is.finite(x)) &&
    (if (open[1]) x > lower else x >= lower) && (if (open[2]) x < upper else x <= upper)
  if (!ok) stop("`", name, "` must be ", number_range(lower, upper, open, finite),
                call. = FALSE)
  invisible(x)
}

# the range of `check_number()` in words
number_range <- function(lower, upper, open, finite) {
  kind <- if (finite) "one finite number" else "one number"
  low <- if (open[1]) "greater than " else ""
  high <- if (open[2]) "less than " else ""
  if (is.finite(lower) && is.finite(upper)) {
    if (!any(open)) return(paste0(kind, " between ", lower, " and ", upper))
    return(paste0(kind, " ", if (open[1]) "greater than " else "from ", lower, " and ",
                  if (open[2]) "less than " else "up to ", upper))
  }
  if (is.finite(lower)) return(paste0(kind, ", ", if (open[1]) low else "", lower,
                                      if (open[1]) "" else " or more"))
  if (is.finite(upper)) return(paste0(kind, ", ", if (open[2]) high else "", upper,
                                      if (open[2]) "" else " or less"))
  kind
}

# stops unless `x` is one TRUE or FALSE
check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE", call. = FALSE)
  }
  invisible(x)
}

# stops unless `x` is one character string (not NA); NULL passes when `null_ok` is TRUE
check_string <- function(x, name, null_ok = FALSE) {
  if (null_ok && is.null(x)) return(invisible(x))
  if (!is.character(x) || length(x) != 1 || is.na(x)) {
    stop("`", name, "` must be one character string", call. = FALSE)
  }
  invisible(x)
}

# TRUE for columns (or a vector) whose spread is zero up to rounding error, which `var(x) == 0`
# misses when the values carry floating-point noise
no_variance <- function(x) {
  x <- as.matrix(x)
  apply(x, 2, function(column) {
    scale <- max(1, abs(column), na.rm = TRUE)
    diff(range(column, na.rm = TRUE)) <= sqrt(.Machine$double.eps) * scale
  })
}

check_env <- function(env, comm, name = "`env`", collinear = FALSE) {
  if (is.null(env)) stop(name, " is NULL", call. = FALSE)
  check_numeric_columns(env, name, paste("Move names of sites to row names and code factors as",
                                         "dummy variables, for example with `model.matrix()`."))
  env <- as.matrix(env)
  if (!is.numeric(env) && !is.logical(env)) stop(name, " is not numeric", call. = FALSE)
  env <- check_sites(env * 1, comm, name)
  if (anyNA(env)) {
    stop(name, " has missing values in ", sum(!stats::complete.cases(env)), " sites. Remove them ",
         "from ", name, " and `comm`", call. = FALSE)
  }
  if (any(is.infinite(env))) stop(name, " has infinite values", call. = FALSE)
  labels <- colnames(env) %||% seq_len(ncol(env))
  constant <- no_variance(env)
  if (any(constant)) {
    stop("predictors without variance in ", name, ": ", paste(labels[constant], collapse = ", "),
         call. = FALSE)
  }
  if (!collinear && qr(sweep(env, 2, colMeans(env)))$rank < ncol(env)) {
    stop("predictors are collinear in ", name, call. = FALSE)
  }
  env
}

# one row per site, matched by name when both have names of sites
check_sites <- function(x, comm, name) {
  if (nrow(x) != nrow(comm)) {
    stop(name, " has ", nrow(x), " rows and `comm` has ", nrow(comm), " sites", call. = FALSE)
  }
  sites <- rownames(comm)
  if (is.null(sites) || is.null(rownames(x)) || anyDuplicated(sites) > 0 ||
        anyDuplicated(rownames(x)) > 0) {
    return(x)
  }
  different <- union(setdiff(rownames(x), sites), setdiff(sites, rownames(x)))
  if (length(different) > 0) {
    stop("sites of ", name, " and `comm` have different names: ", species_list(different),
         call. = FALSE)
  }
  x[sites, , drop = FALSE]
}

check_numeric_columns <- function(x, name, advice) {
  if (!is.data.frame(x)) return(invisible())
  bad <- names(x)[!vapply(x, function(column) is.numeric(column) || is.logical(column), logical(1))]
  if (length(bad) > 0) {
    stop("columns of ", name, " that are not numeric: ", paste(bad, collapse = ", "), ". ", advice,
         call. = FALSE)
  }
}

check_permutations <- function(permutations, comm, name = "`site_permutations`") {
  if (is.null(permutations)) return(invisible())
  n <- nrow(comm)
  if (!is.matrix(permutations) || ncol(permutations) != n) {
    stop(name, " needs one column per site (", n, ") and one row per permutation",
         call. = FALSE)
  }
  if (!all(apply(permutations, 1, function(row) all(sort(row) == seq_len(n))))) {
    stop("each row of ", name, " must hold each site once", call. = FALSE)
  }
}

# a loop with a progress bar and an estimate of the time left, which appears only in
# interactive sessions with purrr installed and when the loop takes more than a few seconds
progress_map <- function(x, f, name, parallel = TRUE) {
  if (parallel && future_workers() > 1) {
    return(future.apply::future_lapply(x, f, future.seed = TRUE))
  }
  if (interactive() && requireNamespace("purrr", quietly = TRUE)) {
    return(purrr::map(x, f, .progress = name)) # nocov
  }
  lapply(x, f)
}

# workers of the plan the user chose with future::plan(), 1 without one
future_workers <- function(installed = requireNamespace("future.apply", quietly = TRUE)) {
  if (!installed) return(1)
  future::nbrOfWorkers()
}

# what the null models used, kept in the results and shown when they are printed
permutation_record <- function(n_perm, site_permutations = NULL, area_permutations = NULL,
                               null = "permutation") {
  restricted <- !is.null(site_permutations)
  given <- is.list(null)
  record <- list(species = n_perm, sites = if (given) length(null) else if (restricted)
                   nrow(site_permutations) else n_perm,
                 null = if (given) "surrogates" else if (null != "permutation") null else
                   if (restricted) "restricted" else
                   "free")
  if (!is.null(area_permutations)) record$areas <- vapply(area_permutations, nrow, 1L)
  record
}

permutation_line <- function(record) {
  sites <- if (record$null %in% c("variogram", "parametric", "surrogates")) {
    paste0(record$sites, " ", if (record$null != "surrogates") paste0(record$null, " "),
           "surrogates of the predictors")
  } else {
    paste0(record$sites, " ", record$null, " of sites")
  }
  cat("permutations: ", sites, ", ", record$species, " of species",
      if (!is.null(record$areas)) {
        paste0(", of areas ", paste(names(record$areas), record$areas, collapse = ", "))
      },
      "; smallest possible p: ",
      format(1 / (min(record$sites, record$species, record$areas) + 1), digits = 2), "\n",
      sep = "")
}

# `alpha` also names the competition matrix, so significance levels are called `level`
deprecated_level <- function(level, alpha) {
  if (is.null(alpha)) return(level)
  warning("`alpha` is deprecated, use `level`", call. = FALSE)
  alpha
}

# sites with individuals. With restricted permutations, empty sites stay in the landscape, so
# draws take the predictors of any site and keep those that land on occupied sites. Row r of
# the result gives the rows of the predictors for the occupied sites
site_draws <- function(comm, n_perm, site_permutations) {
  occupied <- rowSums(comm) > 0
  if (is.null(site_permutations)) {
    return(list(occupied = occupied, rows = t(replicate(n_perm, sample.int(nrow(comm))))))
  }
  list(occupied = occupied, rows = site_permutations[, occupied, drop = FALSE])
}
