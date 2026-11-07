# Checks shared by the analysis functions, with messages that say what to do.

# numeric community without missing or negative values
check_comm <- function(comm) {
  check_numeric_columns(comm, "`comm`", "Move names of sites to row names.")
  comm <- as.matrix(comm) * 1
  if (anyNA(comm)) stop("`comm` has missing values", call. = FALSE)
  if (any(comm < 0)) stop("`comm` has negative values", call. = FALSE)
  comm
}

# numeric predictors in the order of the sites of `comm`, each with variance and, when
# `collinear` is FALSE, none a combination of the others
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
  labels <- colnames(env) %||% seq_len(ncol(env))
  constant <- apply(env, 2, stats::var) == 0
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

check_permutations <- function(permutations, comm) {
  if (is.null(permutations)) return(invisible())
  n <- nrow(comm)
  if (!is.matrix(permutations) || ncol(permutations) != n) {
    stop("`site_permutations` needs one column per site (", n, ") and one row per permutation",
         call. = FALSE)
  }
  if (!all(apply(permutations, 1, function(row) all(sort(row) == seq_len(n))))) {
    stop("each row of `site_permutations` must hold each site once", call. = FALSE)
  }
}

# a loop with a progress bar and an estimate of the time left, which appears only in
# interactive sessions and when the loop takes more than a few seconds
progress_map <- function(x, f, name) {
  purrr::map(x, f, .progress = if (interactive()) name else FALSE)
}

# what the null models used, kept in the results and shown when they are printed
permutation_record <- function(n_perm, site_permutations = NULL) {
  restricted <- !is.null(site_permutations)
  list(species = n_perm, sites = if (restricted) nrow(site_permutations) else n_perm,
       null = if (restricted) "restricted" else "free")
}

permutation_line <- function(record) {
  cat("permutations:", record$sites, record$null, "of sites,", record$species,
      "of species; smallest possible p:",
      format(1 / (min(record$sites, record$species) + 1), digits = 2), "\n")
}
