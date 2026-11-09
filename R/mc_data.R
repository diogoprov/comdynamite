#' Metacommunity data in one object
#'
#' Keeps the community and what describes its species and sites together, with names checked
#' once. The analyses accept an `mc_data` object as their first argument and take from it what
#' the call leaves out, so `pcps_test(d)` uses `d$tree` and `d$env`. Simulations convert with
#' [as_mc_data()].
#'
#' Species are matched by name between `comm`, `tree` and `traits`, and those missing from any
#' of them are dropped with a warning, as in the analyses. Species without individuals stay,
#' since they belong to the regional pool. Sites are matched by name when `env`, `coords`,
#' `areas` or `series` have names, and by position otherwise.
#'
#' @param comm Matrix or data frame with one row per site and one column per species, named by
#'   species. Sites without names are called `site1`, `site2` and so on.
#' @param tree Optional `phylo` object.
#' @param traits Optional data frame, matrix or named vector with one row per species, named by
#'   species.
#' @param env Optional data frame, matrix or vector with one row per site.
#' @param coords Optional matrix or data frame of coordinates with one row per site.
#' @param areas Optional vector with the area or biogeographic event of each site, for
#'   [biogeography_matrix()].
#' @param series Optional array of sites x species x steps, such as the `history` of a
#'   simulation.
#' @param x,object An `mc_data` object, or for [as_mc_data()] an object to convert.
#' @param i,j Sites and species to keep, by position, name or logical vector.
#' @param drop,row.names,optional Not used.
#' @param ... Passed to [mc_data()] by [as_mc_data()]; not used elsewhere.
#' @return `mc_data()` and `as_mc_data()` give a list of class `mc_data` with `comm`, `tree`,
#'   `traits` (data frame), `env` (data frame), `coords` (matrix), `areas` and `series`, each
#'   `NULL` when absent. `[` gives the same for a subset, `summary()` gives [summary_stats()]
#'   and `as.data.frame()` a long table with `site`, `species`, `step` (`NA` without `series`)
#'   and `abundance`, without absences.
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(12)
#' traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100), rmax_range = c(0, 30))
#' env <- seq(0, 100, length.out = 20)
#' comm <- simulate_static(traits, env, incidence = TRUE)$abundance
#' d <- mc_data(comm, tree, traits["optimum"], env = data.frame(env = env))
#' d
#' d[1:10, ]
#' pcps_test(d[rowSums(d$comm) > 0, ], n_perm = 9)
#' @export
mc_data <- function(comm, tree = NULL, traits = NULL, env = NULL, coords = NULL, areas = NULL,
                    series = NULL) {
  comm <- check_comm(comm)
  if (is.null(colnames(comm))) stop("`comm` needs species names in its columns", call. = FALSE)
  if (is.null(rownames(comm))) rownames(comm) <- paste0("site", seq_len(nrow(comm)))
  if (anyDuplicated(rownames(comm))) stop("`comm` has repeated site names", call. = FALSE)
  species <- colnames(comm)
  if (!is.null(tree)) {
    if (!inherits(tree, "phylo")) stop("`tree` must be a `phylo` object", call. = FALSE)
    species <- intersect(species, match_species(comm[, species, drop = FALSE], tree$tip.label,
                                                "`tree`"))
  }
  if (!is.null(traits)) {
    traits <- as_species_table(traits)
    species <- intersect(species, match_species(comm[, species, drop = FALSE], rownames(traits),
                                                "`traits`"))
  }
  sites <- rownames(comm)
  data <- list(comm = comm[, species, drop = FALSE],
               tree = if (!is.null(tree)) ape::keep.tip(tree, species),
               traits = if (!is.null(traits)) traits[species, , drop = FALSE],
               env = if (!is.null(env)) as.data.frame(site_rows(env, sites, "`env`")),
               coords = if (!is.null(coords)) {
                 as.matrix(check_coords(site_rows(coords, sites, "`coords`")))
               },
               areas = if (!is.null(areas)) site_rows(areas, sites, "`areas`"),
               series = if (!is.null(series)) series_of(series, sites, species))
  structure(data, class = "mc_data")
}

#' @rdname mc_data
#' @export
as_mc_data <- function(x, ...) UseMethod("as_mc_data")

#' @rdname mc_data
#' @export
as_mc_data.default <- function(x, ...) mc_data(x, ...)

#' @rdname mc_data
#' @export
as_mc_data.mc_data <- function(x, ...) x

#' @rdname mc_data
#' @export
as_mc_data.mc_sim <- function(x, ...) {
  env <- x$env %||% x$landscape$env
  if (is.list(env) && !is.data.frame(env)) env <- as.data.frame(env)
  if (is.numeric(env) && is.null(dim(env))) env <- data.frame(env = env)
  traits <- if (is.data.frame(x$traits)) x$traits
  if (sim_niche(x)$n_axes > 1) {  # one column with the optima of each axis
    traits <- as.data.frame(lapply(x$traits, function(t) t$optimum),
                            row.names = rownames(x$traits[[1]]))
    names(traits) <- paste0("optimum_", names(env))
  }
  mc_data(x$abundance, tree = x$tree, traits = traits, env = env, coords = x$landscape$coords,
          series = x$history, ...)
}

#' @rdname mc_data
#' @export
print.mc_data <- function(x, ...) {
  cat("<mc_data>", nrow(x$comm), "sites x", ncol(x$comm), "species\n")
  parts <- c(tree = if (!is.null(x$tree)) paste(ape::Ntip(x$tree), "tips"),
             traits = if (!is.null(x$traits)) paste(names(x$traits), collapse = ", "),
             env = if (!is.null(x$env)) paste(names(x$env), collapse = ", "),
             coords = if (!is.null(x$coords)) paste(ncol(x$coords), "columns"),
             areas = if (!is.null(x$areas)) paste(length(unique(x$areas)), "values"),
             series = if (!is.null(x$series)) paste(dim(x$series)[3], "steps"))
  for (part in names(parts)) cat(part, ": ", parts[[part]], "\n", sep = "")
  invisible(x)
}

#' @rdname mc_data
#' @export
`[.mc_data` <- function(x, i, j, drop = FALSE) {
  sites <- if (missing(i)) rownames(x$comm) else rownames(x$comm[i, , drop = FALSE])
  species <- if (missing(j)) colnames(x$comm) else colnames(x$comm[, j, drop = FALSE])
  if (length(species) < 2) stop("keep at least two species", call. = FALSE)
  at <- match(sites, rownames(x$comm))
  x$comm <- x$comm[at, species, drop = FALSE]
  if (!is.null(x$tree)) x$tree <- ape::keep.tip(x$tree, species)
  if (!is.null(x$traits)) x$traits <- x$traits[species, , drop = FALSE]
  if (!is.null(x$env)) x$env <- x$env[at, , drop = FALSE]
  if (!is.null(x$coords)) x$coords <- x$coords[at, , drop = FALSE]
  if (!is.null(x$areas)) x$areas <- x$areas[at]
  if (!is.null(x$series)) x$series <- x$series[at, species, , drop = FALSE]
  x
}

#' @rdname mc_data
#' @export
summary.mc_data <- function(object, ...) summary_stats(object, ...)

#' @rdname mc_data
#' @export
as.data.frame.mc_data <- function(x, row.names = NULL, optional = FALSE, ...) {
  values <- x$series %||% x$comm
  steps <- if (is.null(x$series)) NA_integer_ else seq_len(dim(values)[3])
  if (!is.null(x$series) && !is.null(dimnames(values)[[3]])) {
    steps <- suppressWarnings(as.integer(sub("^step", "", dimnames(values)[[3]])))
  }
  out <- expand.grid(site = rownames(x$comm), species = colnames(x$comm), step = steps,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  out$abundance <- as.vector(values)
  out <- out[out$abundance > 0, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# a part of `data` that a call left out; stops when the object does not have it
data_part <- function(data, part, argument = part) {
  value <- data[[part]]
  if (is.null(value)) {
    stop("the `mc_data` object has no `", part, "`: give `", argument, "`", call. = FALSE)
  }
  value
}

# one predictor or trait from a data frame of an mc_data object
# environment and traits of one niche axis of a simulation, and the number of axes
sim_niche <- function(x, axis = 1) {
  several <- is.list(x$traits) && !is.data.frame(x$traits) && length(x$traits) > 0 &&
    all(vapply(x$traits, is.data.frame, logical(1)))
  list(env = if (is.data.frame(x$env)) x$env[[axis]] else x$landscape$env %||% x$env,
       traits = if (several) x$traits[[axis]] else x$traits,
       n_axes = if (several) length(x$traits) else 1L)
}

one_column <- function(table, part, prefer = NULL) {
  if (is.null(table)) return(NULL)
  numeric <- names(table)[vapply(table, is.numeric, logical(1))]
  column <- if (!is.null(prefer) && prefer %in% numeric) prefer else if (length(numeric) == 1) {
    numeric
  } else {
    stop("the `mc_data` object has ", length(numeric), " numeric columns in `", part,
         "`: give the one to use", call. = FALSE)
  }
  stats::setNames(table[[column]], rownames(table))
}

as_species_table <- function(traits) {
  if (is.null(dim(traits))) {
    if (is.null(names(traits))) stop("`traits` needs species names", call. = FALSE)
    traits <- data.frame(trait = unname(traits), row.names = names(traits))
  }
  traits <- as.data.frame(traits)
  if (.row_names_info(traits) < 0) stop("`traits` needs species names in its rows", call. = FALSE)
  traits
}

# rows of `x` in the order of `sites`: by name when `x` has names, by position otherwise
site_rows <- function(x, sites, name) {
  is_table <- !is.null(dim(x))
  n <- if (is_table) nrow(x) else length(x)
  labels <- if (is_table) {
    if (is.data.frame(x) && .row_names_info(x) < 0) NULL else rownames(x)
  } else {
    names(x)
  }
  if (is.null(labels)) {
    if (n != length(sites)) {
      stop(name, " has ", n, " rows and `comm` ", length(sites), " sites", call. = FALSE)
    }
    at <- seq_len(n)
  } else {
    at <- match(sites, labels)
    if (anyNA(at)) {
      stop(name, " has no row for the sites ", species_list(sites[is.na(at)]), call. = FALSE)
    }
  }
  out <- if (is_table) x[at, , drop = FALSE] else x[at]
  if (is_table) rownames(out) <- sites else names(out) <- sites
  out
}

check_coords <- function(coords) {
  check_numeric_columns(coords, "`coords`", "Coordinates must be numbers.")
  if (!is.numeric(as.matrix(coords))) stop("`coords` must be numbers", call. = FALSE)
  coords
}

series_of <- function(series, sites, species) {
  if (length(dim(series)) != 3) stop("`series` must be an array of sites x species x steps",
                                     call. = FALSE)
  labels <- dimnames(series)
  if (is.null(labels[[2]]) || !all(species %in% labels[[2]])) {
    stop("`series` needs the species of `comm` as names of its second dimension", call. = FALSE)
  }
  rows <- if (is.null(labels[[1]])) seq_along(sites) else match(sites, labels[[1]])
  if (dim(series)[1] != length(sites) || anyNA(rows)) {
    stop("`series` needs the sites of `comm` in its first dimension", call. = FALSE)
  }
  out <- series[rows, species, , drop = FALSE]
  dimnames(out)[[1]] <- sites
  out
}
