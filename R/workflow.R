#' Run a grid of simulation scenarios
#'
#' Runs every row of `scenarios` in every replicate and returns one table. Each replicate draws
#' what `setup` returns once, so scenarios are compared on the same tree and landscape. Each
#' cell of the design has its own seed, which comes from `seed`, its replicate and its row, so
#' results do not depend on the order of the cells nor on the number of workers. Cells run in
#' parallel under a plan of the future package, if you choose one. The random number generator
#' of the session is left as it was.
#'
#' @param scenarios Data frame with one scenario per row and one parameter per column.
#' @param simulate Function of `scenario`, a list with one row of `scenarios`, and `shared`, what
#'   `setup` returned, which gives a simulation.
#' @param setup Optional function without arguments that draws what the scenarios of one
#'   replicate share, such as a list with `tree` and `landscape`.
#' @param summarise Function of `sim` and `shared` that gives a named numeric vector. The default
#'   is [summary_stats()], with `shared$tree` if there is one.
#' @param n_rep Number of replicates.
#' @param seed Integer that sets all the seeds.
#' @return Data frame with the columns of `scenarios`, `replicate` and the statistics.
#' @family simulation studies
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' scenarios <- expand.grid(dispersal = c(0.001, 0.1), lambda = c(0, 1))
#' shared_pool <- function() {
#'   list(tree = treeshapesim::sim_tree_shape(10), landscape = landscape_implicit(8))
#' }
#' one_run <- function(scenario, shared) {
#'   simulate_metacommunity(shared$landscape, species_traits(shared$tree, scenario$lambda),
#'                          competition_matrix(shared$tree, "none"),
#'                          dispersal = scenario$dispersal, n_steps = 30, seed_steps = 1:3)
#' }
#' out <- run_scenarios(scenarios, one_run, shared_pool, n_rep = 2)
#' aggregate(alpha ~ dispersal + lambda, out, mean)
#' @export
run_scenarios <- function(scenarios, simulate, setup = NULL,
                          summarise = function(sim, shared) summary_stats(sim, shared$tree),
                          n_rep = 1, seed = 1) {
  stopifnot(is.data.frame(scenarios), nrow(scenarios) >= 1, is.function(simulate), n_rep >= 1,
            length(seed) == 1)
  old <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, envir = globalenv()))
  n <- nrow(scenarios)
  shared <- lapply(seq_len(n_rep), function(r) {
    cell_seed(seed + r - 1)
    if (is.null(setup)) list() else setup()
  })
  cells <- expand.grid(row = seq_len(n), replicate = seq_len(n_rep))
  rows <- progress_map(seq_len(nrow(cells)), function(i) {
    row <- cells$row[i]
    replicate <- cells$replicate[i]
    cell_seed(seed + n_rep + i)
    summarise(simulate(as.list(scenarios[row, , drop = FALSE]), shared[[replicate]]),
              shared[[replicate]])
  }, "Scenarios")
  statistics <- as.data.frame(do.call(rbind, rows))
  cbind(scenarios[cells$row, , drop = FALSE], replicate = cells$replicate, statistics,
        row.names = NULL)
}

# with its generator named, since workers of the future package use another one
cell_seed <- function(seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
}

#' Check whether a simulation has settled
#'
#' Compares the first and the second half of the recorded steps. With phylogenetic competition
#' and high dispersal, regional exclusion can go on for thousands of steps, and the last step of
#' such a run is not an equilibrium. A small change does not prove equilibrium either, since
#' slow exclusion can fall below `tolerance` over a short record: record steps over a long
#' stretch of the run.
#'
#' @param x An `mc_sim` object with `history`, from `record_steps` of
#'   [simulate_metacommunity()], with at least four recorded steps.
#' @param tolerance Largest relative change in a statistic that counts as settled.
#' @return Data frame with one row for `gamma`, the species in all sites, `alpha`, the mean
#'   species per site, and `abundance`, the mean individuals per site, with their means in the
#'   two halves (`early`, `late`), the relative `change` and whether it is `settled`, invisibly.
#'   It gives a warning if any statistic has not settled.
#' @family simulation studies
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(15)
#' sim <- simulate_metacommunity(landscape_implicit(10), species_traits(tree),
#'                               competition_matrix(tree, "phylo"), n_steps = 200,
#'                               seed_steps = 1:5, record_steps = seq(100, 200, by = 10))
#' suppressWarnings(equilibrium_check(sim))
#' @export
equilibrium_check <- function(x, tolerance = 0.05) {
  stopifnot(inherits(x, "mc_sim"), tolerance >= 0)
  steps <- if (is.null(x$history)) 0 else dim(x$history)[3]
  if (steps < 4) {
    stop("`x` needs at least four recorded steps: use `record_steps`", call. = FALSE)
  }
  present <- x$history > 0
  series <- rbind(gamma = apply(apply(present, c(2, 3), any), 2, sum),
                  alpha = apply(present, 3, sum) / dim(present)[1],
                  abundance = apply(x$history, 3, sum) / dim(present)[1])
  half <- seq_len(steps %/% 2)
  early <- rowMeans(series[, half, drop = FALSE])
  late <- rowMeans(series[, -half, drop = FALSE])
  change <- ifelse(early > 0, (late - early) / early, 0)
  out <- data.frame(statistic = rownames(series), early = unname(early), late = unname(late),
                    change = unname(change), settled = unname(abs(change) <= tolerance))
  if (!all(out$settled)) {
    warning("still changing: ", paste(out$statistic[!out$settled], collapse = ", "),
            ". Run more steps", call. = FALSE)
  }
  invisible(out)
}

#' Simulation results as a long table
#'
#' @param x An `mc_sim` object.
#' @param row.names,optional Not used.
#' @param zeros Keep absences?
#' @param ... Not used.
#' @return Data frame with `site`, `species`, `step` and `abundance`, from `history`, or from
#'   `abundance` with `step` as `NA` if no steps were recorded. `tibble::as_tibble()` turns it
#'   into a tibble.
#' @family simulation studies
#' @examples
#' set.seed(1)
#' sim <- simulate_metacommunity(landscape_implicit(4), species_traits(5),
#'                               competition_matrix(5, "none"), n_steps = 10,
#'                               seed_steps = 1:2, record_steps = c(5, 10))
#' head(as.data.frame(sim))
#' @export
as.data.frame.mc_sim <- function(x, row.names = NULL, optional = FALSE, zeros = FALSE, ...) {
  values <- x$history %||% x$abundance
  labels <- dimnames(values)
  steps <- if (is.null(x$history)) NA_integer_ else as.integer(sub("^step", "", labels[[3]]))
  out <- expand.grid(site = labels[[1]], species = labels[[2]], step = steps,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  out$abundance <- as.vector(values)
  if (!zeros) out <- out[out$abundance > 0, , drop = FALSE]
  rownames(out) <- NULL
  out
}
