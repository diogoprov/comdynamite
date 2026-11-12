#' Imperfect detection of a community
#'
#' Turns the true abundances of a simulation into what a survey would record. Each individual
#' is detected in each of `visits` visits with probability `detection`, independently, so the
#' count of a species at a site in a visit is binomial and a species with more individuals is
#' more likely to be found, as in the models of Royle & Nichols (2003) and Royle (2004). The
#' record keeps the largest count over visits, or only whether the species was detected. The
#' same applies to every step of `series`, with new visits at each step. Use it inside
#' `simulate` in [abc_table()], with `detection` drawn from a prior, so that the reference table
#' carries the same observation error as the data.
#'
#' @param x An `mc_sim` object or an [mc_data()] object with true abundances.
#' @param detection Probability that an individual is detected in a visit: one number, or one
#'   per species, named by species.
#' @param visits Number of visits to each site.
#' @param record `"counts"` for the largest count over visits or `"presence"` for detected at
#'   least once (1) or not (0).
#' @return An [mc_data()] object with the recorded `comm` and `series`, and the attribute
#'   `detection`.
#' @references
#' Royle, J. A., & Nichols, J. D. (2003). Estimating abundance from repeated presence–absence
#' data or point counts. *Ecology*, 84, 777–790. \doi{10.1890/0012-9658(2003)084[0777:EAFRPA]2.0.CO;2}
#'
#' Royle, J. A. (2004). N-mixture models for estimating population size from spatially
#' replicated counts. *Biometrics*, 60, 108–115. \doi{10.1111/j.0006-341X.2004.00142.x}
#' @family inference
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' tree <- treeshapesim::sim_tree_shape(8)
#' sim <- simulate_metacommunity(landscape_grid(side = 4), species_traits(tree, 0.5),
#'                               competition_matrix(tree, "none"), dispersal = 0.05,
#'                               n_steps = 40)
#' seen <- observe_detection(sim, detection = 0.1, visits = 3)
#' c(true = sum(sim$abundance > 0), recorded = sum(seen$comm > 0))
#' @export
observe_detection <- function(x, detection, visits = 3, record = c("counts", "presence")) {
  record <- match.arg(record)
  data <- if (inherits(x, "mc_data")) x else as_mc_data(x)
  species <- colnames(data$comm)
  stopifnot(is.numeric(detection), all(detection >= 0 & detection <= 1), length(visits) == 1,
            visits >= 1)
  if (length(detection) == 1 && is.null(names(detection))) {
    detection <- stats::setNames(rep(detection, length(species)), species)
  } else {
    if (is.null(names(detection)) || !all(species %in% names(detection))) {
      stop("`detection` needs one value per species, named by species.", call. = FALSE)
    }
    detection <- detection[species]
  }
  if (any(data$comm != round(data$comm) | data$comm < 0)) {
    stop("Detection needs abundances as counts of individuals.", call. = FALSE)
  }
  survey <- function(abundance) {
    p <- matrix(detection, nrow(abundance), ncol(abundance), byrow = TRUE)
    best <- matrix(0, nrow(abundance), ncol(abundance), dimnames = dimnames(abundance))
    for (v in seq_len(visits)) {
      counts <- stats::rbinom(length(abundance), as.vector(abundance), as.vector(p))
      best <- pmax(best, matrix(counts, nrow(abundance)))
    }
    if (record == "presence") best <- (best > 0) * 1
    best
  }
  data$comm <- survey(data$comm)
  if (!is.null(data$series)) {
    for (t in seq_len(dim(data$series)[3])) data$series[, , t] <- survey(data$series[, , t])
  }
  attr(data, "detection") <- list(detection = detection, visits = visits, record = record)
  data
}
