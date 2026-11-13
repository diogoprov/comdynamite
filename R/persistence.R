#' Metapopulation persistence
#'
#' Estimates, for each species, the rates of colonisation \eqn{c}{c} and local extirpation
#' \eqn{e}{e} of a two-state Markov model of presence in sites over time, and from them the two
#' axes of metapopulation persistence of Januario et al. (2025): the equilibrium occupancy
#' \eqn{\epsilon = \mathrm{logit}(c / (c + e))}{epsilon = logit(c / (c + e))} and the
#' resilience \eqn{\tau = \log(c + e)}{tau = log(c + e)}, the log of the rate of turnover. The
#' estimates are the maximum likelihood ones: \eqn{c}{c} is the share of the site-intervals in
#' which an absent species arrived, and \eqn{e}{e} the share of those in which a present one was
#' lost.
#'
#' For [simulate_evolving_metacommunity()], transitions come from its `series`, between
#' consecutive generations, for every lineage (extinct ones included, with the sites lost in
#' the generation they died out); a site that a mother passed to a new species counts neither
#' as an extirpation nor as exposure. For [simulate_metacommunity()], they come from
#' the abundances at `record_steps` (`history`), between consecutive records, which should be
#' evenly spaced (seeding from the pool, `seed_steps`, counts as colonisation). An array of sites x species x times can be given directly.
#'
#' @param x An `mc_evolving` or `mc_sim` object, or an array (sites x species x times) of
#'   abundances or presences.
#' @param window Number of the last intervals (generations, or intervals between records) to
#'   use; `NULL` for all.
#' @param ... Unused.
#' @return A data frame with one row per species (lineage, for `mc_evolving`): `species`, the
#'   number of `intervals` used, `colonisation`, `extirpation`, the equilibrium occupancy
#'   \eqn{c / (c + e)}{c / (c + e)} as `occupancy`, `epsilon` and `tau`. For `mc_evolving`, also
#'   `extant` and the good species each extant lineage belongs to (`good_species`). Rates are
#'   `NA` without any interval in the corresponding state, and `epsilon` is infinite when one of
#'   the rates is 0.
#' @references Januario, M., Pinsky, M. L., & Rabosky, D. L. (2025). The metapopulation bridge
#'   to macroevolutionary speciation rates: a conceptual framework and empirical test. *Ecology
#'   Letters*, 28, e70021. \doi{10.1111/ele.70021}
#' @family simulate metacommunities
#' @examples
#' set.seed(1)
#' ls <- landscape_grid(5, disp_scale = 0.5)
#' ancestor <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
#' sim <- simulate_evolving_metacommunity(ls, ancestor, dispersal = 0.01, extirpation = 0.02,
#'                                        n_steps = 200,
#'                                        speciation = speciation_control("fission", rate = 0.02))
#' head(persistence(sim, window = 50))
#' @export
persistence <- function(x, window = NULL, ...) UseMethod("persistence")

#' @export
persistence.mc_evolving <- function(x, window = NULL, ...) {
  s <- x$series
  if (is.null(s)) {
    stop("`x` has no `series`: rerun it with this version of the package", call. = FALSE)
  }
  if (!is.null(window)) {
    stopifnot(length(window) == 1, window >= 1)
    s <- s[s$generation > x$generations - window, , drop = FALSE]
  }
  s <- s[!is.na(s$present_before), , drop = FALSE]
  lineage <- factor(s$lineage, unique(x$species$species))
  sums <- function(v) tapply(s[[v]], lineage, sum)
  out <- persistence_table(levels(lineage), as.vector(table(lineage)), sums("colonised"),
                           sums("absent_before"), sums("extirpated"), sums("present_before"))
  out$extant <- x$species$fate[match(out$species, x$species$species)] == "extant"
  out$good_species <- unname(x$lineage_species[out$species])
  out
}

#' @export
persistence.mc_sim <- function(x, window = NULL, ...) {
  if (is.null(x$history) || dim(x$history)[3] < 2) {
    stop("`x` needs at least two recorded steps: use `record_steps`", call. = FALSE)
  }
  persistence.default(x$history, window)
}

#' @export
persistence.default <- function(x, window = NULL, ...) {
  if (!is.array(x) || length(dim(x)) != 3 || dim(x)[3] < 2) {
    stop("`x` must be a simulation or an array of sites x species x times (two or more)",
         call. = FALSE)
  }
  times <- dim(x)[3]
  if (!is.null(window)) {
    stopifnot(length(window) == 1, window >= 1)
    x <- x[, , max(1, times - window):times, drop = FALSE]
    times <- dim(x)[3]
  }
  pres <- x > 0
  a <- pres[, , -times, drop = FALSE]
  b <- pres[, , -1, drop = FALSE]
  each <- function(y) apply(y, 2, sum)
  species <- dimnames(x)[[2]] %||% paste0("sp", seq_len(dim(x)[2]))
  persistence_table(species, rep(times - 1, dim(x)[2]), each(!a & b), each(!a),
                    each(a & !b), each(a))
}

persistence_table <- function(species, intervals, colonised, absent, extirpated, present) {
  col <- ifelse(absent > 0, colonised / absent, NA_real_)
  ext <- ifelse(present > 0, extirpated / present, NA_real_)
  occ <- col / (col + ext)
  data.frame(species = species, intervals = as.integer(intervals),
             colonisation = as.vector(col), extirpation = as.vector(ext),
             occupancy = as.vector(occ), epsilon = as.vector(stats::qlogis(occ)),
             tau = as.vector(log(col + ext)), stringsAsFactors = FALSE)
}
