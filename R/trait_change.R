#' Eco-evolutionary partition of trait change over time
#'
#' Splits the change of the community-weighted mean trait between recorded steps into the
#' change in the relative abundance of species (sorting), the arrival and loss of species
#' (turnover), the change within species and their covariance, and splits the change within
#' species into a genetic and a plastic part when the genetic means are known. For a simulation
#' with evolution, it can also split the change in the growth rate of each population into the
#' effects of evolution, plasticity and the environment, following the
#' genotype-phenotype-environment partition of Ellner et al. (2011).
#'
#' @details
#' **Community level.** In each site, with \eqn{p_j}{p_j} the relative abundance of species
#' \eqn{j}{j} and \eqn{z_j}{z_j} its mean phenotype, the community mean is
#' \eqn{\bar z = \sum_j p_j z_j}{zbar = sum_j p_j z_j}. Between two steps (0 and 1), its change
#' is split without residual into
#' * `sorting`: \eqn{\sum_S \Delta p_j z_{j0}}{sum_S dp_j z_j0}, over the species \eqn{S}{S}
#'   present at both steps;
#' * `turnover`: \eqn{\sum_{gained} p_{j1} z_{j1} - \sum_{lost} p_{j0} z_{j0}}{sum_gained p_j1
#'   z_j1 - sum_lost p_j0 z_j0}, the species that arrived or were lost;
#' * `within`: \eqn{\sum_S p_{j0} \Delta z_j}{sum_S p_j0 dz_j}, the change within species;
#' * `covariance`: \eqn{\sum_S \Delta p_j \Delta z_j}{sum_S dp_j dz_j}, positive when the
#'   species that increased are those whose trait changed in the same direction as the mean.
#'
#' This is the temporal counterpart of the partition among sites of [trait_flex()]
#' (Lepš et al. 2011). Among the species present at both steps it follows the Price (1970)
#' equation, with the relative abundance in place of fitness: `sorting` is its selection term
#' (the covariance between relative growth and trait) and `within` plus `covariance` its
#' transmission term (the change within species, weighted by the abundances at step 1), as in
#' the continuous-time form used by Ellner et al. (2011). As the phenotype is the genetic
#' mean plus the plastic part, \eqn{z_j = g_j + b_j (E - E_0)}{z_j = g_j + b_j (E - E0)},
#' `within` and `covariance` are
#' each the sum of a genetic part (`genetic`, `cov_genetic`, from \eqn{\Delta g_j}{dg_j}) and a
#' plastic part (`plastic`, `cov_plastic`, from the change in the plastic part). Without
#' genetic means (a data frame without `genetic`), these parts are `NA`. Sites empty at
#' either step give `NA`.
#'
#' In a simulation, the phenotype of a population is its genetic mean plus its plastic
#' response to the environment of the site at that step, as in its growth; the `phenotype`
#' recorded by [simulate_metacommunity()] after dispersal, which mixes in immigrants still
#' carrying the response to their site of origin, is not used.
#'
#' **Population level (`growth = TRUE`).** For each population present at both steps, the
#' response is its per-capita growth rate, \eqn{X = \log \lambda}{X = log(lambda)}, with
#' \eqn{\lambda}{lambda} the Gaussian niche response of the simulation (the growth rate of a
#' population of mean phenotype \eqn{z}{z} and phenotypic variance \eqn{G + V_E}{G + VE} in
#' environment \eqn{E}{E}). Ellner et al. (2011, eqns 5 and 6) split the rate of change of
#' such a response into evolution (the change of the genetic mean), non-heritable trait change
#' (here plasticity) and the direct effect of the environment, through its partial
#' derivatives in continuous time. For changes between samples, Hairston et al. (2005, eqn 8)
#' compute the effect of evolution and of ecological change on \eqn{X}{X} by moving one
#' variable at a time from its value at step 0 to its value at step 1 and averaging the two
#' orders (evolution first or ecology first). The same is done here with three components
#' (genetic mean, plastic part, environment): the contribution of each is its effect averaged
#' over the six orders of change (the Shapley value), which reduces to their eqn 8 with two
#' components, so that `genetic`, `plastic` and `environment` add up to the change in
#' \eqn{X}{X} exactly. The phenotypic variance is held at
#' \eqn{G + V_E}{G + VE} of the settings, so that \eqn{X}{X} depends on the trait only through
#' its mean (the "population mean assumption" of Ellner et al. 2011); with `mutation_var`,
#' `loci` or `variance_selection`, where the variance of each population changes, the
#' partition ignores that change. It needs a simulation with one evolving trait and a table
#' of traits with `sigma` and `r_max`.
#'
#' @param x A simulation with `evolution` and `record_steps` (from [simulate_metacommunity()]
#'   or [simulate_evolving_metacommunity()]), or a data frame with one row per population and
#'   step and the columns named by `site`, `species`, `step`, `abundance`, `trait` and, optionally,
#'   `genetic`.
#' @param steps The steps to compare, in increasing order; consecutive pairs are compared. At
#'   least two. `NULL` (the default) takes all recorded steps.
#' @param growth Also split the change in the growth rate of each population? `TRUE` or
#'   `FALSE`; `NULL` (the default) is `TRUE` for a simulation with one evolving trait and
#'   `FALSE` otherwise. Not available for data frames.
#' @param trait For a data frame, the name of the column with the mean phenotype of each
#'   population. For a simulation with several evolving traits, the name or number of the
#'   trait; ignored with one trait.
#' @param genetic For a data frame, the name of a column with the genetic mean of each
#'   population (from a common garden, for example), or `NULL` to leave the change within
#'   species unsplit. One character string or `NULL`.
#' @param site,species,step,abundance For a data frame, the names of its columns with the
#'   site, the species, the step (numbers) and the abundance of each population. One
#'   character string each.
#' @return A list of class `trait_change` with `community`, a data frame with one row per
#'   site and pair of steps (`site`, `from`, `to`, `mean_from`, `mean_to`, `total`, `sorting`,
#'   `turnover`, `within`, `covariance`, `genetic`, `plastic`, `cov_genetic`, `cov_plastic`),
#'   and `growth`, a data frame with one row per population present at both steps of each
#'   pair (`site`, `species`, `from`, `to`, `total`, `genetic`, `plastic`, `environment`) or
#'   `NULL`. [summary()][summary.trait_change] averages them over sites and populations, and
#'   [plot()][plot.trait_change] draws the averages.
#' @references Ellner, S. P., Geber, M. A., & Hairston, N. G. (2011). Does rapid evolution
#'   matter? Measuring the rate of contemporary evolution and its impacts on ecological
#'   dynamics. *Ecology Letters*, 14, 603--614. \doi{10.1111/j.1461-0248.2011.01616.x}
#'
#'   Hairston, N. G., Ellner, S. P., Geber, M. A., Yoshida, T., & Fox, J. A. (2005). Rapid
#'   evolution and the convergence of ecological and evolutionary time. *Ecology Letters*, 8,
#'   1114--1127. \doi{10.1111/j.1461-0248.2005.00812.x}
#'
#'   Lepš, J., de Bello, F., Šmilauer, P., & Doležal, J. (2011). Community trait response to
#'   environment: disentangling species turnover vs intraspecific trait variability effects.
#'   *Ecography*, 34, 856--863. \doi{10.1111/j.1600-0587.2010.06904.x}
#'
#'   Price, G. R. (1970). Selection and covariance. *Nature*, 227, 520--521.
#'   \doi{10.1038/227520a0}
#' @seealso [trait_flex()] for the partition among sites, [evolution_control()].
#' @family trait variation among individuals
#' @examples
#' set.seed(8)
#' ls <- landscape_grid(4)
#' tr <- ape::rcoal(4)
#' warming <- outer(ls$env, 0.5 * seq_len(60), "+")
#' sim <- simulate_metacommunity(ls, species_traits(tr), competition_matrix(tr, "phylo"),
#'                               n_steps = 60, record_steps = c(20, 40, 60),
#'                               env_time = warming,
#'                               evolution = evolution_control(1, plasticity = 0.3))
#' res <- trait_change(sim)
#' summary(res)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(res)
#'
#' # field data: phenotypes over time, without genetic means
#' d <- data.frame(site = "a", species = rep(c("x", "y"), 2), step = rep(1:2, each = 2),
#'                 abundance = c(10, 10, 5, 15), trait = c(1, 3, 1.5, 3))
#' trait_change(d)$community
#' @export
trait_change <- function(x, steps = NULL, growth = NULL, trait = "trait", genetic = NULL,
                         site = "site", species = "species", step = "step",
                         abundance = "abundance") {
  cite_method("trait_change")
  if (!is.null(growth)) check_flag(growth, "growth")
  if (is.data.frame(x)) {
    if (isTRUE(growth)) {
      stop("`growth = TRUE` needs a simulation with evolution, not a data frame",
           call. = FALSE)
    }
    state <- change_data_state(x, trait, genetic, site, species, step, abundance)
    growth <- FALSE
  } else {
    state <- change_sim_state(x, trait)
    if (is.null(growth)) growth <- state$univariate && change_has_niche(x)
    if (growth && !state$univariate) {
      stop("`growth = TRUE` needs a simulation with one evolving trait", call. = FALSE)
    }
    if (growth && !change_has_niche(x)) {
      stop("`growth = TRUE` needs the traits of the simulation (`sigma` and `r_max`)",
           call. = FALSE)
    }
  }
  steps <- change_steps(steps, state$steps)
  pairs <- seq_len(length(steps) - 1)
  community <- do.call(rbind, lapply(pairs, function(i) {
    change_community(state$at(steps[i]), state$at(steps[i + 1]), steps[i], steps[i + 1])
  }))
  growth_tab <- if (growth) {
    niche <- change_niche(x, state)
    do.call(rbind, lapply(pairs, function(i) {
      change_growth(state$at(steps[i]), state$at(steps[i + 1]), steps[i], steps[i + 1], niche)
    }))
  }
  rownames(community) <- NULL
  if (!is.null(growth_tab)) rownames(growth_tab) <- NULL
  structure(list(community = community, growth = growth_tab, steps = steps,
                 split = state$split),
            class = "trait_change")
}

#' the steps to compare: all recorded ones, or those asked for, which must be recorded
#' @noRd
change_steps <- function(steps, recorded) {
  if (is.null(steps)) steps <- recorded
  if (!is.numeric(steps) || anyNA(steps) || length(steps) < 2) {
    stop("`steps` must be at least two numbers", call. = FALSE)
  }
  if (is.unsorted(steps, strictly = TRUE)) {
    stop("`steps` must be in increasing order, without repeats", call. = FALSE)
  }
  missing <- setdiff(steps, recorded)
  if (length(missing)) {
    stop("step(s) ", paste(missing, collapse = ", "), " not recorded; recorded: ",
         paste(recorded, collapse = ", "), call. = FALSE)
  }
  steps
}

#' the state of a simulation at its recorded steps: a function of the step giving the
#' abundances, genetic means, plastic parts (sites x species) and environment of the sites
#' @noRd
change_sim_state <- function(x, trait) {
  ev <- x$evolution
  if (!is.list(x) || is.null(ev)) {
    stop("`x` must be a simulation with `evolution` or a data frame", call. = FALSE)
  }
  if (is.null(x$history) || is.null(ev$history) || is.null(ev$environment_history)) {
    stop("`x` has no recorded steps: run the simulation with `record_steps`", call. = FALSE)
  }
  univariate <- length(dim(ev$genetic)) == 2
  axis <- 1
  if (!univariate) {
    axes <- dimnames(ev$genetic)[[3]]
    axis <- if (is.numeric(trait)) trait else match(trait, axes)
    if (length(axis) != 1 || is.na(axis) || axis < 1 || axis > length(axes)) {
      stop("`trait` must name one of the traits of the simulation (",
           paste(axes, collapse = ", "), ") or give its number", call. = FALSE)
    }
  }
  labels <- dimnames(ev$history)[[length(dim(ev$history))]]
  labels <- intersect(labels, dimnames(x$history)[[3]])
  steps <- as.numeric(sub("^step", "", labels))
  at <- function(s) {
    a <- sampled_abundance(x, s)
    label <- step_label(s)
    if (univariate) {
      st <- evolved_at(ev, s)
      g <- st$genetic
      q <- st$plastic
      e <- unname(ev$environment_history[, label])
    } else {
      g <- ev$history[, , axis, label]
      g <- matrix(g, dim(ev$history)[1], dimnames = dimnames(ev$genetic)[1:2])
      e <- unname(ev$environment_history[, axis, label])
      b <- ev$plasticity[colnames(g), axis]
      q <- outer(e - ev$reference[axis], b)
      q[is.na(g)] <- NA_real_
      dimnames(q) <- dimnames(g)
    }
    sp <- colnames(a)
    list(abundance = a, genetic = g[rownames(a), sp, drop = FALSE],
         plastic = q[rownames(a), sp, drop = FALSE], environment = e)
  }
  list(at = at, steps = steps, univariate = univariate, split = TRUE)
}

#' the state of a long data frame of populations, as change_sim_state()
#' @noRd
change_data_state <- function(x, trait, genetic, site, species, step, abundance) {
  check_string(trait, "trait")
  check_string(genetic, "genetic", null_ok = TRUE)
  check_string(site, "site")
  check_string(species, "species")
  check_string(step, "step")
  check_string(abundance, "abundance")
  check_data_columns(x, c(site, species, step, abundance, trait, genetic), name = "x")
  if (!is.numeric(x[[step]]) || anyNA(x[[step]])) {
    stop("column `", step, "` must be numbers without missing values", call. = FALSE)
  }
  a <- x[[abundance]]
  if (!is.numeric(a) || anyNA(a) || any(a < 0) || any(!is.finite(a))) {
    stop("column `", abundance, "` must be finite non-negative numbers", call. = FALSE)
  }
  present <- a > 0
  for (col in c(trait, genetic)) {
    v <- x[[col]]
    if (!is.numeric(v)) stop("column `", col, "` must be numeric", call. = FALSE)
    if (any(!is.finite(v[present]))) {
      stop("column `", col, "` must be finite where the abundance is above 0", call. = FALSE)
    }
  }
  keys <- paste(x[[site]], x[[species]], x[[step]], sep = "\r")
  if (anyDuplicated(keys)) {
    stop("`x` must have one row per site, species and step", call. = FALSE)
  }
  sites <- unique(as.character(x[[site]]))
  spp <- unique(as.character(x[[species]]))
  steps <- sort(unique(x[[step]]))
  fill <- function(values, s, empty) {
    m <- matrix(empty, length(sites), length(spp), dimnames = list(sites, spp))
    rows <- x[[step]] == s
    idx <- cbind(match(as.character(x[[site]][rows]), sites),
                 match(as.character(x[[species]][rows]), spp))
    m[idx] <- values[rows]
    m
  }
  at <- function(s) {
    ab <- fill(a, s, 0)
    z <- fill(x[[trait]], s, NA_real_)
    z[ab <= 0] <- NA_real_
    if (is.null(genetic)) {
      list(abundance = ab, genetic = z, plastic = NULL)
    } else {
      g <- fill(x[[genetic]], s, NA_real_)
      g[ab <= 0] <- NA_real_
      list(abundance = ab, genetic = g, plastic = z - g)
    }
  }
  list(at = at, steps = steps, univariate = TRUE, split = !is.null(genetic))
}

#' the community-level partition between two states, one row per site
#' @noRd
change_community <- function(s0, s1, from, to) {
  a0 <- s0$abundance
  a1 <- s1$abundance
  sites <- rownames(a0)
  spp <- union(colnames(a0), colnames(a1))
  widen <- function(m, empty) {
    out <- matrix(empty, length(sites), length(spp), dimnames = list(sites, spp))
    if (!is.null(m)) out[, colnames(m)] <- m[sites, , drop = FALSE]
    out
  }
  split <- !is.null(s0$plastic)
  a0 <- widen(a0, 0)
  a1 <- widen(a1, 0)
  g0 <- widen(s0$genetic, NA_real_)
  g1 <- widen(s1$genetic, NA_real_)
  q0 <- if (split) widen(s0$plastic, NA_real_) else 0 * g0
  q1 <- if (split) widen(s1$plastic, NA_real_) else 0 * g1
  p0 <- a0 / rowSums(a0)
  p1 <- a1 / rowSums(a1)
  z0 <- g0 + q0
  z1 <- g1 + q1
  in0 <- a0 > 0
  in1 <- a1 > 0
  both <- in0 & in1
  # weighted sums over a set of species, 0 where the set is empty
  wsum <- function(w, v, keep) rowSums(ifelse(keep, w * v, 0))
  dp <- p1 - p0
  out <- data.frame(
    site = sites, from = from, to = to,
    mean_from = wsum(p0, z0, in0), mean_to = wsum(p1, z1, in1),
    total = NA_real_,
    sorting = wsum(dp, z0, both),
    turnover = wsum(p1, z1, in1 & !in0) - wsum(p0, z0, in0 & !in1),
    within = wsum(p0, z1 - z0, both),
    covariance = wsum(dp, z1 - z0, both),
    genetic = if (split) wsum(p0, g1 - g0, both) else NA_real_,
    plastic = if (split) wsum(p0, q1 - q0, both) else NA_real_,
    cov_genetic = if (split) wsum(dp, g1 - g0, both) else NA_real_,
    cov_plastic = if (split) wsum(dp, q1 - q0, both) else NA_real_
  )
  out$total <- out$mean_to - out$mean_from
  empty <- rowSums(a0) <= 0 | rowSums(a1) <= 0
  out[empty, -(1:3)] <- NA_real_
  out
}

#' does the simulation carry a table of traits with the niche of each species?
#' @noRd
change_has_niche <- function(x) {
  is.data.frame(x$traits) && all(c("sigma", "r_max") %in% names(x$traits))
}

#' the niche breadth, maximum growth and phenotypic variance of each species of a simulation
#' with one evolving trait, named by species
#' @noRd
change_niche <- function(x, state) {
  ev <- x$evolution
  species <- colnames(state$at(state$steps[1])$abundance)
  missing <- setdiff(species, rownames(x$traits))
  if (length(missing)) {
    stop("the traits of the simulation lack species ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  per <- function(value) {
    if (!is.null(names(value))) return(unname(value[species]))
    evolution_each(value, length(species))
  }
  g <- per(ev$genetic_var %||% ev$settings$genetic_var)
  ve <- per(ev$env_var %||% ev$settings$env_var %||% 0)
  list(sigma = stats::setNames(x$traits[species, "sigma"], species),
       r_max = stats::setNames(x$traits[species, "r_max"], species),
       v = stats::setNames(g + ve, species))
}

#' log growth rate of mean phenotypes p (sites x species) in environment e (one value per site)
#' under the Gaussian niche of the simulation (log of evolution_growth(), without underflow)
#' @noRd
change_log_growth <- function(p, e, sigma, r_max, v) {
  width <- sigma^2 + v
  shrink <- ifelse(is.finite(sigma), 0.5 * log(sigma^2 / width), 0)
  out <- -sweep((e - p)^2, 2, 2 * width, "/")
  sweep(out, 2, log(r_max) + shrink, "+")
}

#' the population-level partition between two states: Shapley values of the genetic mean, the
#' plastic part and the environment for the log growth rate, one row per population present
#' at both steps
#' @noRd
change_growth <- function(s0, s1, from, to, niche) {
  spp <- intersect(colnames(s0$abundance), colnames(s1$abundance))
  both <- s0$abundance[, spp, drop = FALSE] > 0 & s1$abundance[, spp, drop = FALSE] > 0 &
    !is.na(s0$genetic[, spp, drop = FALSE]) & !is.na(s1$genetic[, spp, drop = FALSE])
  g <- list(s0$genetic[, spp, drop = FALSE], s1$genetic[, spp, drop = FALSE])
  q <- list(s0$plastic[, spp, drop = FALSE], s1$plastic[, spp, drop = FALSE])
  e <- list(s0$environment, s1$environment)
  x_at <- function(i, j, k) {
    change_log_growth(g[[i]] + q[[j]], e[[k]], niche$sigma[spp], niche$r_max[spp],
                      niche$v[spp])
  }
  # X at the eight corners: index 1 for step 0 and 2 for step 1 of each component
  corner <- array(list(), c(2, 2, 2))
  for (i in 1:2) for (j in 1:2) for (k in 1:2) corner[[i, j, k]] <- x_at(i, j, k)
  shapley <- function(comp) {
    out <- 0
    for (u in 1:2) for (w in 1:2) {
      # the other two components at (u, w); weight 1/3 when both stay or both moved, 1/6 else
      weight <- if (u == w) 1 / 3 else 1 / 6
      new <- switch(comp, g = corner[[2, u, w]], q = corner[[u, 2, w]], e = corner[[u, w, 2]])
      old <- switch(comp, g = corner[[1, u, w]], q = corner[[u, 1, w]], e = corner[[u, w, 1]])
      out <- out + weight * (new - old)
    }
    out
  }
  parts <- list(genetic = shapley("g"), plastic = shapley("q"), environment = shapley("e"))
  total <- corner[[2, 2, 2]] - corner[[1, 1, 1]]
  idx <- which(both, arr.ind = TRUE)
  if (!nrow(idx)) {
    return(data.frame(site = character(), species = character(), from = numeric(),
                      to = numeric(), total = numeric(), genetic = numeric(),
                      plastic = numeric(), environment = numeric()))
  }
  data.frame(site = rownames(both)[idx[, 1]], species = spp[idx[, 2]], from = from, to = to,
             total = total[idx], genetic = parts$genetic[idx], plastic = parts$plastic[idx],
             environment = parts$environment[idx])
}

#' @export
print.trait_change <- function(x, digits = 3, ...) {
  s <- summary(x)
  cat("Eco-evolutionary partition of trait change over", length(x$steps), "recorded steps (",
      paste(x$steps, collapse = ", "), ")\n\n")
  cat("Community-weighted mean, averaged over sites:\n")
  print(s$community, digits = digits, row.names = FALSE)
  if (!is.null(s$growth)) {
    cat("\nGrowth rate (log lambda) of populations, averaged over populations:\n")
    print(s$growth, digits = digits, row.names = FALSE)
  }
  invisible(x)
}

#' Summarise an eco-evolutionary partition of trait change
#'
#' Averages the partition of [trait_change()] over sites (community level) and over
#' populations (population level) for each pair of steps.
#'
#' @param object Output of [trait_change()].
#' @param ... Not used.
#' @return A list with `community`, a data frame with one row per pair of steps, the number
#'   of sites with both steps occupied (`sites`) and the mean of each component over those
#'   sites, and `growth`, the same over populations (with `populations` and the mean absolute
#'   value of each component, `abs_genetic`, `abs_plastic` and `abs_environment`) or `NULL`.
#' @examples
#' d <- data.frame(site = "a", species = rep(c("x", "y"), 2), step = rep(1:2, each = 2),
#'                 abundance = c(10, 10, 5, 15), trait = c(1, 3, 1.5, 3))
#' summary(trait_change(d))
#' @export
summary.trait_change <- function(object, ...) {
  comp <- c("total", "sorting", "turnover", "within", "covariance", "genetic", "plastic",
            "cov_genetic", "cov_plastic")
  com <- object$community
  ok <- !is.na(com$total)
  by_pair <- split(com[ok, ], factor(paste(com$from, com$to)[ok],
                                     unique(paste(com$from, com$to))))
  community <- do.call(rbind, lapply(by_pair, function(d) {
    data.frame(from = d$from[1], to = d$to[1], sites = nrow(d),
               t(colMeans(d[comp])), check.names = FALSE)
  }))
  rownames(community) <- NULL
  growth <- NULL
  if (!is.null(object$growth) && nrow(object$growth)) {
    gr <- object$growth
    by_pair <- split(gr, factor(paste(gr$from, gr$to), unique(paste(gr$from, gr$to))))
    growth <- do.call(rbind, lapply(by_pair, function(d) {
      data.frame(from = d$from[1], to = d$to[1], populations = nrow(d),
                 total = mean(d$total), genetic = mean(d$genetic), plastic = mean(d$plastic),
                 environment = mean(d$environment), abs_genetic = mean(abs(d$genetic)),
                 abs_plastic = mean(abs(d$plastic)),
                 abs_environment = mean(abs(d$environment)))
    }))
    rownames(growth) <- NULL
  }
  list(community = community, growth = growth)
}

#' Plot an eco-evolutionary partition of trait change
#'
#' Draws, for each pair of steps, the components of the change of the community-weighted
#' mean trait (averaged over sites) as stacked bars, with the total change as a point; or the
#' same for the growth rate of populations. Components that lower the mean are drawn below
#' zero. Needs the package ggplot2.
#'
#' @param x Output of [trait_change()].
#' @param level `"community"` (the default) or `"growth"`.
#' @param colours Colours of the components, in the order in which they are drawn: sorting,
#'   turnover, genetic, plastic and covariance (or within and covariance without genetic
#'   means) for `"community"`; genetic, plastic and environment for `"growth"`. `NULL` uses a
#'   palette safe for colour blindness.
#' @param ... Not used.
#' @return A ggplot object, which can be changed further with ggplot2 functions.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' d <- data.frame(site = "a", species = rep(c("x", "y"), 3), step = rep(1:3, each = 2),
#'                 abundance = c(10, 10, 5, 15, 2, 18), trait = c(1, 3, 1.5, 3, 2, 3.5))
#' plot(trait_change(d))
#' @export
plot.trait_change <- function(x, level = c("community", "growth"), colours = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot() of trait_change() needs the package ggplot2", call. = FALSE) # nocov
  }
  level <- match.arg(level)
  s <- summary(x)
  if (level == "growth") {
    if (is.null(s$growth)) {
      stop("`x` has no partition of the growth rate: use `growth = TRUE`", call. = FALSE)
    }
    tab <- s$growth
    parts <- c("genetic", "plastic", "environment")
    labels <- c("Evolution", "Plasticity", "Environment")
    y_lab <- "Mean change in log growth rate"
  } else {
    tab <- s$community
    if (isTRUE(x$split)) {
      tab$covariance <- tab$cov_genetic + tab$cov_plastic
      parts <- c("sorting", "turnover", "genetic", "plastic", "covariance")
      labels <- c("Sorting", "Turnover", "Genetic", "Plastic", "Covariance")
    } else {
      parts <- c("sorting", "turnover", "within", "covariance")
      labels <- c("Sorting", "Turnover", "Within species", "Covariance")
    }
    y_lab <- "Mean change in community-weighted trait"
  }
  if (is.null(colours)) {
    colours <- c("#0072B2", "#56B4E9", "#D55E00", "#E69F00", "grey60")[seq_along(parts)]
  }
  if (length(colours) < length(parts)) {
    stop("`colours` needs ", length(parts), " colours", call. = FALSE)
  }
  tab$interval <- factor(paste(tab$from, tab$to, sep = "\u2013"),
                         unique(paste(tab$from, tab$to, sep = "\u2013")))
  long <- data.frame(interval = rep(tab$interval, length(parts)),
                     part = factor(rep(labels, each = nrow(tab)), labels),
                     value = unlist(tab[parts], use.names = FALSE))
  .data <- ggplot2::.data
  ggplot2::ggplot(long, ggplot2::aes(x = .data$interval, y = .data$value)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_col(ggplot2::aes(fill = .data$part), width = 0.7) +
    ggplot2::geom_point(data = tab, ggplot2::aes(x = .data$interval, y = .data$total),
                        shape = 21, size = 2.5, fill = "white", colour = "black") +
    ggplot2::scale_fill_manual(values = stats::setNames(colours[seq_along(parts)], labels),
                               name = NULL) +
    ggplot2::labs(x = "Steps", y = y_lab) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), legend.position = "top")
}
