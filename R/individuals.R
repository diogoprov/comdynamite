#' Individuals from a simulated metacommunity
#'
#' Draws individuals from the populations of a simulated metacommunity and gives each one a
#' trait value, so that intraspecific variation and local adaptation can be studied with the
#' T-statistics of Violle et al. (2012) ([tstats()]). The trait is the environmental optimum of
#' [species_traits()], the one that sets growth in [simulate_metacommunity()]. The value of
#' individual `k` of species `j` at site `i` is
#'
#' `z_ijk = theta_j + h_j (E_i - theta_j) + d_ij + e_ijk`,
#'
#' where `theta_j` is the optimum of the species, `E_i` the environment of the site,
#' `h_j` (`local_adaptation`) the fraction of the gap between them that the population closes,
#' `d_ij ~ N(0, population_sd^2)` a deviation of the population that does not track the
#' environment (drift, founder effects), and `e_ijk ~ N(0, intraspecific_sd^2)` the variation
#' among individuals of a population. With `local_adaptation = 0` and `population_sd = 0` all
#' populations of a species share its mean; with `local_adaptation = 1` every population sits on
#' the environment of its site and species no longer differ within a site.
#'
#' The traits of individuals are drawn after the simulation and do not feed back on the
#' dynamics: abundances are those of the species means. The function describes the pattern
#' that local adaptation would leave in samples of individuals, not its demographic
#' consequences.
#'
#' @param x An `mc_sim` object from [simulate_metacommunity()].
#' @param n Number of individuals measured in each population, or all of them when the
#'   population is smaller. Ignored when `fraction` is given.
#' @param fraction Optional probability that each individual is measured, so that the number
#'   measured in a population is Binomial(abundance, `fraction`), proportional to abundance on
#'   average.
#' @param local_adaptation Fraction (0 to 1) of the gap between the optimum of the species and
#'   the environment of the site that each population closes; one value or one per species
#'   (named by species, or in their order). The same holds for the two standard deviations.
#' @param population_sd Standard deviation of the deviations of populations from the species
#'   mean that do not track the environment; one value or one per species.
#' @param intraspecific_sd Standard deviation of trait values among the individuals of a
#'   population; one value or one per species.
#' @param step Recorded step to sample, as in `record_steps`; `NULL` (the default) uses the
#'   abundances at the last step.
#' @param env Optional environment of each site at the sampled step. The default is the
#'   environment of the landscape; give it for simulations with `env_time`.
#' @return Data frame with one row per individual and the columns `site`, `species`,
#'   `individual` (its number within the population), `abundance` (of the population), `env`
#'   (of the site), `optimum` (of the species), `population_mean` (the expected value of the
#'   population, `theta_j + h_j (E_i - theta_j) + d_ij`) and `trait`. Populations with no
#'   individual measured are absent.
#' @references Violle, C., Enquist, B. J., McGill, B. J., Jiang, L., Albert, C. H., Hulshof,
#'   C., Jung, V., & Messier, J. (2012). The return of the variance: intraspecific variability
#'   in community ecology. *Trends in Ecology & Evolution*, 27, 244--252.
#'   \doi{10.1016/j.tree.2011.11.014}
#' @seealso [tstats()] for the T-statistics.
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(15)
#' sim <- simulate_metacommunity(landscape_implicit(10), species_traits(tree),
#'                               competition_matrix(tree, "none"), n_steps = 100,
#'                               seed_steps = 1:5)
#' ind <- sample_individuals(sim, n = 10, local_adaptation = 0.5, intraspecific_sd = 2)
#' head(ind)
#' @export
sample_individuals <- function(x, n = 10, fraction = NULL, local_adaptation = 0,
                               population_sd = 0, intraspecific_sd = 1, step = NULL,
                               env = NULL) {
  stopifnot(inherits(x, "mc_sim"))
  abundance <- sampled_abundance(x, step)
  sites <- rownames(abundance)
  species <- colnames(abundance)
  k <- length(species)
  if (is.null(env)) env <- x$landscape$env
  if (!is.numeric(env) || length(env) != length(sites) || anyNA(env)) {
    stop("`env` needs one number per site (", length(sites), ")", call. = FALSE)
  }
  h <- per_species(local_adaptation, species, "local_adaptation")
  if (any(h < 0 | h > 1)) stop("`local_adaptation` must be between 0 and 1", call. = FALSE)
  pop_sd <- per_species(population_sd, species, "population_sd")
  ind_sd <- per_species(intraspecific_sd, species, "intraspecific_sd")
  if (any(pop_sd < 0) || any(ind_sd < 0)) {
    stop("`population_sd` and `intraspecific_sd` must not be negative", call. = FALSE)
  }
  if (is.null(fraction)) {
    stopifnot(length(n) == 1, n >= 1, n == round(n))
    measured <- pmin(abundance, n)
  } else {
    stopifnot(length(fraction) == 1, fraction > 0, fraction <= 1)
    measured <- matrix(stats::rbinom(length(abundance), abundance, fraction),
                       nrow(abundance), dimnames = dimnames(abundance))
  }
  theta <- x$traits[species, "optimum"]
  pop_mean <- matrix(theta, length(sites), k, byrow = TRUE) +
    sweep(outer(env, theta, "-"), 2, h, "*") +
    sweep(matrix(stats::rnorm(length(sites) * k), length(sites), k), 2, pop_sd, "*")
  occupied <- which(measured > 0)
  site_of <- (occupied - 1) %% length(sites) + 1
  species_of <- (occupied - 1) %/% length(sites) + 1
  by_site <- order(site_of, species_of)
  occupied <- occupied[by_site]
  site_of <- site_of[by_site]
  species_of <- species_of[by_site]
  counts <- as.vector(measured)[occupied]
  rows <- rep(seq_along(occupied), counts)
  out <- data.frame(site = sites[site_of][rows], species = species[species_of][rows],
                    individual = sequence(counts),
                    abundance = as.vector(abundance)[occupied][rows],
                    env = env[site_of][rows], optimum = theta[species_of][rows],
                    population_mean = pop_mean[occupied][rows], stringsAsFactors = FALSE)
  out$trait <- out$population_mean + stats::rnorm(nrow(out)) * ind_sd[species_of][rows]
  out
}

#' T-statistics of intraspecific trait variation
#'
#' Computes the three variance ratios of Violle et al. (2012) for each site (community) and
#' trait, from trait values of individuals such as those of [sample_individuals()]:
#'
#' * `T_IP.IC`, the mean variance among individuals within the populations of the site over the
#'   variance among all individuals of the site. Low values mean that the populations of
#'   co-occurring species overlap little (niche differentiation within the community).
#' * `T_IC.IR`, the variance among individuals of the site over that among all individuals of
#'   the region. Low values mean that the site holds a narrow part of the regional trait range
#'   (environmental filtering of individuals).
#' * `T_PC.PR`, the variance among population means within the site over that among all
#'   population means of the region. Low values mean filtering of populations.
#'
#' A population is a species at a site. Variances use `n - 1`, means of the variances within
#' populations are not weighted, and populations of a single individual have no variance and
#' are left out of `T_IP.IC`, as in [cati::Tstats()].
#'
#' `method = "base"` computes the observed values in base R, with no null model. It is fast and
#' needs no other package, which suits simulation studies with many replicates, where the
#' truth is known and each statistic is compared among scenarios rather than with a null
#' model. `method = "cati"` calls [cati::Tstats()] (Taudiere & Violle 2016), which adds null
#' models: individuals permuted within each site for `T_IP.IC`, individuals drawn from the
#' regional pool for `T_IC.IR`, and population means permuted for `T_PC.PR`. Both give the
#' same observed values. [cati::Tstats()] finds the populations of a site by matching its name
#' as a regular expression, so a site called `site1` would also take the populations of
#' `site10`; the sites are therefore passed under labels of equal width (`site01`, `site02`,
#' ...) and given back their names in the result.
#'
#' @section Local adaptation:
#' With `local_adaptation` in [sample_individuals()], population means move toward the
#' environment of their site, so co-occurring species become more alike and the region spreads
#' along the environment. The variance among populations within a site falls, so `T_PC.PR` and
#' `T_IC.IR` fall, and `T_IP.IC` rises, because the populations that remain overlap more. A
#' fall in `T_IP.IC` with local adaptation would need populations to diverge from co-occurring
#' species, as in character displacement, which [sample_individuals()] does not simulate.
#'
#' @param individuals Data frame with one row per individual, such as from
#'   [sample_individuals()].
#' @param traits Names of the columns with trait values.
#' @param site,species Names of the columns with the site and the species of each individual.
#' @param method `"base"` (the default) for the observed values only, or `"cati"` for
#'   [cati::Tstats()] with null models; it needs the package cati.
#' @param nperm Number of randomisations of the null models of `method = "cati"`; `NULL` for
#'   observed values only.
#' @param ... Further arguments to [cati::Tstats()], such as `reg.pool` or
#'   `independantTraits`.
#' @return With `method = "base"`, a data frame with one row per site and trait and the columns
#'   `site`, `trait`, `n_populations`, `T_IP.IC`, `T_IC.IR` and `T_PC.PR`; statistics that
#'   cannot be computed, such as `T_PC.PR` for a site with one species, are `NA`. With
#'   `method = "cati"`, the `Tstats` object of [cati::Tstats()], with the site names of
#'   `individuals`, for its methods [cati::ses.Tstats()], `plot()` and `summary()`.
#' @references Taudiere, A., & Violle, C. (2016). cati: an R package using functional traits
#'   to detect and quantify multi-level community assembly processes. *Ecography*, 39,
#'   699--708. \doi{10.1111/ecog.01433}
#'
#'   Violle, C., Enquist, B. J., McGill, B. J., Jiang, L., Albert, C. H., Hulshof, C., Jung,
#'   V., & Messier, J. (2012). The return of the variance: intraspecific variability in
#'   community ecology. *Trends in Ecology & Evolution*, 27, 244--252.
#'   \doi{10.1016/j.tree.2011.11.014}
#' @seealso [sample_individuals()].
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(20)
#' sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
#'                               competition_matrix(tree, "none"), n_steps = 100,
#'                               seed_steps = 1:5)
#' none <- tstats(sample_individuals(sim, intraspecific_sd = 2))
#' local <- tstats(sample_individuals(sim, local_adaptation = 0.8, intraspecific_sd = 2))
#' colMeans(none[, c("T_IP.IC", "T_IC.IR", "T_PC.PR")], na.rm = TRUE)
#' colMeans(local[, c("T_IP.IC", "T_IC.IR", "T_PC.PR")], na.rm = TRUE)
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE) && requireNamespace("cati", quietly = TRUE)
#' # with the null models of cati; ses.Tstats() lists T_IP.IC, T_IC.IR and T_PC.PR in order
#' res <- tstats(sample_individuals(sim, intraspecific_sd = 2), method = "cati", nperm = 19)
#' round(cati::ses.Tstats(res)[[1]]$ses[, 1], 1)
#' @export
tstats <- function(individuals, traits = "trait", site = "site", species = "species",
                   method = c("base", "cati"), nperm = 99, ...) {
  method <- match.arg(method)
  stopifnot(is.data.frame(individuals), is.character(traits), length(traits) >= 1)
  missing_columns <- setdiff(c(traits, site, species), names(individuals))
  if (length(missing_columns) > 0) {
    stop("`individuals` has no column ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  values <- as.matrix(individuals[traits])
  if (!is.numeric(values)) stop("the columns in `traits` must be numeric", call. = FALSE)
  site_id <- as.character(individuals[[site]])
  species_id <- as.character(individuals[[species]])
  if (method == "cati") return(tstats_cati(values, site_id, species_id, nperm, ...))
  site_levels <- unique(site_id)
  population <- paste(site_id, species_id, sep = "\r")
  pop_site <- site_id[!duplicated(population)]
  out <- lapply(traits, function(trait) {
    z <- values[, trait]
    pop_mean <- tapply(z, population, mean, na.rm = TRUE)[unique(population)]
    pop_var <- tapply(z, population, stats::var, na.rm = TRUE)[unique(population)]
    var_ic <- tapply(z, site_id, stats::var, na.rm = TRUE)[site_levels]
    var_pc <- tapply(pop_mean, pop_site, stats::var, na.rm = TRUE)[site_levels]
    mean_ip <- tapply(pop_var, pop_site, mean, na.rm = TRUE)[site_levels]
    data.frame(site = site_levels, trait = trait,
               n_populations = as.vector(table(pop_site)[site_levels]),
               T_IP.IC = as.vector(mean_ip / var_ic),
               T_IC.IR = as.vector(var_ic / stats::var(z, na.rm = TRUE)),
               T_PC.PR = as.vector(var_pc / stats::var(pop_mean, na.rm = TRUE)),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  ratios <- c("T_IP.IC", "T_IC.IR", "T_PC.PR")
  out[ratios] <- lapply(out[ratios], function(v) ifelse(is.finite(v), v, NA_real_))
  rownames(out) <- NULL
  out
}

tstats_cati <- function(values, site_id, species_id, nperm, ...) {
  if (!cati_available()) {
    stop("`method = \"cati\"` needs the package cati: install.packages(\"cati\"). ",
         "`method = \"base\"` gives the observed values without null models.", call. = FALSE)
  }
  original <- unique(site_id)
  safe <- sprintf("site%0*d", nchar(length(original)), seq_along(original))
  site_map <- stats::setNames(original, safe)
  plots <- factor(safe[match(site_id, original)], levels = safe)
  res <- cati::Tstats(values, ind.plot = plots, sp = factor(species_id), nperm = nperm,
                      printprogress = FALSE, ...)
  restore_site_names(res, site_map)
}

cati_available <- function() requireNamespace("cati", quietly = TRUE)

# gives back the names of sites in every dimnames and factor of a cati result
restore_site_names <- function(x, site_map) {
  if (is.factor(x)) {
    if (all(levels(x) %in% names(site_map))) levels(x) <- unname(site_map[levels(x)])
    return(x)
  }
  if (is.list(x) && !is.data.frame(x)) {
    for (i in seq_along(x)) {
      if (!is.null(x[[i]])) x[[i]] <- restore_site_names(x[[i]], site_map)
    }
    return(x)
  }
  labels <- dimnames(x)
  if (!is.null(labels)) {
    dimnames(x) <- lapply(labels, function(d) {
      if (!is.null(d) && all(d %in% names(site_map))) unname(site_map[d]) else d
    })
  }
  x
}

sampled_abundance <- function(x, step) {
  if (is.null(step)) return(x$abundance)
  stopifnot(length(step) == 1)
  if (is.null(x$history)) stop("`x` has no recorded steps: use `record_steps`", call. = FALSE)
  label <- paste0("step", step)
  if (!label %in% dimnames(x$history)[[3]]) {
    stop("step ", step, " was not recorded", call. = FALSE)
  }
  matrix(x$history[, , label], dim(x$history)[1], dimnames = dimnames(x$history)[1:2])
}

per_species <- function(value, species, name) {
  if (!is.numeric(value) || anyNA(value)) stop("`", name, "` must be numbers", call. = FALSE)
  if (length(value) == 1 && is.null(names(value))) return(rep(value, length(species)))
  if (!is.null(names(value))) {
    if (!all(species %in% names(value))) {
      stop("`", name, "` has no value for ", species_list(setdiff(species, names(value))),
           call. = FALSE)
    }
    return(unname(value[species]))
  }
  if (length(value) != length(species)) {
    stop("`", name, "` needs one value or one per species (", length(species), ")",
         call. = FALSE)
  }
  value
}
