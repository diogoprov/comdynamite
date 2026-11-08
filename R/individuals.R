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
#' With `displacement > 0`, the populations of each site also diverge from one another, as in
#' ecological character displacement. Their means are the equilibrium of
#'
#' `mu_ij = m_ij + c sum_k p_ik (mu_ij - mu_ik) exp(-(mu_ij - mu_ik)^2 / (2 w^2))`,
#'
#' where `m_ij` is the mean above (`theta_j + h_j (E_i - theta_j) + d_ij`), `c` is
#' `displacement`, `p_ik` the share of species `k` among the individuals of the site and `w`
#' (`displacement_width`) the width of a Gaussian competition kernel. Each population is pushed
#' away from the populations close to it in trait, more so from abundant ones, and held near
#' `m_ij`, so `c` sets the strength of competition relative to the selection that keeps the
#' population near `m_ij`. The push from one competitor is largest at a distance `w` and at
#' most `c w exp(-1/2)`, about `0.6 c w`. The equilibrium is found by damped fixed-point
#' iteration from `m_ij`, with steps of `1 / (1 + c)`. Populations with the same `m_ij`, as
#' with `local_adaptation = 1` and no `population_sd`, are not pushed apart, since the push is
#' symmetric.
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
#' @param displacement Strength of character displacement among the populations of a site
#'   (`c` in Details), relative to the selection that keeps each population near its mean
#'   without competitors; 0 (the default) for none.
#' @param displacement_width Width `w` of the Gaussian competition kernel of character
#'   displacement, in units of the trait. The default is the mean niche breadth (`sigma`) of
#'   the species.
#' @param step Recorded step to sample, as in `record_steps`; `NULL` (the default) uses the
#'   abundances at the last step.
#' @param env Optional environment of each site at the sampled step. The default is the
#'   environment of the landscape; give it for simulations with `env_time`.
#' @return Data frame with one row per individual and the columns `site`, `species`,
#'   `individual` (its number within the population), `abundance` (of the population), `env`
#'   (of the site), `optimum` (of the species), `population_mean` (the expected value of the
#'   population, `theta_j + h_j (E_i - theta_j) + d_ij`, after character displacement) and
#'   `trait`. Populations with no
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
                               population_sd = 0, intraspecific_sd = 1, displacement = 0,
                               displacement_width = NULL, step = NULL, env = NULL) {
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
  stopifnot(is.numeric(displacement), length(displacement) == 1, displacement >= 0)
  if (displacement > 0) {
    width <- displacement_width %||% mean(x$traits$sigma)
    stopifnot(is.numeric(width), length(width) == 1, width > 0)
    for (i in seq_along(sites)) {
      present <- which(abundance[i, ] > 0)
      if (length(present) < 2) next
      pop_mean[i, present] <- displace(pop_mean[i, present], abundance[i, present],
                                       displacement, width)
    }
  }
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
#' are left out of `T_IP.IC`, as in [cati::Tstats()]. Individuals with a missing value of a
#' trait are left out for that trait.
#'
#' @section Null models:
#' With `nperm > 0`, `method = "base"` compares each statistic with a null distribution, from
#' the three null models of Taudiere & Violle (2016), implemented in base R from their
#' description:
#'
#' * `T_IP.IC`: trait values are permuted among the individuals of each site, which keeps the
#'   variance of the site and breaks the structure of populations within it;
#' * `T_IC.IR`: each site receives as many individuals as it has, drawn without replacement from
#'   all individuals of the region;
#' * `T_PC.PR`: each individual takes the mean of its population, and each site then receives
#'   as many of these values as it has individuals, drawn without replacement from the region;
#'   population means are computed again from the values drawn.
#'
#' The standardised effect size is `(observed - mean of null) / sd of null`, and the p-value
#' is `(1 + r) / (1 + nperm)`, where `r` is the number of null values at or below the observed
#' one (`alternative = "less"`), at or above it (`"greater"`), or twice the smaller of the two
#' (`"two.sided"`, at most 1). Low values of the three statistics are the usual hypotheses
#' (niche differentiation within communities for `T_IP.IC`, environmental filtering for the
#' others), so `"less"` tests them. Each trait is randomised on its own.
#'
#' `method = "cati"` calls [cati::Tstats()] (Taudiere & Violle 2016) instead and returns its
#' object. Both methods give the same observed values. Two problems of cati 0.99.6 affect it:
#' [cati::Tstats()] finds the populations of a site by matching its name as a regular
#' expression, so a site called `site1` would also take the populations of `site10`
#' (\url{https://github.com/adrientaudiere/cati/issues/11}); the sites are therefore passed
#' under labels of equal width (`site01`, `site02`, ...) and given back their names in the
#' result. And its null model for `T_IP.IC` puts together the permutations of different sites
#' in the wrong order, so that the null values of a site come from the individuals of other
#' sites; its standardised effect sizes for `T_IP.IC` are then wrong. `method = "base"` is
#' recommended.
#'
#' @section Local adaptation and character displacement:
#' With `local_adaptation` in [sample_individuals()], population means move toward the
#' environment of their site, so co-occurring species become more alike and the region spreads
#' along the environment. The variance among populations within a site falls, so `T_PC.PR` and
#' `T_IC.IR` fall, and `T_IP.IC` rises, because the populations that remain overlap more. With
#' `displacement`, co-occurring populations diverge from one another instead, which raises the
#' variance among populations within sites and lowers `T_IP.IC`.
#'
#' @param individuals Data frame with one row per individual, such as from
#'   [sample_individuals()].
#' @param traits Names of the columns with trait values.
#' @param site,species Names of the columns with the site and the species of each individual.
#' @param method `"base"` (the default), computed in base R, or `"cati"` for [cati::Tstats()];
#'   it needs the package cati.
#' @param nperm Number of randomisations of the null models; 0 or `NULL` for the observed
#'   values only, much faster in simulation studies with many replicates.
#' @param alternative Direction of the test of `method = "base"`: `"less"` (the default),
#'   `"greater"` or `"two.sided"`.
#' @param ... Further arguments to [cati::Tstats()], such as `reg.pool` or
#'   `independantTraits`.
#' @return With `method = "base"`, a data frame with one row per site and trait and the columns
#'   `site`, `trait`, `n_populations`, `T_IP.IC`, `T_IC.IR` and `T_PC.PR`; statistics that
#'   cannot be computed, such as `T_PC.PR` for a site with one species, are `NA`. With
#'   `nperm > 0`, also `ses_` and `p_` columns for each statistic, and the attribute `null`,
#'   an array (randomisation x row x statistic) with the null values. With `method = "cati"`,
#'   the `Tstats` object of [cati::Tstats()], with the site names of `individuals`.
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
#' ratios <- c("T_IP.IC", "T_IC.IR", "T_PC.PR")
#' none <- tstats(sample_individuals(sim, intraspecific_sd = 2), nperm = 0)
#' local <- tstats(sample_individuals(sim, local_adaptation = 0.8, intraspecific_sd = 2),
#'                 nperm = 0)
#' apart <- tstats(sample_individuals(sim, displacement = 2, intraspecific_sd = 2), nperm = 0)
#' rbind(none = colMeans(none[ratios], na.rm = TRUE),
#'       local = colMeans(local[ratios], na.rm = TRUE),
#'       apart = colMeans(apart[ratios], na.rm = TRUE))
#' # standardised effect sizes and p-values from the null models
#' res <- tstats(sample_individuals(sim, intraspecific_sd = 2), nperm = 99)
#' head(res[c("site", "T_IP.IC", "ses_T_IP.IC", "p_T_IP.IC")])
#' @export
tstats <- function(individuals, traits = "trait", site = "site", species = "species",
                   method = c("base", "cati"), nperm = 99,
                   alternative = c("less", "greater", "two.sided"), ...) {
  method <- match.arg(method)
  alternative <- match.arg(alternative)
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
  nperm <- if (is.null(nperm)) 0 else nperm
  stopifnot(length(nperm) == 1, nperm >= 0, nperm == round(nperm))
  ratios <- c("T_IP.IC", "T_IC.IR", "T_PC.PR")
  site_levels <- unique(site_id)
  per_trait <- lapply(traits, function(trait) {
    keep <- !is.na(values[, trait])
    z <- values[keep, trait]
    groups <- tstats_groups(site_id[keep], species_id[keep], site_levels)
    observed <- t_ratios(z, groups)
    null <- if (nperm > 0) tstats_null(z, groups, nperm) else NULL
    list(table = data.frame(site = site_levels, trait = trait,
                            n_populations = tabulate(groups$pop_site, length(site_levels)),
                            observed, stringsAsFactors = FALSE),
         null = null)
  })
  out <- do.call(rbind, lapply(per_trait, `[[`, "table"))
  rownames(out) <- NULL
  out[ratios] <- lapply(out[ratios], finite_or_na)
  if (nperm == 0) return(out)
  null <- do.call(abind_rows, lapply(per_trait, `[[`, "null"))
  for (stat in ratios) {
    sims <- null[, , stat, drop = FALSE][, , 1]
    sims <- matrix(sims, nperm)
    obs <- out[[stat]]
    valid <- colSums(is.finite(sims))
    null_mean <- colMeans(sims, na.rm = TRUE)
    null_sd <- apply(sims, 2, stats::sd, na.rm = TRUE)
    below <- colSums(sims <= rep(obs, each = nperm), na.rm = TRUE)
    above <- colSums(sims >= rep(obs, each = nperm), na.rm = TRUE)
    p <- switch(alternative, less = (1 + below) / (1 + valid),
                greater = (1 + above) / (1 + valid),
                two.sided = pmin(1, 2 * pmin(1 + below, 1 + above) / (1 + valid)))
    out[[paste0("ses_", stat)]] <- finite_or_na((obs - null_mean) / null_sd)
    out[[paste0("p_", stat)]] <- ifelse(is.na(obs) | valid == 0, NA_real_, p)
  }
  attr(out, "null") <- null
  out
}

# integer codes of sites, populations and the site of each population
tstats_groups <- function(site_id, species_id, site_levels) {
  population <- paste(match(site_id, site_levels), species_id, sep = "\r")
  pop_levels <- unique(population)
  site <- match(site_id, site_levels)
  pop <- match(population, pop_levels)
  list(site = site, pop = pop, pop_site = site[match(seq_along(pop_levels), pop)],
       n_sites = length(site_levels), n_pops = length(pop_levels),
       by_site = split(seq_along(site), factor(site, seq_along(site_levels))))
}

# the three ratios per site, from trait values and the groups of tstats_groups()
t_ratios <- function(z, g) {
  z <- z - mean(z)
  pop_mean <- group_mean(z, g$pop, g$n_pops)
  pop_var <- group_var(z, g$pop, g$n_pops)
  var_ic <- group_var(z, g$site, g$n_sites)
  has_var <- !is.na(pop_var)
  mean_ip <- group_sum(ifelse(has_var, pop_var, 0), g$pop_site, g$n_sites) /
    group_sum(as.numeric(has_var), g$pop_site, g$n_sites)
  var_pc <- group_var(pop_mean, g$pop_site, g$n_sites)
  data.frame(T_IP.IC = mean_ip / var_ic, T_IC.IR = var_ic / stats::var(z),
             T_PC.PR = var_pc / stats::var(pop_mean))
}

# null values (randomisation x site x statistic) of the three null models
tstats_null <- function(z, g, nperm) {
  n <- length(z)
  positions <- unlist(g$by_site, use.names = FALSE)
  sizes <- lengths(g$by_site)
  carried <- group_mean(z, g$pop, g$n_pops)[g$pop]
  out <- array(NA_real_, c(nperm, g$n_sites, 3),
               list(NULL, NULL, c("T_IP.IC", "T_IC.IR", "T_PC.PR")))
  for (i in seq_len(nperm)) {
    within <- z
    for (ix in g$by_site) within[ix] <- z[ix[sample.int(length(ix))]]
    out[i, , "T_IP.IC"] <- t_ratios(within, g)$T_IP.IC
    regional <- numeric(n)
    regional[positions] <- unlist(lapply(sizes, function(m) z[sample.int(n, m)]))
    out[i, , "T_IC.IR"] <- t_ratios(regional, g)$T_IC.IR
    populations <- numeric(n)
    populations[positions] <- unlist(lapply(sizes, function(m) carried[sample.int(n, m)]))
    out[i, , "T_PC.PR"] <- t_ratios(populations, g)$T_PC.PR
  }
  out
}

group_sum <- function(x, g, n_groups) {
  as.vector(rowsum(x, factor(g, seq_len(n_groups)), reorder = TRUE))
}

group_mean <- function(x, g, n_groups) group_sum(x, g, n_groups) / tabulate(g, n_groups)

group_var <- function(x, g, n_groups) {
  n <- tabulate(g, n_groups)
  mean_g <- group_sum(x, g, n_groups) / n
  out <- group_sum((x - mean_g[g])^2, g, n_groups) / (n - 1)
  out[n < 2] <- NA_real_
  out
}

finite_or_na <- function(v) ifelse(is.finite(v), v, NA_real_)

# binds arrays (randomisation x row x statistic) along their rows
abind_rows <- function(...) {
  parts <- list(...)
  if (length(parts) == 1) return(parts[[1]])
  dims <- dim(parts[[1]])
  rows <- sum(vapply(parts, function(p) dim(p)[2], numeric(1)))
  out <- array(NA_real_, c(dims[1], rows, dims[3]), list(NULL, NULL, dimnames(parts[[1]])[[3]]))
  start <- 0
  for (p in parts) {
    out[, start + seq_len(dim(p)[2]), ] <- p
    start <- start + dim(p)[2]
  }
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

# equilibrium of mu_j = target_j + c sum_k p_k (mu_j - mu_k) exp(-(mu_j - mu_k)^2 / (2 w^2)),
# found by damped fixed-point iteration from the targets; the step keeps the iteration stable,
# since the slope of the push is between -0.45 c and c
displace <- function(target, abundance, strength, width, tol = 1e-8, max_iter = 1e5) {
  p <- abundance / sum(abundance)
  step <- 1 / (1 + strength)
  mu <- target
  for (iter in seq_len(max_iter)) {
    gap <- outer(mu, mu, "-")
    push <- strength * as.vector((gap * exp(-gap^2 / (2 * width^2))) %*% p)
    change <- step * (target + push - mu)
    mu <- mu + change
    if (max(abs(change)) < tol * width * step) return(mu)
  }
  warning("character displacement did not converge at a site; ",
          "a smaller `displacement` may help", call. = FALSE)
  mu
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
