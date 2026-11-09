# Estimation of local adaptation and character displacement from the means of populations,
# with diagnostics of the width of the competition kernel and a parametric bootstrap.

#' Local adaptation and character displacement from population means
#'
#' Estimates how far the populations of a species follow the environment of their site and
#' how far co-occurring species push them apart, by fitting both terms of the model of
#' [sample_individuals()] at once to the means of the populations:
#'
#' \deqn{m_{ij} = a_j + h E_i + c \sum_k p_{ik} (m_{ij} - m_{ik})
#'   \exp\left(-\frac{(m_{ij} - m_{ik})^2}{2 w^2}\right),}{m_ij = a_j + h E_i +
#'   c sum_k p_ik (m_ij - m_ik) exp(-(m_ij - m_ik)^2 / (2 w^2)),}
#'
#' where \eqn{m_{ij}}{m_ij} is the mean of species \eqn{j} at site \eqn{i}, \eqn{a_j}{a_j} a
#' constant of the species, \eqn{E_i}{E_i} the environment of the site and
#' \eqn{p_{ik}}{p_ik} the relative abundance of species \eqn{k} at the site. In
#' [sample_individuals()], \eqn{a_j = (1 - h) \theta_j}{a_j = (1 - h) theta_j}, `h` is
#' `local_adaptation`, `c` is `displacement` and `w` is `displacement_width`. The constants of
#' the species are removed by centring each variable within species, and `h` and `c` are the
#' coefficients of a regression through the origin; `w`, when not given, is the value of a
#' grid that leaves the smallest residual sum of squares.
#'
#' Local adaptation alone moves the populations of a species along the environment, and
#' displacement alone moves them away from the species they meet. The two leave different
#' signatures, but each term alone picks up part of the other: the slope of population means on
#' the environment falls when displacement pulls co-occurring populations apart, and the
#' spread of a site rises with displacement and falls with local adaptation. Fitting both
#' terms together separates them when populations of the same species meet different
#' environments and different competitors. The approach follows the advice of Goldberg & Lande
#' (2006) and Adams & Collyer (2007) to judge divergence among sympatric species after the
#' variation of the trait along the environment is accounted for, here with the strength of
#' competition at each site in place of a contrast between sympatry and allopatry.
#'
#' With several columns in `env`, each environmental variable gets its own coefficient
#' (\eqn{h E_i} becomes \eqn{\sum_v h_v E_{iv}}{sum_v h_v E_iv}). With `species_slopes = TRUE`
#' each species gets its own slope on each variable; the table then gives their unweighted mean
#' and the attribute `species_slopes` the slope of each species.
#'
#' @section Competitors that were not measured:
#' By default the competitors of a population are the other populations measured at its site,
#' weighted by `abundance`. Species that occur at a site but were not measured are then left
#' out of the push and bias `c`. With `community`, a matrix of abundances (sites in rows,
#' species in columns, with names as in `individuals`), every species present at a site pushes,
#' with its share of the abundance of the site. A species not measured at a site enters with
#' the mean of its populations measured elsewhere, moved along the gradient by the fitted slopes
#' times the difference between the environment of the site and the mean environment of the
#' sites where the species was measured (the model is fitted once with the plain means, then
#' again with the moved ones). Species never measured have no trait value and are left out,
#' with a warning, of both the push and the total abundance of the site that sets the shares.
#' Populations measured at a site where `community` gives
#' them a smaller abundance than the number of individuals measured (such as zero, when the
#' census and the measurements come from different visits) take that number as abundance.
#'
#' @section Sampling error:
#' The mean of a population enters both sides of the equation, so its sampling error would
#' correlate the push of competitors with the response and inflate `c`. The individuals of each
#' population are therefore split at random in two halves (so the estimates change slightly
#' with the seed of the random number generator): the response is the
#' mean of one half, and the own mean in the push is that of the other half (competitors enter
#' with their full means). Both assignments are stacked in one fit. Only populations with at
#' least two individuals, of species with at least two such populations, enter as responses;
#' all populations enter as competitors. Sampling error left in the push still biases `c`
#' toward zero, more so with few individuals per population and with strong displacement:
#' strongly displaced populations sit in the tail of the kernel of their competitors, where the
#' push is small but its sampling error is not.
#'
#' Deviations of whole populations from the model (from genetic drift or unmeasured
#' environmental variables) are not removed by the split: they enter the own mean in the push
#' and the response alike, and bias `c` downward, below zero when there is no displacement. The
#' bootstrap reproduces this bias, so `p_displacement` and `displacement_corrected` account for
#' it, but the raw `displacement` does not.
#'
#' @section Width of the kernel:
#' Without displacement the width is not identified: the residual sum of squares changes little
#' along the grid, and a narrow kernel, which acts on few pairs of populations, can give a value
#' of `c` far from zero, of either sign. The grid runs from the larger of one eighth of the
#' median standard deviation of population means within sites and the pooled standard
#' deviation within populations (a kernel narrower than the spread of individuals mostly fits
#' sampling noise), but no higher than that median, to four times the median, and is refined by
#' [stats::optimize()] around its best value. A warning is given when the best width stays at
#' one of the ends of the grid after the refinement or when the residual sum of squares changes
#' by less than 5% along the grid; `plot()` shows the profile of
#' the fit and of the estimates along the width. Fix `width` from independent knowledge of the
#' niche widths when it is not identified.
#'
#' @section Bootstrap:
#' With `nboot > 0`, new data are simulated with the design of the data (the same sites,
#' species, environment, abundances and numbers of individuals per population), with the
#' within-population standard deviation of the data, a population-level standard deviation
#' estimated from the residuals of the fit without displacement, and the fitted species
#' constants and slopes, and refitted as the data were (with the same `width`, estimated or
#' fixed). Two sets are simulated:
#' * without displacement (`c = 0`), which gives `p_displacement`, the two-sided p-value of
#'   `c` with the width chosen anew in each simulation, so that the choice of the width is
#'   accounted for. The p-value is the share of simulations at least as far from the mean of
#'   the simulations as the estimate, since that mean is not zero when sampling error or
#'   deviations of populations bias `c`;
#' * with the fitted values, \eqn{c_0 = \max(\hat c, 0)}{c_0 = max(c_hat, 0)} (the model of
#'   [sample_individuals()] has no attraction), which give a first-order bias correction of
#'   `c`, `displacement_corrected` (\eqn{\hat c - (\bar c^* - c_0)}{c_hat - (mean(c*) - c_0)}),
#'   and basic bootstrap intervals for `c` and for the slopes (\eqn{\hat\theta -
#'   (q^*_{1 - \alpha/2} - \theta_0)}{theta_hat - (q*(1 - alpha/2) - theta_0)} to
#'   \eqn{\hat\theta - (q^*_{\alpha/2} - \theta_0)}{theta_hat - (q*(alpha/2) - theta_0)}, with
#'   \eqn{\theta_0}{theta_0} the value that generated the simulations).
#'
#' The correction and the intervals assume that the bias does not change with `c` near the
#' estimate; it does change with strong displacement and with deviations of whole populations,
#' so the intervals of `c` cover the true value less often than `level`. The correction removes
#' most of the attenuation of `c` by sampling error when `w` is well identified; when it is not,
#' the intervals are wide and the p-value is the better guide. Results change with the seed of
#' the random number generator. The simulated values are in the attribute `bootstrap`.
#'
#' @param individuals Data frame with one row per individual, such as from
#'   [sample_individuals()].
#' @param trait,site,species Names of the columns with the trait value, the site and the
#'   species.
#' @param env Names of one or more columns with the environment of the site.
#' @param abundance Name of the column with the abundance of the population; not used with
#'   `community`.
#' @param width Width `w` of the competition kernel; `NULL` (the default) chooses it from the
#'   grid described in the section on the width. `c` is estimated at that width, and the two
#'   trade off when data are few.
#' @param species_slopes Give each species its own slope on each environmental variable.
#' @param community Optional matrix or data frame of abundances, sites in rows and species in
#'   columns, so that species present but not measured also push (see the section on
#'   competitors).
#' @param nboot Number of parametric bootstrap simulations (0, the default, for none).
#' @param level Confidence level of the bootstrap intervals.
#' @return An object of class `population_shifts`: a data frame with one row and the columns
#'   `local_adaptation` (`h`; with several variables, `local_adaptation_<variable>`),
#'   `displacement` (`c`), `width` (`w`), `r_squared` (of the fit within species),
#'   `n_populations` and `n_species` (that entered as responses). With `nboot > 0` it also has
#'   `displacement_corrected`, `displacement_low`, `displacement_high`, `p_displacement` and,
#'   for each slope, `_low` and `_high` columns. The attribute `profile` holds, for each width
#'   tried, the residual sum of squares, \eqn{R^2} and the estimates.
#' @references Adams, D. C., & Collyer, M. L. (2007). Analysis of character divergence along
#'   environmental gradients and other covariates. *Evolution*, 61, 510--515.
#'   \doi{10.1111/j.1558-5646.2007.00063.x}
#'
#'   Goldberg, E. E., & Lande, R. (2006). Ecological and reproductive character displacement on
#'   an environmental gradient. *Evolution*, 60, 1344--1357. \doi{10.1554/05-696.1}
#' @seealso [sample_individuals()], [tstats()] with `null_ip = "species"`.
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(20)
#' sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
#'                               competition_matrix(tree, "none"), n_steps = 100,
#'                               seed_steps = 1:5)
#' ind <- sample_individuals(sim, n = 20, local_adaptation = 0.5, displacement = 2)
#' fit <- population_shifts(ind, nboot = 49)
#' fit
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(fit)
#' @export
population_shifts <- function(individuals, trait = "trait", site = "site", species = "species",
                              env = "env", abundance = "abundance", width = NULL,
                              species_slopes = FALSE, community = NULL, nboot = 0,
                              level = 0.95) {
  stopifnot(is.data.frame(individuals), is.character(env), length(env) >= 1,
            !anyDuplicated(env), is.logical(species_slopes), length(species_slopes) == 1,
            !is.na(species_slopes), is.numeric(nboot), length(nboot) == 1, nboot >= 0,
            nboot == round(nboot), is.numeric(level), length(level) == 1, level > 0,
            level < 1)
  if (!is.null(width)) stopifnot(is.numeric(width), length(width) >= 1, all(width > 0))
  if (is.null(community) && is.null(abundance)) {
    stop("give the `abundance` column or a `community` matrix", call. = FALSE)
  }
  columns <- c(trait, site, species, env, if (is.null(community)) abundance)
  missing_columns <- setdiff(columns, names(individuals))
  if (length(missing_columns) > 0) {
    stop("`individuals` has no column ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  d <- individuals[stats::complete.cases(individuals[columns]), columns]
  env_cols <- paste0(".env", seq_along(env))
  names(d) <- c("trait", "site", "species", env_cols, if (is.null(community)) "abundance")
  numeric_cols <- c("trait", env_cols, if (is.null(community)) "abundance")
  if (!all(vapply(d[numeric_cols], is.numeric, logical(1)))) {
    stop("trait, environment and abundance must be numeric", call. = FALSE)
  }
  d$site <- as.character(d$site)
  d$species <- as.character(d$species)
  if (!is.null(community)) community <- check_community(community, d)
  slope_names <- if (length(env) == 1) "local_adaptation" else
    paste0("local_adaptation_", env)
  prep <- shifts_prepare(d, env_cols, community, warn = TRUE)
  fit <- shifts_fit_full(prep, width, species_slopes, slope_names)
  prep <- fit$prep
  if (is.null(width) && !is.na(fit$table$displacement) && (fit$boundary || fit$flat)) {
    why <- c(if (fit$boundary) "the best width is at an end of the grid",
             if (fit$flat) sprintf("the residual sum of squares changes by %.1f%% along the grid",
                                   100 * fit$rss_change))
    warning(paste(why, collapse = " and "), ": `w` and `c` are poorly identified; see ",
            "plot(), give `width` or use `nboot`", call. = FALSE)
  }
  out <- fit$table
  if (nboot > 0) out <- shifts_bootstrap(out, fit, prep, d, env_cols, community, width,
                                         species_slopes, slope_names, nboot, level)
  attr(out, "profile") <- fit$profile
  if (species_slopes) attr(out, "species_slopes") <- fit$species_slopes
  if (nboot > 0) attr(out, "bootstrap") <- attr(out, "draws")
  attr(out, "draws") <- NULL
  if (is.null(width) && isTRUE(prep$within_sd > prep$spread / 8)) {
    attr(out, "floor") <- prep$floor
  }
  class(out) <- c("population_shifts", class(out))
  out
}

# checks the community matrix and returns it as a numeric matrix
check_community <- function(community, d) {
  community <- as.matrix(community)
  if (!is.numeric(community) || anyNA(community) || any(community < 0)) {
    stop("`community` must hold non-negative numbers without missing values", call. = FALSE)
  }
  if (is.null(rownames(community)) || is.null(colnames(community))) {
    stop("`community` needs site names as row names and species names as column names",
         call. = FALSE)
  }
  absent <- setdiff(unique(d$site), rownames(community))
  if (length(absent) > 0) {
    stop("`community` has no row for site ", paste(absent, collapse = ", "), call. = FALSE)
  }
  absent <- setdiff(unique(d$species), colnames(community))
  if (length(absent) > 0) {
    stop("`community` has no column for species ", paste(absent, collapse = ", "),
         call. = FALSE)
  }
  community
}

# populations, the two random halves, the competitors of each site and the scales of the trait
shifts_prepare <- function(d, env_cols, community, warn = FALSE) {
  label <- paste(d$site, d$species, sep = "\r")
  first <- !duplicated(label)
  pops <- d[first, c("site", "species", env_cols), drop = FALSE]
  pid <- match(label, label[first])
  n_pops <- nrow(pops)
  pops$n <- tabulate(pid, n_pops)
  pops$mean <- group_mean(d$trait, pid, n_pops)
  odd <- stats::ave(stats::runif(length(pid)), pid, FUN = function(u) rank(u) %% 2) == 1
  halves <- cbind(group_mean(d$trait[odd], pid[odd], n_pops),
                  group_mean(d$trait[!odd], pid[!odd], n_pops))
  deviation <- d$trait - pops$mean[pid]
  pooled <- pops$n[pid] > 1
  within_sd <- if (any(pooled)) {
    sqrt(sum(deviation[pooled]^2) / (sum(pooled) - sum(pops$n > 1)))
  } else NA_real_
  species_mean <- tapply(pops$mean, pops$species, mean)
  # mean environment of the sites where each species was measured
  species_env <- matrix(vapply(env_cols, function(v) {
    as.vector(tapply(pops[[v]], pops$species, mean)[names(species_mean)])
  }, numeric(length(species_mean))), length(species_mean),
  dimnames = list(names(species_mean), env_cols))
  sites <- unique(pops$site)
  members <- lapply(sites, function(s) {
    ix <- which(pops$site == s)
    if (is.null(community)) {
      return(list(ix = ix, share = d$abundance[first][ix], ext_mean = numeric(0),
                  ext_base = numeric(0), ext_gap = matrix(0, 0, length(env_cols)),
                  ext_share = numeric(0)))
    }
    row <- community[s, , drop = FALSE][1, ]
    measured <- pmax(row[pops$species[ix]], pops$n[ix])
    others <- setdiff(names(row)[row > 0], pops$species[ix])
    others <- others[others %in% names(species_mean)]
    # environment of the site minus the mean environment of the species elsewhere, so that
    # the mean of an unmeasured competitor can follow the gradient (see shifts_follow())
    site_env <- unlist(pops[ix[1], env_cols])
    gap <- matrix(rep(site_env, each = length(others)), length(others), length(env_cols)) -
      species_env[others, , drop = FALSE]
    list(ix = ix, share = unname(measured), ext_mean = unname(species_mean[others]),
         ext_base = unname(species_mean[others]), ext_gap = unname(gap),
         ext_share = unname(row[others]), ext_species = others)
  })
  dropped <- character(0)
  if (!is.null(community)) {
    present <- colnames(community)[colSums(community[sites, , drop = FALSE]) > 0]
    dropped <- setdiff(present, names(species_mean))
    if (warn && length(dropped) > 0) {
      warning("species in `community` that were never measured are left out of the push: ",
              paste(dropped, collapse = ", "), call. = FALSE)
    }
  }
  for (k in seq_along(members)) {
    total <- sum(members[[k]]$share) + sum(members[[k]]$ext_share)
    members[[k]]$share <- members[[k]]$share / total
    members[[k]]$ext_share <- members[[k]]$ext_share / total
  }
  spread <- stats::median(vapply(members, function(m) {
    if (length(m$ix) < 2) NA_real_ else stats::sd(pops$mean[m$ix])
  }, numeric(1)), na.rm = TRUE)
  if (!is.finite(spread) || spread <= 0) spread <- stats::sd(pops$mean)
  floor <- min(max(spread / 8, if (is.finite(within_sd)) within_sd else 0), spread)
  use <- pops$n >= 2
  use <- use & stats::ave(use, pops$species, FUN = sum) >= 2
  if (!any(use)) {
    stop("no species has two populations with at least two individuals", call. = FALSE)
  }
  list(pops = pops, halves = halves, members = members, use = use, env_cols = env_cols,
       within_sd = within_sd, spread = spread, floor = floor, n_pops = n_pops)
}

# push of the competitors of each population, with `own` as the population's own mean
shifts_push <- function(prep, own, w) {
  out <- numeric(prep$n_pops)
  means <- prep$pops$mean
  for (m in prep$members) {
    ix <- m$ix
    if (length(ix) + length(m$ext_mean) < 2) next
    gap <- outer(own[ix], means[ix], "-")
    diag(gap) <- 0
    push <- as.vector((gap * exp(-gap^2 / (2 * w^2))) %*% m$share)
    if (length(m$ext_mean) > 0) {
      ext <- outer(own[ix], m$ext_mean, "-")
      push <- push + as.vector((ext * exp(-ext^2 / (2 * w^2))) %*% m$ext_share)
    }
    out[ix] <- push
  }
  out
}

# moves the means of unmeasured competitors along the gradient with the fitted slopes
shifts_follow <- function(prep, slopes) {
  for (k in seq_along(prep$members)) {
    m <- prep$members[[k]]
    if (length(m$ext_base) > 0) {
      prep$members[[k]]$ext_mean <- m$ext_base + as.vector(m$ext_gap %*% slopes)
    }
  }
  prep
}

# fits the model; with unmeasured competitors, fits again after moving their means along the
# gradient with the slopes of the first fit
shifts_fit_full <- function(prep, width, species_slopes, slope_names) {
  fit <- shifts_fit(prep, width, species_slopes, slope_names)
  external <- any(vapply(prep$members, function(m) length(m$ext_base) > 0, logical(1)))
  slopes <- unlist(fit$table[slope_names])
  if (external && all(is.finite(slopes))) {
    prep <- shifts_follow(prep, slopes)
    fit <- shifts_fit(prep, width, species_slopes, slope_names)
  }
  fit$prep <- prep
  fit
}

# fits the model over the grid of widths
shifts_fit <- function(prep, width, species_slopes, slope_names) {
  use <- prep$use
  pops <- prep$pops
  sp <- rep(pops$species[use], 2)
  centre <- function(v) v - stats::ave(v, sp)
  y <- centre(c(prep$halves[use, 1], prep$halves[use, 2]))
  e <- vapply(prep$env_cols, function(v) centre(rep(pops[[v]][use], 2)), numeric(length(y)))
  e <- matrix(e, length(y))
  colnames(e) <- slope_names
  slope_species <- sort(unique(sp))
  if (species_slopes) {
    e <- do.call(cbind, lapply(seq_along(slope_names), function(v) {
      x <- vapply(slope_species, function(s) e[, v] * (sp == s), numeric(length(y)))
      x <- matrix(x, length(y))
      colnames(x) <- paste(slope_names[v], slope_species, sep = "\r")
      x
    }))
  }
  fit_width <- function(w) {
    push <- centre(c(shifts_push(prep, prep$halves[, 2], w)[use],
                     shifts_push(prep, prep$halves[, 1], w)[use]))
    x <- cbind(e, displacement = push)
    keep <- colSums(x^2) > 0
    coef <- stats::setNames(rep(NA_real_, ncol(x)), colnames(x))
    fit <- stats::lm.fit(x[, keep, drop = FALSE], y)
    coef[keep] <- fit$coefficients
    list(coef = coef, rss = sum(fit$residuals^2))
  }
  grid <- if (is.null(width)) {
    exp(seq(log(prep$floor), log(4 * prep$spread), length.out = 21))
  } else width
  fits <- lapply(grid, fit_width)
  rss <- vapply(fits, `[[`, numeric(1), "rss")
  boundary <- FALSE
  flat <- FALSE
  rss_change <- NA_real_
  if (is.null(width)) {
    at <- which.min(rss)
    boundary <- at %in% c(1, length(grid))
    rss_change <- (max(rss) - min(rss)) / min(rss)
    flat <- isTRUE(rss_change < 0.05)
    bounds <- log(grid[c(max(at - 1, 1), min(at + 1, length(grid)))])
    if (bounds[1] < bounds[2]) {
      refined <- exp(stats::optimize(function(lw) fit_width(exp(lw))$rss, bounds)$minimum)
    } else refined <- grid[at]
    # still at the end when the refinement did not move away from the end of the grid
    step <- diff(log(grid[1:2]))
    ends <- log(grid[c(1, length(grid))])
    boundary <- boundary && min(abs(log(refined) - ends)) < step / 10
    grid <- c(grid, refined)
    fits <- c(fits, list(fit_width(refined)))
    rss <- c(rss, fits[[length(fits)]]$rss)
  }
  best <- which.min(rss)
  slopes_of <- function(coef) {
    if (!species_slopes) return(coef[slope_names])
    vapply(slope_names, function(v) {
      mean(coef[startsWith(names(coef), paste0(v, "\r"))], na.rm = TRUE)
    }, numeric(1))
  }
  total <- sum(y^2)
  profile <- data.frame(width = grid, rss = rss, r_squared = 1 - rss / total,
                        displacement = vapply(fits, function(f) f$coef[["displacement"]],
                                              numeric(1)))
  profile <- cbind(profile, do.call(rbind, lapply(fits, function(f) slopes_of(f$coef))))
  coef <- fits[[best]]$coef
  table <- data.frame(as.list(slopes_of(coef)), displacement = coef[["displacement"]],
                      width = grid[best], r_squared = 1 - rss[best] / total,
                      n_populations = sum(use), n_species = length(slope_species),
                      check.names = FALSE)
  species_table <- NULL
  if (species_slopes) {
    species_table <- do.call(rbind, lapply(slope_names, function(v) {
      data.frame(species = slope_species, variable = v,
                 slope = unname(coef[paste(v, slope_species, sep = "\r")]))
    }))
  }
  list(table = table, coef = coef, profile = profile, boundary = boundary, flat = flat,
       rss_change = rss_change, species_slopes = species_table, slope_species = slope_species)
}

# parametric bootstrap with the design of the data
shifts_bootstrap <- function(out, fit, prep, d, env_cols, community, width, species_slopes,
                             slope_names, nboot, level) {
  pops <- prep$pops
  env <- as.matrix(pops[env_cols])
  # slopes of each population's species (the mean slope for species without their own)
  slope_matrix <- matrix(unlist(out[slope_names]), nrow(pops), length(slope_names),
                         byrow = TRUE)
  if (species_slopes) {
    for (v in seq_along(slope_names)) {
      s <- fit$coef[paste(slope_names[v], pops$species, sep = "\r")]
      slope_matrix[!is.na(s), v] <- s[!is.na(s)]
    }
  }
  along <- rowSums(env * slope_matrix)
  constant <- stats::ave(pops$mean - along, pops$species)
  residual <- (pops$mean - along - constant)[prep$use]
  n_used <- pops$n[prep$use]
  sampling <- if (is.finite(prep$within_sd)) prep$within_sd^2 / n_used else 0
  df <- max(sum(prep$use) - length(fit$slope_species) - length(slope_names), 1)
  pop_sd <- sqrt(max(0, sum(residual^2) / df - mean(sampling)))
  within_sd <- if (is.finite(prep$within_sd)) prep$within_sd else 0
  c_hat <- out$displacement
  w_hat <- out$width
  simulate <- function(strength) {
    target <- constant + along + stats::rnorm(nrow(pops), 0, pop_sd)
    m <- target
    if (strength > 0) {
      for (mb in prep$members) {
        k <- length(mb$ix)
        if (k + length(mb$ext_mean) < 2) next
        # competitors that were not measured are displaced too, but only measured ones are kept
        moved <- suppressWarnings(displace(c(target[mb$ix], mb$ext_mean),
                                           c(mb$share, mb$ext_share), strength, w_hat))
        m[mb$ix] <- moved[seq_len(k)]
      }
    }
    rows <- rep(seq_len(nrow(pops)), pops$n)
    sim <- data.frame(trait = m[rows] + stats::rnorm(length(rows), 0, within_sd),
                      site = pops$site[rows], species = pops$species[rows],
                      stringsAsFactors = FALSE)
    for (v in env_cols) sim[[v]] <- pops[[v]][rows]
    if (is.null(community)) {
      sim$abundance <- d$abundance[!duplicated(paste(d$site, d$species, sep = "\r"))][rows]
    }
    sim_prep <- shifts_prepare(sim, env_cols, community)
    unlist(shifts_fit_full(sim_prep, width, species_slopes,
                           slope_names)$table[c(slope_names, "displacement")])
  }
  draws <- function(strength) {
    t(vapply(seq_len(nboot), function(b) simulate(strength), numeric(length(slope_names) + 1)))
  }
  null <- draws(0)
  fitted_strength <- if (is.na(c_hat)) 0 else max(c_hat, 0)
  fitted <- if (fitted_strength > 0) draws(fitted_strength) else null
  alpha <- 1 - level
  # basic interval: the estimate minus the error of the simulations about the value that
  # generated them
  basic <- function(estimate, values, truth = estimate) {
    q <- stats::quantile(values, c(1 - alpha / 2, alpha / 2), na.rm = TRUE, names = FALSE)
    estimate - (q - truth)
  }
  if (is.na(c_hat)) {
    out$displacement_corrected <- out$displacement_low <- out$displacement_high <- NA_real_
    out$p_displacement <- NA_real_
  } else {
    c_null <- null[, "displacement"]
    c_fitted <- fitted[, "displacement"]
    out$displacement_corrected <- c_hat - (mean(c_fitted, na.rm = TRUE) - fitted_strength)
    interval <- basic(c_hat, c_fitted, fitted_strength)
    out$displacement_low <- interval[1]
    out$displacement_high <- interval[2]
    # centred on the mean of the null draws, which need not be zero
    centre <- mean(c_null, na.rm = TRUE)
    out$p_displacement <- (1 + sum(abs(c_null - centre) >= abs(c_hat - centre),
                                   na.rm = TRUE)) / (1 + sum(!is.na(c_null)))
  }
  for (v in slope_names) {
    interval <- basic(out[[v]], fitted[, v])
    out[[paste0(v, "_low")]] <- interval[1]
    out[[paste0(v, "_high")]] <- interval[2]
  }
  attr(out, "draws") <- rbind(data.frame(set = "null", null, check.names = FALSE),
                              data.frame(set = "fitted", fitted, check.names = FALSE))
  out
}

#' Plot the profile of a fit of population shifts
#'
#' Draws, along the width `w` of the competition kernel (on a log scale), the \eqn{R^2} of the
#' fit, the estimate of displacement `c` and the slopes `h` of [population_shifts()], with the
#' chosen width as a dashed line and the lower end of the default grid (the pooled standard
#' deviation within populations, when it sets it) as a dotted line. A flat \eqn{R^2} means that
#' `w` is not identified, and estimates that change along it are then not to be trusted. With a
#' bootstrap, the interval of `c` is drawn at the chosen width.
#'
#' @param x Result of [population_shifts()] with the width estimated (more than one width in
#'   its profile).
#' @param ... Not used.
#' @return A ggplot object.
#' @export
plot.population_shifts <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot() for population_shifts() needs the package ggplot2", call. = FALSE) # nocov
  }
  profile <- attr(x, "profile")
  if (is.null(profile) || nrow(profile) < 2) {
    stop("the width was fixed, so there is no profile to draw", call. = FALSE)
  }
  profile <- profile[order(profile$width), , drop = FALSE]
  value_cols <- setdiff(names(profile), c("width", "rss"))
  long <- do.call(rbind, lapply(value_cols, function(v) {
    data.frame(width = profile$width, quantity = v, value = profile[[v]])
  }))
  long$quantity <- factor(long$quantity, value_cols)
  .data <- ggplot2::.data
  p <- ggplot2::ggplot(long, ggplot2::aes(.data$width, .data$value)) +
    ggplot2::geom_hline(data = data.frame(quantity = factor(setdiff(value_cols, "r_squared"),
                                                            value_cols)),
                        ggplot2::aes(yintercept = 0), colour = "grey70") +
    ggplot2::geom_line(colour = "#0072B2", na.rm = TRUE) +
    ggplot2::geom_point(colour = "#0072B2", size = 1.2, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = x$width, linetype = "dashed", colour = "#D55E00")
  floor <- attr(x, "floor")
  if (!is.null(floor)) {
    p <- p + ggplot2::geom_vline(xintercept = floor, linetype = "dotted", colour = "grey40")
  }
  if (!is.null(x$displacement_low)) {
    band <- data.frame(width = x$width, quantity = factor("displacement", value_cols),
                       low = x$displacement_low, high = x$displacement_high)
    p <- p + ggplot2::geom_errorbar(data = band, ggplot2::aes(x = .data$width,
                                                              ymin = .data$low,
                                                              ymax = .data$high),
                                    inherit.aes = FALSE, width = 0.05, colour = "#D55E00")
  }
  p +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(ggplot2::vars(.data$quantity), scales = "free_y", ncol = 1) +
    ggplot2::labs(x = "Width of the competition kernel (w)", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
