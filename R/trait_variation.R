# Analyses of trait variation among individuals, populations and communities, written from the
# published definitions as replacements for the functions of the package cati (Taudiere &
# Violle 2016), several of which give wrong results in version 0.99.6.

#' Summary of the null-model tests of T-statistics
#'
#' Summarises the result of [tstats()] with null models: for each trait and statistic, the share
#' of sites below and above the null distribution, the mean standardised effect size and a test
#' of the mean over sites, and the sites that deviate.
#'
#' Each site is tested in both directions from the null values kept by [tstats()]: the p-value
#' of a low value is `(1 + r) / (1 + n)`, where `r` is the number of null values at or below the
#' observed one and `n` the number of null values, and that of a high value counts null values
#' at or above it. A site deviates when one of the two is below `level / 2`, so that the two
#' tests together have size `level`. A site cannot be tested when its observed value is missing
#' or when every null value equals the observed one, as for `T_IP.IC` at a site with one
#' species, where the statistic is 1 whatever the permutation; such sites get `NA` and are left
#' out of the shares, instead of a small p-value in both directions.
#'
#' The test of the mean compares the mean of the observed statistic over the testable sites with
#' the mean over the same sites of each randomisation, so that the null distribution of the
#' mean keeps the dependence among sites within a randomisation.
#'
#' @param x Data frame from [tstats()] with `nperm > 0`.
#' @param level Size of the two-sided test of each site and of the mean. A site can deviate
#'   only if `1 / (1 + nperm) < level / 2`, so `nperm` must be at least `2 / level`; a warning
#'   is given otherwise.
#' @return A list with two data frames. `traits` has one row per trait and statistic, with the
#'   columns `trait`, `statistic`, `n_sites` (testable sites), `percent_low`, `percent_high`,
#'   `mean_ses`, `p_mean_low`, `p_mean_high`, `sites_low` and `sites_high` (site names separated
#'   by spaces). `sites` has one row per site, trait and statistic, with `observed`, `ses`,
#'   `p_low`, `p_high` and `deviation` (`"low"`, `"high"`, `"none"` or `NA`).
#' @seealso [tstats()].
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(20)
#' sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
#'                               competition_matrix(tree, "none"), n_steps = 100,
#'                               seed_steps = 1:5)
#' res <- tstats(sample_individuals(sim, intraspecific_sd = 2), nperm = 99)
#' tstats_summary(res)$traits
#' @export
tstats_summary <- function(x, level = 0.05) {
  null <- attr(x, "null")
  if (!is.data.frame(x) || is.null(null)) {
    stop("`x` must come from tstats() with nperm > 0", call. = FALSE)
  }
  stopifnot(is.numeric(level), length(level) == 1, level > 0, level < 1)
  ratios <- c("T_IP.IC", "T_IC.IR", "T_PC.PR")
  if (1 / (1 + dim(null)[1]) >= level / 2) {
    warning("with ", dim(null)[1], " randomisations no p-value can fall below level / 2 = ",
            level / 2, "; use nperm >= ", ceiling(2 / level), call. = FALSE)
  }
  sites <- list()
  traits <- list()
  for (stat in ratios) {
    sims <- matrix(null[, , stat], dim(null)[1])
    sims[!is.finite(sims)] <- NA
    obs <- x[[stat]]
    tested <- site_tests(obs, sims, level)
    testable <- tested$testable
    p_low <- tested$p_low
    p_high <- tested$p_high
    deviation <- tested$deviation
    ses <- tested$ses
    sites[[stat]] <- data.frame(site = x$site, trait = x$trait, statistic = stat, observed = obs,
                                ses = ses, p_low = p_low, p_high = p_high,
                                deviation = deviation, stringsAsFactors = FALSE)
    for (tr in unique(x$trait)) {
      rows <- which(x$trait == tr & testable)
      if (length(rows) == 0) {
        mean_low <- mean_high <- NA_real_
      } else {
        null_mean <- rowMeans(sims[, rows, drop = FALSE])
        obs_mean <- mean(obs[rows])
        ok <- !is.na(null_mean)
        mean_low <- (1 + sum(null_mean[ok] <= obs_mean)) / (1 + sum(ok))
        mean_high <- (1 + sum(null_mean[ok] >= obs_mean)) / (1 + sum(ok))
      }
      dev <- deviation[rows]
      traits[[paste(stat, tr)]] <- data.frame(
        trait = tr, statistic = stat, n_sites = length(rows),
        percent_low = if (length(rows)) 100 * mean(dev == "low") else NA_real_,
        percent_high = if (length(rows)) 100 * mean(dev == "high") else NA_real_,
        mean_ses = if (length(rows)) mean(ses[rows]) else NA_real_,
        p_mean_low = mean_low, p_mean_high = mean_high,
        sites_low = paste(x$site[rows][dev == "low"], collapse = " "),
        sites_high = paste(x$site[rows][dev == "high"], collapse = " "),
        stringsAsFactors = FALSE)
    }
  }
  traits <- do.call(rbind, traits)
  traits <- traits[order(match(traits$trait, unique(x$trait))), ]
  rownames(traits) <- NULL
  sites <- do.call(rbind, sites)
  rownames(sites) <- NULL
  list(traits = traits, sites = sites)
}

#' Intraspecific variability and species turnover in community trait means
#'
#' Splits the variation among sites of the community-weighted mean of a trait into the parts due
#' to species turnover, to intraspecific variability and to their covariation, following Lepš et
#' al. (2011), from the trait values of individuals.
#'
#' At each site, the "specific" mean weights the mean trait of each population (a species at
#' the site) by its abundance, and the "fixed" mean weights the mean of the species over all
#' individuals in the data by the same abundances. Their difference is the part of the
#' community mean due to intraspecific variability. The sums of squares among sites of the
#' specific mean (total), of the fixed mean (turnover) and of the difference (intraspecific)
#' are related by `SS_total = SS_turnover + SS_intraspecific + SS_covariation`; a positive
#' covariation means that turnover and intraspecific variability move the community mean the
#' same way. With a `formula`, each of the three site means is regressed on the variables of
#' `env`, and the sums of squares are split into the parts explained by the model and the
#' residuals, with the p-value of each regression.
#'
#' Abundances are those of the `abundance` column when given (one value per population, as in
#' [sample_individuals()]), otherwise the number of individuals measured. Individuals with a
#' missing trait value are left out of the population means, the species means and the counts
#' of that trait together, so that both community means use the same populations and weights.
#'
#' @param individuals Data frame with one row per individual.
#' @param traits Names of the columns with trait values.
#' @param site,species Names of the columns with the site and the species.
#' @param abundance Name of a column with the abundance of each population (the same value
#'   for all its individuals), or `NULL` to weight by the number of individuals measured.
#' @param formula One-sided formula with variables of `env`, with an intercept; `~ 1` (the
#'   default) gives the total sums of squares only. With several terms, `explained` is the sum
#'   of squares of the whole model and `p` the p-value of its F test.
#' @param env Data frame with one row per site, a column named like `site` and the variables of
#'   `formula`; it is matched to the sites by name. Sites with missing values in these
#'   variables are left out of all sums of squares, with a warning.
#' @return A list with `decomposition`, a data frame with one row per trait and component
#'   (`turnover`, `intraspecific`, `covariation`, `total`) and the columns `explained`,
#'   `residual`, `total` (sums of squares) and `p` (of the regression, `NA` for the
#'   covariation and without a model), and `sites`, a data frame with the specific and fixed
#'   means and their difference for each site and trait. The list has class `trait_flex`,
#'   with a [plot()][plot.trait_flex] method.
#' @references Lepš, J., de Bello, F., Šmilauer, P., & Doležal, J. (2011). Community trait
#'   response to environment: disentangling species turnover vs intraspecific trait
#'   variability effects. *Ecography*, 34, 856--863. \doi{10.1111/j.1600-0587.2010.06904.x}
#' @seealso [tstats()], [nested_variance()].
#' @examples
#' set.seed(1)
#' ind <- data.frame(site = rep(paste0("s", 1:6), each = 20),
#'                   species = sample(c("a", "b", "c"), 120, replace = TRUE))
#' env <- data.frame(site = paste0("s", 1:6), moisture = 1:6)
#' ind$trait <- c(a = 10, b = 20, c = 30)[ind$species] +
#'   env$moisture[match(ind$site, env$site)] + rnorm(120)
#' res <- trait_flex(ind, formula = ~ moisture, env = env)
#' res$decomposition
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(res)
#' @export
trait_flex <- function(individuals, traits = "trait", site = "site", species = "species",
                       abundance = NULL, formula = ~1, env = NULL) {
  stopifnot(is.data.frame(individuals), is.character(traits), length(traits) >= 1)
  needed <- c(traits, site, species, abundance)
  missing_columns <- setdiff(needed, names(individuals))
  if (length(missing_columns) > 0) {
    stop("`individuals` has no column ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  stopifnot(inherits(formula, "formula"), length(formula) == 2)
  if (attr(stats::terms(formula), "intercept") == 0) {
    stop("`formula` must keep the intercept: the sums of squares are taken around the mean",
         call. = FALSE)
  }
  site_id <- as.character(individuals[[site]])
  species_id <- as.character(individuals[[species]])
  site_levels <- unique(site_id)
  model_vars <- all.vars(formula)
  if (length(model_vars) > 0) {
    if (!is.data.frame(env) || !site %in% names(env)) {
      stop("`env` must be a data frame with a column `", site, "`", call. = FALSE)
    }
    missing_vars <- setdiff(model_vars, names(env))
    if (length(missing_vars) > 0) {
      stop("`env` has no column ", paste(missing_vars, collapse = ", "), call. = FALSE)
    }
    rows <- match(site_levels, as.character(env[[site]]))
    if (anyNA(rows)) stop("`env` lacks some sites of `individuals`", call. = FALSE)
    env_sites <- env[rows, model_vars, drop = FALSE]
    env_ok <- stats::complete.cases(env_sites)
    if (!all(env_ok)) {
      warning("sites with missing values in `env` are left out of the decomposition: ",
              paste(site_levels[!env_ok], collapse = ", "), call. = FALSE)
    }
    model_formula <- stats::reformulate(attr(stats::terms(formula), "term.labels"), ".y")
  } else {
    env_sites <- data.frame(row.names = seq_along(site_levels))
    env_ok <- rep(TRUE, length(site_levels))
  }
  if (!is.null(abundance)) {
    a <- individuals[[abundance]]
    pops <- paste(site_id, species_id, sep = "\r")
    if (anyNA(a)) stop("`abundance` has missing values", call. = FALSE)
    if (any(tapply(a, pops, function(v) length(unique(v))) > 1)) {
      stop("`abundance` must be the same for all individuals of a population", call. = FALSE)
    }
  }
  decomposition <- list()
  site_means <- list()
  for (trait in traits) {
    z <- individuals[[trait]]
    if (!is.numeric(z)) stop("trait `", trait, "` must be numeric", call. = FALSE)
    keep <- !is.na(z)
    pop <- paste(site_id, species_id, sep = "\r")[keep]
    first <- !duplicated(pop)
    pop_site <- site_id[keep][first]
    pop_species <- species_id[keep][first]
    pop_mean <- tapply(z[keep], factor(pop, unique(pop)), mean)
    weight <- if (is.null(abundance)) {
      as.vector(table(factor(pop, unique(pop))))
    } else {
      as.numeric(individuals[[abundance]][keep][first])
    }
    species_mean <- tapply(z[keep], species_id[keep], mean)[pop_species]
    weighted <- function(v) {
      sapply(site_levels, function(s) {
        i <- pop_site == s
        if (!any(i)) NA_real_ else sum(weight[i] * v[i]) / sum(weight[i])
      })
    }
    specific <- weighted(as.vector(pop_mean))
    fixed <- weighted(as.vector(species_mean))
    values <- list(turnover = fixed, intraspecific = specific - fixed, total = specific)
    ss <- lapply(values, function(y) {
      ok <- !is.na(y) & env_ok
      total <- sum((y[ok] - mean(y[ok]))^2)
      if (length(model_vars) == 0) return(c(explained = 0, residual = total, total = total,
                                             p = NA_real_))
      fit <- stats::lm(model_formula, data = cbind(.y = y, env_sites)[ok, , drop = FALSE])
      residual <- sum(stats::residuals(fit)^2)
      f <- summary(fit)$fstatistic
      p <- if (is.null(f)) NA_real_ else stats::pf(f[1], f[2], f[3], lower.tail = FALSE)
      c(explained = total - residual, residual = residual, total = total, p = unname(p))
    })
    ss$covariation <- ss$total - ss$turnover - ss$intraspecific
    ss$covariation["p"] <- NA_real_
    ss <- ss[c("turnover", "intraspecific", "covariation", "total")]
    decomposition[[trait]] <- data.frame(trait = trait, component = names(ss),
                                         do.call(rbind, ss), row.names = NULL)
    site_means[[trait]] <- data.frame(site = site_levels, trait = trait, specific = specific,
                                      fixed = fixed, intraspecific = specific - fixed,
                                      row.names = NULL, stringsAsFactors = FALSE)
  }
  structure(list(decomposition = do.call(rbind, unname(decomposition)),
                 sites = do.call(rbind, unname(site_means))), class = "trait_flex")
}

#' Plot the decomposition of community trait means
#'
#' Draws, for each trait, the sums of squares among sites of the community mean due to species
#' turnover, intraspecific variability, their covariation and the total, from [trait_flex()].
#' With an environmental model each bar is split into the part explained by the model and the
#' residual part. A negative covariation (turnover and intraspecific variability moving the
#' community mean in opposite directions) is drawn below zero. Needs the package ggplot2.
#'
#' @param x Result of [trait_flex()].
#' @param relative If `TRUE` (the default), sums of squares are divided by the total sum of
#'   squares of each trait, so that traits on different scales can be compared.
#' @param colours Two colours, for the explained and the residual parts (only the first is
#'   used without a model). The default is the blue of the Okabe-Ito palette and a light grey,
#'   as in [plot.nested_variance()].
#' @param ... Not used.
#' @return A ggplot object, which can be changed further with ggplot2 functions.
#' @seealso [trait_flex()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' set.seed(1)
#' ind <- data.frame(site = rep(paste0("s", 1:8), each = 20),
#'                   species = sample(c("a", "b", "c"), 160, replace = TRUE))
#' env <- data.frame(site = paste0("s", 1:8), moisture = 1:8)
#' shift <- env$moisture[match(ind$site, env$site)]
#' ind$height <- c(a = 10, b = 20, c = 30)[ind$species] + shift + rnorm(160)
#' ind$leaf <- c(a = 3, b = 2, c = 1)[ind$species] - 0.3 * shift + rnorm(160)
#' plot(trait_flex(ind, traits = c("height", "leaf"), formula = ~ moisture, env = env))
#' @export
plot.trait_flex <- function(x, relative = TRUE, colours = c("#0072B2", "grey80"), ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot() of trait_flex() needs the package ggplot2", call. = FALSE) # nocov
  }
  stopifnot(is.logical(relative), length(relative) == 1, !is.na(relative))
  check_colours(colours, 2)
  dec <- x$decomposition
  components <- c("turnover", "intraspecific", "covariation", "total")
  scale <- if (relative) {
    dec$total[dec$component == "total"][match(dec$trait, dec$trait[dec$component == "total"])]
  } else {
    1
  }
  has_model <- any(dec$explained != 0, na.rm = TRUE)
  parts <- if (has_model) c(explained = "Explained by the model", residual = "Residual") else
    c(total = "Sum of squares")
  long <- do.call(rbind, lapply(names(parts), function(part) {
    data.frame(trait = factor(dec$trait, unique(dec$trait)),
               component = factor(dec$component, components),
               part = factor(parts[[part]], rev(parts)), value = dec[[part]] / scale)
  }))
  fill <- stats::setNames(colours[seq_along(parts)], unname(parts))
  .data <- ggplot2::.data
  ggplot2::ggplot(long, ggplot2::aes(x = .data$component, y = .data$value, fill = .data$part,
                                     group = .data$part)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_wrap(ggplot2::vars(.data$trait), scales = if (relative) "fixed" else "free_y") +
    ggplot2::scale_fill_manual(values = fill, breaks = unname(parts), name = NULL,
                               guide = if (has_model) "legend" else "none") +
    ggplot2::scale_y_continuous(labels = if (relative) function(v) paste0(100 * v, "%") else
      ggplot2::waiver()) +
    ggplot2::labs(x = NULL, y = if (relative) "Share of the total sum of squares" else
      "Sum of squares") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), legend.position = "top")
}

#' Variance components of traits across nested levels
#'
#' Splits the variance of each trait among nested levels, such as regions, sites, species and
#' populations, with a linear mixed model of random intercepts nested from the broadest level
#' to the finest, fitted by restricted maximum likelihood (Messier et al. 2010).
#'
#' The levels must be nested: each level of a factor must lie within a single level of the
#' factor before it. Two checks guard against the common mistakes. A factor that does not split
#' any group of the factors before it (for instance species after populations, since each
#' population holds a single species) adds no level of variation, and the function stops,
#' because the order of the factors is then most likely reversed. A factor whose levels occur
#' under several groups of the factors before it (for instance genera, found at several sites,
#' after sites) is crossed with them rather than nested; the model then treats each combination
#' as a separate group, which is rarely what is meant, and the function warns.
#'
#' @param data Data frame with the traits and the factors.
#' @param traits Names of the columns with trait values.
#' @param factors Names of the columns with the nested factors, from the broadest level to the
#'   finest.
#' @return A data frame with one row per factor and one for the variance within the finest
#'   level (`within`), and one column per trait with the share of variance of each level. The
#'   attribute `variance` holds the variances. It has class `nested_variance`, with a
#'   [plot()][plot.nested_variance] method.
#' @references Messier, J., McGill, B. J., & Lechowicz, M. J. (2010). How do traits vary across
#'   ecological scales? A case for trait-based ecology. *Ecology Letters*, 13, 838--848.
#'   \doi{10.1111/j.1461-0248.2010.01476.x}
#' @seealso [trait_flex()], [tstats()].
#' @examplesIf requireNamespace("nlme", quietly = TRUE)
#' set.seed(1)
#' d <- data.frame(region = rep(c("north", "south"), each = 60),
#'                 site = rep(paste0("s", 1:6), each = 20),
#'                 species = sample(c("a", "b", "c"), 120, replace = TRUE))
#' d$population <- paste(d$site, d$species)
#' d$trait <- c(a = 10, b = 20, c = 30)[d$species] + rnorm(120)
#' res <- nested_variance(d, "trait", c("region", "site", "population"))
#' res
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(res)
#' @export
nested_variance <- function(data, traits, factors) {
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("nested_variance() needs the package nlme", call. = FALSE) # nocov
  }
  stopifnot(is.data.frame(data), is.character(traits), is.character(factors),
            length(factors) >= 1)
  missing_columns <- setdiff(c(traits, factors), names(data))
  if (length(missing_columns) > 0) {
    stop("`data` has no column ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  single <- factors[vapply(factors, function(f) length(unique(stats::na.omit(data[[f]]))) < 2,
                           logical(1))]
  if (length(single) > 0) {
    stop("factor ", paste0("`", single, "`", collapse = ", "), " has a single level",
         call. = FALSE)
  }
  groups <- lapply(seq_along(factors), function(k) {
    do.call(paste, c(lapply(factors[seq_len(k)], function(f) as.character(data[[f]])),
                     sep = "\r"))
  })
  for (k in seq_along(factors)[-1]) {
    if (length(unique(groups[[k]])) == length(unique(groups[[k - 1]]))) {
      stop("`", factors[k], "` does not split any level of `", factors[k - 1],
           "`; order `factors` from the broadest level to the finest", call. = FALSE)
    }
    parents <- tapply(groups[[k - 1]], as.character(data[[factors[k]]]),
                      function(v) length(unique(v)))
    if (any(parents > 1)) {
      warning("labels of `", factors[k], "` recur under several levels of `",
              paste(factors[seq_len(k - 1)], collapse = "/"), "`. If they are the same units, ",
              "the factors are crossed and the nested model treats each combination as a ",
              "separate group; if the labels are only reused (site 1, 2, ... in each region), ",
              "the results are right and this warning can be ignored", call. = FALSE)
    }
  }
  model_data <- data.frame(lapply(data[factors], function(f) factor(as.character(f))))
  random <- stats::as.formula(paste("~ 1 |", paste(factors, collapse = "/")))
  variances <- sapply(traits, function(trait) {
    model_data$.y <- data[[trait]]
    fit <- nlme::lme(.y ~ 1, random = random, data = model_data, na.action = stats::na.omit,
                     method = "REML")
    v <- ape::varcomp(fit, scale = FALSE)
    stats::setNames(as.vector(v), c(factors, "within"))
  })
  variances <- matrix(variances, ncol = length(traits),
                      dimnames = list(c(factors, "within"), traits))
  out <- as.data.frame(sweep(variances, 2, colSums(variances), "/"))
  out <- data.frame(level = rownames(out), out, row.names = NULL, check.names = FALSE)
  attr(out, "variance") <- variances
  class(out) <- c("nested_variance", class(out))
  out
}

#' Plot the shares of variance across nested levels
#'
#' Draws one stacked bar per trait with the share of variance of each level from
#' [nested_variance()], the broadest level at the bottom and the variance within the finest
#' level (`within`, in grey) at the top. Needs the package ggplot2.
#'
#' The levels take colours of the Okabe-Ito palette (Okabe & Ito 2008), which stay distinct for
#' the common forms of colour blindness, in an order that alternates darker and lighter colours
#' so that neighbouring segments also differ in grey scale: blue, orange, bluish green, reddish
#' purple, vermillion and sky blue. Yellow is left out because it hardly differs from the grey
#' of `within` and from a white background. Percentages are written in white on dark colours
#' and in black on light ones.
#'
#' @param x Result of [nested_variance()].
#' @param min_label Smallest share labelled with its percentage inside the bar; 1 or more
#'   removes the labels.
#' @param colours Colours of the levels, one per row of `x` (the last for `within`), or `NULL`
#'   (the default) for the Okabe-Ito colours above and a light grey for `within`. With more
#'   than six factors, `colours` must be given.
#' @param ... Not used.
#' @return A ggplot object, which can be changed further with ggplot2 functions.
#' @references Okabe, M., & Ito, K. (2008). *Color Universal Design (CUD): Colorblind Barrier
#'   Free*. \url{https://jfly.uni-koeln.de/color/}
#' @seealso [nested_variance()], [plot.trait_flex()].
#' @examplesIf requireNamespace("nlme", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
#' set.seed(1)
#' d <- data.frame(region = rep(c("north", "south"), each = 60),
#'                 site = rep(paste0("s", 1:6), each = 20),
#'                 species = sample(c("a", "b", "c"), 120, replace = TRUE))
#' d$population <- paste(d$site, d$species)
#' d$size <- c(a = 10, b = 20, c = 30)[d$species] + rnorm(120)
#' d$shape <- rnorm(120, rep(c(0, 2), each = 60))
#' plot(nested_variance(d, c("size", "shape"), c("region", "site", "population")))
#' @export
plot.nested_variance <- function(x, min_label = 0.05, colours = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot() of nested_variance() needs the package ggplot2", call. = FALSE) # nocov
  }
  stopifnot(is.numeric(min_label), length(min_label) == 1)
  levels <- x$level
  traits <- setdiff(names(x), "level")
  long <- data.frame(level = factor(rep(levels, length(traits)), levels),
                     trait = factor(rep(traits, each = length(levels)), traits),
                     share = unlist(lapply(traits, function(t) x[[t]]), use.names = FALSE))
  long$label <- ifelse(long$share >= min_label, sprintf("%.0f%%", 100 * long$share), "")
  n_factors <- length(levels) - 1
  if (is.null(colours)) {
    okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9")
    if (n_factors > length(okabe_ito)) {
      stop("more than ", length(okabe_ito), " factors: give one colour per level in `colours`",
           call. = FALSE)
    }
    colours <- c(okabe_ito[seq_len(n_factors)], "grey80")
  }
  check_colours(colours, length(levels))
  fill <- stats::setNames(colours, levels)
  dark <- colSums(grDevices::col2rgb(fill) * c(0.299, 0.587, 0.114)) < 140
  long$text <- ifelse(dark[as.character(long$level)], "white", "black")
  .data <- ggplot2::.data
  ggplot2::ggplot(long, ggplot2::aes(x = .data$trait, y = .data$share, fill = .data$level,
                                     group = .data$level)) +
    ggplot2::geom_col(position = ggplot2::position_stack(reverse = TRUE), width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label, colour = .data$text),
                       position = ggplot2::position_stack(vjust = 0.5, reverse = TRUE),
                       size = 3.5, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = fill, name = "Level") +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_y_continuous(labels = function(v) paste0(100 * v, "%"),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = NULL, y = "Share of variance") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
}

#' Trait indices of communities against null models
#'
#' Computes any index of the trait values of the individuals of each site (community), such as
#' the variance, the range or a nearest-neighbour distance, and compares it with randomised
#' communities.
#'
#' The null models follow Taudiere & Violle (2016), written from their description:
#'
#' * `"local"`: trait values are permuted among the individuals of each site; the values of the
#'   site are kept and only their link with species is broken, so it is useful for indices that
#'   use `species`;
#' * `"regional.ind"`: each site receives as many trait values as it has individuals, drawn
#'   without replacement from all individuals of the data;
#' * `"regional.pop"`: each individual takes the mean of its population, and each site receives
#'   as many of these values as it has individuals, drawn without replacement from all
#'   individuals. The observed index uses the values of the individuals, so the null removes
#'   the variation within populations as well as the species composition of the site;
#' * `"regional.pop.prab"`: the index is computed on the means of the populations of each site
#'   (one value per population, observed and null alike), and each site receives as many
#'   population means as it has populations, drawn without replacement from all populations.
#'
#' In all cases the individuals (or populations) keep their species, so an index with an
#' argument named `species`, such as `function(x, species)`, receives the species of the
#' values of the site. With `multivariate = TRUE` the index
#' receives the matrix of all traits of the site, and randomisations move whole rows, so that
#' the traits of an individual stay together. Sites and species are matched by name, so the
#' order of the rows does not matter. An index that the null model cannot change (for example
#' `stats::var` under `"local"`) gets `ses = NA`.
#'
#' @param individuals Data frame with one row per individual.
#' @param index Function of the trait values of a site (its first argument: a vector, or a
#'   matrix with `multivariate = TRUE`) and optionally of their species (an argument named
#'   `species`), returning a single number, such as `stats::var`, [mnnd()] or
#'   `function(x, species) stats::var(tapply(x, species, mean))`; or a named list of such
#'   functions, all computed on the same observed and randomised communities, such as
#'   `list(variance = stats::var, MNND = mnnd, range = trait_range)`.
#' @param traits Names of the columns with trait values.
#' @param site,species Names of the columns with the site and the species.
#' @param null Null model (see Details).
#' @param nperm Number of randomisations; 0 for the observed values only.
#' @param alternative Direction of the test: `"two.sided"` (the default), `"less"` or
#'   `"greater"`.
#' @param multivariate If `TRUE`, the index receives all traits at once.
#' @return A data frame of class `community_index` with one row per site, trait (or
#'   `"multivariate"` with `multivariate = TRUE`) and index, and the columns `site`, `trait`,
#'   `index` (the name in the list, or the name of the function), `observed`, and with
#'   `nperm > 0` also `ses` and `p`; the attribute `null` holds the null values
#'   (randomisation x row), and [plot()][plot.community_index] draws the tests. Individuals with a missing value are left out for that trait
#'   (for any trait with `multivariate = TRUE`).
#' @references Taudiere, A., & Violle, C. (2016). cati: an R package using functional traits
#'   to detect and quantify multi-level community assembly processes. *Ecography*, 39,
#'   699--708. \doi{10.1111/ecog.01433}
#' @seealso [tstats()].
#' @examples
#' set.seed(1)
#' ind <- data.frame(site = rep(paste0("s", 1:5), each = 20),
#'                   species = sample(letters[1:4], 100, replace = TRUE))
#' ind$trait <- c(a = 1, b = 2, c = 3, d = 4)[ind$species] + rnorm(100, sd = 0.3)
#' # range of trait values within sites, against draws from the whole region
#' community_index(ind, function(x) diff(range(x)), null = "regional.ind", nperm = 99)
#' # variance among the population means of a site
#' pop_var <- function(x, species) stats::var(tapply(x, species, mean))
#' community_index(ind, pop_var, null = "local", nperm = 99)
#' # several indices at once, on the same randomised communities
#' community_index(ind, list(variance = stats::var, MNND = mnnd, kurtosis = trait_kurtosis),
#'                 null = "regional.ind", nperm = 99)
#' @export
community_index <- function(individuals, index, traits = "trait", site = "site",
                            species = "species",
                            null = c("local", "regional.ind", "regional.pop",
                                     "regional.pop.prab"),
                            nperm = 99, alternative = c("two.sided", "less", "greater"),
                            multivariate = FALSE) {
  null <- match.arg(null)
  alternative <- match.arg(alternative)
  if (is.function(index)) {
    label <- paste(deparse(substitute(index)), collapse = "")
    if (nchar(label) > 30 || grepl("^function", label)) label <- "index"
    index <- stats::setNames(list(index), label)
  }
  if (!is.list(index) || length(index) == 0 || !all(vapply(index, is.function, logical(1))) ||
      is.null(names(index)) || any(names(index) == "") || anyDuplicated(names(index))) {
    stop("`index` must be a function or a list of functions with distinct names",
         call. = FALSE)
  }
  stopifnot(is.data.frame(individuals), is.character(traits),
            length(traits) >= 1, isTRUE(multivariate) || isFALSE(multivariate))
  missing_columns <- setdiff(c(traits, site, species), names(individuals))
  if (length(missing_columns) > 0) {
    stop("`individuals` has no column ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  nperm <- if (is.null(nperm)) 0 else nperm
  stopifnot(length(nperm) == 1, nperm >= 0, nperm == round(nperm))
  index_names <- names(index)
  n_index <- length(index)
  uses_species <- vapply(index, function(f) "species" %in% names(formals(f)), logical(1))
  apply_index <- function(x, sp) {
    vapply(seq_len(n_index), function(k) {
      value <- if (uses_species[k]) index[[k]](x, species = sp) else index[[k]](x)
      if (!is.numeric(value) || length(value) != 1) {
        stop("index `", index_names[k], "` must return a single number", call. = FALSE)
      }
      as.numeric(value)
    }, numeric(1))
  }
  site_all <- as.character(individuals[[site]])
  species_all <- as.character(individuals[[species]])
  site_levels <- unique(site_all)
  sets <- if (multivariate) list(multivariate = traits) else stats::setNames(as.list(traits), traits)
  per_set <- lapply(names(sets), function(label) {
    values <- as.matrix(individuals[sets[[label]]])
    if (!is.numeric(values)) stop("the columns in `traits` must be numeric", call. = FALSE)
    keep <- stats::complete.cases(values)
    values <- values[keep, , drop = FALSE]
    site_id <- site_all[keep]
    species_id <- species_all[keep]
    if (null == "regional.pop" || null == "regional.pop.prab") {
      pop <- paste(site_id, species_id, sep = "\r")
      pop_means <- rowsum(values, pop, reorder = FALSE) / as.vector(table(pop)[unique(pop)])
    }
    if (null == "regional.pop.prab") {
      # one row per population: the index sees population means
      first <- !duplicated(pop)
      values <- pop_means[pop[first], , drop = FALSE]
      site_id <- site_id[first]
      species_id <- species_id[first]
    } else if (null == "regional.pop") {
      carried <- pop_means[pop, , drop = FALSE]
    }
    rows_of <- split(seq_along(site_id), factor(site_id, site_levels))
    # one value per site and index, sites varying fastest
    site_index <- function(v) {
      as.vector(t(matrix(vapply(rows_of, function(i) {
        if (length(i) == 0) return(rep(NA_real_, n_index))
        x <- v[i, , drop = FALSE]
        apply_index(if (multivariate) x else x[, 1], species_id[i])
      }, numeric(n_index)), n_index)))
    }
    observed <- site_index(values)
    sims <- matrix(NA_real_, nperm, length(site_levels) * n_index)
    n <- nrow(values)
    for (r in seq_len(nperm)) {
      shuffled <- values
      if (null == "local") {
        for (i in rows_of) shuffled[i, ] <- values[i[sample.int(length(i))], , drop = FALSE]
      } else {
        pool <- if (null == "regional.pop") carried else values
        for (i in rows_of) {
          shuffled[i, ] <- pool[sample.int(n, length(i)), , drop = FALSE]
        }
      }
      sims[r, ] <- site_index(shuffled)
    }
    list(table = data.frame(site = rep(site_levels, n_index), trait = label,
                            index = rep(index_names, each = length(site_levels)),
                            observed = observed, row.names = NULL, stringsAsFactors = FALSE),
         null = sims)
  })
  out <- do.call(rbind, lapply(per_set, `[[`, "table"))
  rownames(out) <- NULL
  class(out) <- c("community_index", class(out))
  if (nperm == 0) return(out)
  sims <- do.call(cbind, lapply(per_set, `[[`, "null"))
  compared <- null_compare(out$observed, sims, alternative)
  out$ses <- compared$ses
  out$p <- compared$p
  attr(out, "null") <- sims
  out
}

# stops unless `colours` holds n valid colours
check_colours <- function(colours, n) {
  ok <- is.character(colours) && length(colours) == n &&
    !inherits(try(grDevices::col2rgb(colours), silent = TRUE), "try-error")
  if (!ok) stop("`colours` must be ", n, " valid colours", call. = FALSE)
  invisible(colours)
}

# one-sided tests of each site against its null values (columns of sims), as in
# tstats_summary(): sites with a missing observed value, no valid null value or every null value
# tied with the observed one cannot be tested and get NA
site_tests <- function(obs, sims, level) {
  sims[!is.finite(sims)] <- NA
  nperm <- nrow(sims)
  valid <- colSums(!is.na(sims))
  ties <- colSums(sims == rep(obs, each = nperm), na.rm = TRUE)
  testable <- !is.na(obs) & valid > 0 & ties < valid
  below <- colSums(sims <= rep(obs, each = nperm), na.rm = TRUE)
  above <- colSums(sims >= rep(obs, each = nperm), na.rm = TRUE)
  p_low <- ifelse(testable, (1 + below) / (1 + valid), NA_real_)
  p_high <- ifelse(testable, (1 + above) / (1 + valid), NA_real_)
  deviation <- ifelse(p_low < level / 2, "low", ifelse(p_high < level / 2, "high", "none"))
  ses <- null_compare(obs, sims)$ses
  ses[!testable] <- NA_real_
  list(testable = testable, p_low = p_low, p_high = p_high, deviation = deviation, ses = ses)
}
