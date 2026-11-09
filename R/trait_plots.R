# Plots of null-model tests and of trait distributions of populations, with ggplot2 (Suggests).

#' Null intervals of the tests of each site
#'
#' Summarises, for each site, trait and statistic, the null distribution kept by [tstats()] or
#' [community_index()]: its median, an outer interval (`level`) and an inner one (`inner`),
#' with the observed value, the standardised effect size and the one-sided tests of
#' [tstats_summary()]. This is the table drawn by the `plot()` methods of these objects, and
#' can be used to build other figures.
#'
#' @param x Result of [tstats()] or [community_index()] with `n_perm > 0`.
#' @param level Coverage of the outer null interval, and `1 - level` the size of the two-sided
#'   test that defines `deviation`.
#' @param inner Coverage of the inner null interval.
#' @return A data frame with one row per site, trait and statistic (for [community_index()],
#'   `statistic` is the name of the index) and the columns `observed`, `null_median`, `null_low`,
#'   `null_high`, `inner_low`, `inner_high`, `ses`, `p_low`, `p_high` and `deviation`
#'   (`"low"`, `"high"`, `"none"` or `NA` for sites that cannot be tested).
#' @seealso [plot.tstats()], [tstats_summary()].
#' @family trait variation among individuals
#' @examples
#' set.seed(1)
#' ind <- data.frame(site = rep(paste0("s", 1:6), each = 20),
#'                   species = sample(c("a", "b", "c"), 120, replace = TRUE))
#' ind$trait <- c(a = 10, b = 20, c = 30)[ind$species] + rnorm(120)
#' head(null_intervals(tstats(ind, n_perm = 99)))
#' @export
null_intervals <- function(x, level = 0.95, inner = 0.5) {
  null <- attr(x, "null")
  if (!is.data.frame(x) || is.null(null)) {
    stop("`x` must come from tstats() or community_index() with n_perm > 0", call. = FALSE)
  }
  stopifnot(is.numeric(level), length(level) == 1, level > 0, level < 1,
            is.numeric(inner), length(inner) == 1, inner > 0, inner < level)
  parts <- if (length(dim(null)) == 3) {
    lapply(dimnames(null)[[3]], function(stat) {
      list(statistic = stat, obs = x[[stat]], sims = matrix(null[, , stat], dim(null)[1]))
    })
  } else {
    list(list(statistic = x$index, obs = x$observed, sims = matrix(null, nrow(null))))
  }
  quantiles <- function(sims, p) {
    apply(sims, 2, function(v) {
      if (all(is.na(v))) NA_real_ else stats::quantile(v, p, na.rm = TRUE, names = FALSE)
    })
  }
  rows <- lapply(parts, function(part) {
    sims <- part$sims
    sims[!is.finite(sims)] <- NA
    tested <- site_tests(part$obs, sims, 1 - level)
    data.frame(site = x$site, trait = x$trait, statistic = part$statistic,
               observed = part$obs, null_median = quantiles(sims, 0.5),
               null_low = quantiles(sims, (1 - level) / 2),
               null_high = quantiles(sims, 1 - (1 - level) / 2),
               inner_low = quantiles(sims, (1 - inner) / 2),
               inner_high = quantiles(sims, 1 - (1 - inner) / 2),
               ses = tested$ses, p_low = tested$p_low, p_high = tested$p_high,
               deviation = tested$deviation, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Plot the null-model tests of each site
#'
#' Draws the tests of [tstats()] or [community_index()] site by site, in one panel per trait and
#' statistic. With `type = "intervals"` (the default) each site is a row with its null
#' distribution as two nested intervals (thin: `level`; thick: `inner`), the null median as a
#' tick and the observed value as a point coloured by the result of the test: below the null,
#' above it, or within it. Sites that cannot be tested (such as `T_IP.IC` at a site with one
#' species) are drawn as hollow points. With `type = "ses"` the standardised effect sizes are
#' drawn instead, with the band of
#' \eqn{\pm z_{1 - (1 - \mathrm{level})/2}}{+/- qnorm(1 - (1 - level) / 2)}; sites that
#' cannot be tested have no standardised effect size and are left empty. The colours still give
#' the result of the permutation test, which does not assume a normal null distribution. With
#' `hatch = TRUE` (the default) the area expected by chance (the outer null interval of each
#' site, or the band of standardised effect sizes) is hatched.
#'
#' The plot is built from [null_intervals()], and the ggplot object can be changed further with
#' ggplot2: for instance, `+ ggplot2::facet_wrap(~ statistic)`, `+ ggplot2::theme_bw()` or
#' `+ ggplot2::scale_x_log10()`.
#'
#' @param x Result of [tstats()] or [community_index()] with `n_perm > 0`.
#' @param type `"intervals"` or `"ses"`.
#' @param level,inner Coverage of the outer and inner null intervals; `1 - level` is also the
#'   size of the two-sided test (as in [tstats_summary()]).
#' @param statistics,traits,sites Subsets to draw; `NULL` (the default) draws all.
#' @param order Order of the sites within each panel's axis: `"data"` (as in `x`), `"ses"` (by
#'   the mean standardised effect size over the panels drawn) or `"observed"` (by the mean
#'   observed value).
#' @param colours Three colours for sites below the null, within it and above it. The defaults
#'   are blue, grey and vermillion of the Okabe-Ito palette.
#' @param point_size,interval_width Size of the observed points and width of the thick inner
#'   interval.
#' @param hatch Hatch the area expected by chance; with `FALSE` it is drawn as a thin line
#'   (`"intervals"`) or a grey band (`"ses"`).
#' @param hatch_density Number of hatching lines across the width of a panel.
#' @param scales Passed to [ggplot2::facet_grid()]; `"free_x"` (the default) lets each
#'   statistic keep its own range.
#' @param ... Not used.
#' @return A ggplot object.
#' @seealso [null_intervals()], [tstats_summary()], [plot.trait_flex()].
#' @family trait variation among individuals
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' set.seed(1)
#' ind <- data.frame(site = rep(paste0("s", 1:8), each = 20),
#'                   species = sample(c("a", "b", "c"), 160, replace = TRUE))
#' ind$size <- c(a = 10, b = 20, c = 30)[ind$species] + rnorm(160, sd = 3)
#' res <- tstats(ind, traits = "size", n_perm = 199)
#' plot(res)
#' plot(res, type = "ses", order = "ses")
#' # an index of communities
#' plot(community_index(ind, stats::var, traits = "size", null = "regional.ind", n_perm = 199))
#' @export
plot.tstats <- function(x, type = c("intervals", "ses"), level = 0.95, inner = 0.5,
                        statistics = NULL, traits = NULL, sites = NULL,
                        order = c("data", "ses", "observed"),
                        colours = c("#0072B2", "grey55", "#D55E00"), point_size = 2,
                        interval_width = 2.5, hatch = TRUE, hatch_density = 60,
                        scales = "free_x", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot() of null-model tests needs the package ggplot2", call. = FALSE) # nocov
  }
  type <- match.arg(type)
  order <- match.arg(order)
  check_colours(colours, 3)
  d <- null_intervals(x, level, inner)
  for (arg in c("statistics", "traits", "sites")) {
    keep <- get(arg)
    column <- c(statistics = "statistic", traits = "trait", sites = "site")[[arg]]
    if (!is.null(keep)) {
      unknown <- setdiff(keep, d[[column]])
      if (length(unknown) > 0) {
        stop("`", arg, "` has values not in `x`: ", paste(unknown, collapse = ", "),
             call. = FALSE)
      }
      d <- d[d[[column]] %in% keep, , drop = FALSE]
    }
  }
  site_levels <- unique(d$site)
  if (order != "data") {
    key <- tapply(d[[if (order == "ses") "ses" else "observed"]], d$site, mean, na.rm = TRUE)
    key <- key[site_levels]
    site_levels <- site_levels[order(key, na.last = FALSE)]
  }
  d$site <- factor(d$site, site_levels)
  d$trait <- factor(d$trait, unique(d$trait))
  d$statistic <- factor(d$statistic, unique(d$statistic))
  labels <- c(low = "Below the null", none = "Within the null", high = "Above the null")
  d$result <- factor(labels[ifelse(is.na(d$deviation), "none", d$deviation)], labels)
  d$tested <- !is.na(d$deviation)
  fill <- stats::setNames(colours, labels)
  stopifnot(is.logical(hatch), length(hatch) == 1, !is.na(hatch),
            is.numeric(hatch_density), length(hatch_density) == 1, hatch_density > 0)
  .data <- ggplot2::.data
  band <- stats::qnorm(1 - (1 - level) / 2)
  # spacing of the hatching lines: a fixed fraction of the range of each statistic (facet column)
  span_of <- if (type == "intervals") {
    tapply(c(d$null_low, d$null_high, d$observed), rep(d$statistic, 3),
           function(v) if (all(is.na(v))) NA_real_ else diff(range(v, na.rm = TRUE)))
  } else {
    tapply(c(d$ses, rep(c(-band, band), length.out = nrow(d))), rep(d$statistic, 2),
           function(v) if (all(is.na(v))) NA_real_ else diff(range(v, na.rm = TRUE)))
  }
  span_of[!is.finite(span_of) | span_of == 0] <- 1
  if (type == "intervals") {
    rects <- d[!is.na(d$null_low) & !is.na(d$null_high), , drop = FALSE]
    rects$ymin <- as.integer(rects$site) - 0.32
    rects$ymax <- as.integer(rects$site) + 0.32
    rects$xmin <- rects$null_low
    rects$xmax <- rects$null_high
  } else {
    rects <- unique(d[c("trait", "statistic")])
    rects$ymin <- 0.5
    rects$ymax <- nlevels(d$site) + 0.5
    rects$xmin <- -band
    rects$xmax <- band
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(y = .data$site))
  if (hatch) {
    lines <- hatch_lines(rects, span_of[as.character(rects$statistic)] / hatch_density)
    p <- p +
      ggplot2::geom_blank() +
      ggplot2::geom_rect(data = rects, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                                    ymin = .data$ymin, ymax = .data$ymax),
                         inherit.aes = FALSE, fill = NA, colour = "grey70", linewidth = 0.3) +
      ggplot2::geom_segment(data = lines, ggplot2::aes(x = .data$x, xend = .data$xend,
                                                       y = .data$y, yend = .data$yend),
                            inherit.aes = FALSE, colour = "grey70", linewidth = 0.25)
  }
  if (type == "intervals") {
    if (!hatch) {
      p <- p + ggplot2::geom_linerange(ggplot2::aes(xmin = .data$null_low,
                                                    xmax = .data$null_high),
                                       colour = "grey75", linewidth = 0.5, na.rm = TRUE)
    }
    p <- p +
      ggplot2::geom_linerange(ggplot2::aes(xmin = .data$inner_low, xmax = .data$inner_high),
                              colour = "grey75", linewidth = interval_width, na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(x = .data$null_median), shape = "|", size = 3,
                          colour = "grey40", na.rm = TRUE) +
      ggplot2::labs(x = "Observed value and null distribution", y = NULL)
    x_value <- "observed"
  } else {
    if (!hatch) {
      p <- p + ggplot2::annotate("rect", xmin = -band, xmax = band, ymin = -Inf, ymax = Inf,
                                 fill = "grey92")
    }
    p <- p +
      ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
      ggplot2::labs(x = "Standardised effect size", y = NULL)
    x_value <- "ses"
  }
  p +
    ggplot2::geom_point(data = d[d$tested, , drop = FALSE],
                        ggplot2::aes(x = .data[[x_value]], fill = .data$result), shape = 21,
                        size = point_size, colour = "white", stroke = 0.3, na.rm = TRUE,
                        show.legend = TRUE) +
    ggplot2::geom_point(data = d[!d$tested, , drop = FALSE],
                        ggplot2::aes(x = .data[[x_value]]), shape = 21, size = point_size,
                        colour = "grey40", fill = NA, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = fill, limits = names(fill), drop = FALSE, name = NULL) +
    ggplot2::facet_grid(ggplot2::vars(.data$trait), ggplot2::vars(.data$statistic),
                        scales = scales) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), legend.position = "top")
}

#' @rdname plot.tstats
#' @export
plot.community_index <- function(x, ...) plot.tstats(x, ...)

# diagonal hatching lines inside rectangles (columns xmin, xmax, ymin, ymax, and any facet
# columns, which are kept), spaced by `spacing` (one value per rectangle) along x; over one unit
# of y a line moves 2.5 spacings along x, so that the angle is the same in every rectangle
hatch_lines <- function(rects, spacing) {
  rows <- lapply(seq_len(nrow(rects)), function(i) {
    r <- rects[i, , drop = FALSE]
    h <- r$ymax - r$ymin
    shift <- 2.5 * spacing[i] * h
    starts <- seq(r$xmin - shift, r$xmax, by = spacing[i])
    # line: x = c + shift * (y - ymin) / h, clipped to the rectangle
    y_at <- function(x) r$ymin + h * (x - starts) / shift
    y0 <- pmax(r$ymin, y_at(r$xmin))
    y1 <- pmin(r$ymax, y_at(r$xmax))
    ok <- y1 > y0
    if (!any(ok)) return(NULL)
    x_at <- function(y) starts[ok] + shift * (y - r$ymin) / h
    cbind(r[rep(1, sum(ok)), setdiff(names(r), c("xmin", "xmax", "ymin", "ymax")),
            drop = FALSE],
          x = x_at(y0[ok]), xend = x_at(y1[ok]), y = y0[ok], yend = y1[ok])
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(x = numeric(0), xend = numeric(0), y = numeric(0), yend = numeric(0))
  }
  rownames(out) <- NULL
  out
}

#' Trait probability densities of populations and communities
#'
#' Estimates the trait probability density (TPD; Carmona et al. 2016, 2019) of each population
#' with a Gaussian kernel, and that of each community as the sum of the densities of its
#' populations weighted by their relative abundances, so that both integrate to one.
#'
#' All populations share one bandwidth, so that their densities can be compared and summed.
#' By default it follows Silverman's rule of thumb (as [stats::bw.nrd0()]) applied to the
#' deviations of individuals from the mean of their population, pooled over all populations
#' with two or more individuals (the standard deviation with \eqn{n_i - 1}{n_i - 1} degrees of
#' freedom per population, and the interquartile range of the deviations rescaled by
#' \eqn{\sqrt{n_i / (n_i - 1)}}{sqrt(n_i / (n_i - 1))}, so that small populations do not
#' shrink the bandwidth), with the mean number of individuals per population as the sample
#' size. A population of one individual then gets a kernel of that width around its
#' value, instead of no density.
#'
#' One bandwidth suits populations whose spread is about the same. For traits whose spread
#' grows with their mean, such as body size across species of very different sizes, the
#' bandwidth follows the many small populations and the densities of large ones break into one
#' peak per individual; the logarithm of the trait usually evens out the spread. With the
#' default bandwidth, a warning is given when the trait is positive and the standard deviation
#' of populations of three or more individuals rises with their mean (Spearman correlation
#' above 0.7 over at least ten populations).
#'
#' @param individuals Data frame with one row per individual.
#' @param trait Name of the column with the trait.
#' @param site,species Names of the columns with the site and the species.
#' @param abundance Name of a column with the abundance of each population (the same value
#'   for all its individuals), or `NULL` to weight populations by the number of individuals
#'   measured.
#' @param bw Bandwidth of the kernel, or `NULL` for the rule above.
#' @param n Number of points of the grid of trait values.
#' @return A data frame with the columns `site`, `species` (`NA` for communities), `level`
#'   (`"population"` or `"community"`), `weight` (relative abundance of the population in the
#'   site, 1 for communities), `trait` (grid value) and `density`. The bandwidth is in the
#'   attribute `bw`.
#' @references Carmona, C. P., de Bello, F., Mason, N. W. H., & Lepš, J. (2016). Traits without
#'   borders: integrating functional diversity across scales. *Trends in Ecology & Evolution*,
#'   31, 382--394. \doi{10.1016/j.tree.2016.02.003}
#'
#'   Carmona, C. P., de Bello, F., Mason, N. W. H., & Lepš, J. (2019). Trait probability density
#'   (TPD): measuring functional diversity across scales based on TPD with R. *Ecology*, 100,
#'   e02876. \doi{10.1002/ecy.2876}
#' @seealso [plot_population_traits()].
#' @family trait variation among individuals
#' @examples
#' set.seed(1)
#' ind <- data.frame(site = rep(c("s1", "s2"), each = 30),
#'                   species = sample(c("a", "b"), 60, replace = TRUE))
#' ind$size <- c(a = 10, b = 14)[ind$species] + rnorm(60)
#' d <- trait_densities(ind, "size")
#' # each density integrates to one
#' tapply(d$density, paste(d$site, d$species), sum) * diff(d$trait[1:2])
#' @export
trait_densities <- function(individuals, trait = "trait", site = "site", species = "species",
                            abundance = NULL, bw = NULL, n = 512) {
  stopifnot(is.data.frame(individuals), is.character(trait), length(trait) == 1,
            length(n) == 1, n >= 2)
  missing_columns <- setdiff(c(trait, site, species, abundance), names(individuals))
  if (length(missing_columns) > 0) {
    stop("`individuals` has no column ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  z <- individuals[[trait]]
  if (!is.numeric(z)) stop("trait `", trait, "` must be numeric", call. = FALSE)
  keep <- !is.na(z)
  z <- z[keep]
  site_id <- as.character(individuals[[site]])[keep]
  species_id <- as.character(individuals[[species]])[keep]
  pop <- paste(site_id, species_id, sep = "\r")
  pop_levels <- unique(pop)
  first <- match(pop_levels, pop)
  size <- tabulate(match(pop, pop_levels), length(pop_levels))
  weight <- if (is.null(abundance)) {
    size
  } else {
    a <- individuals[[abundance]][keep]
    if (anyNA(a) || any(tapply(a, pop, function(v) length(unique(v))) > 1)) {
      stop("`abundance` must have one value, not missing, for each population", call. = FALSE)
    }
    as.numeric(a[first])
  }
  pop_site <- site_id[first]
  weight <- weight / stats::ave(weight, pop_site, FUN = sum)
  if (is.null(bw)) {
    deviation <- z - stats::ave(z, pop)
    used <- size[match(pop, pop_levels)] > 1
    if (sum(used) < 2) {
      stop("no population has two or more individuals: give `bw`", call. = FALSE)
    }
    # deviations from the population mean underestimate the spread by sqrt((n - 1) / n):
    # the standard deviation is pooled with n - 1 degrees of freedom per population and the
    # interquartile range is taken on deviations rescaled by sqrt(n / (n - 1))
    n_i <- size[match(pop, pop_levels)][used]
    d <- deviation[used]
    pooled_sd <- sqrt(sum(d^2) / (length(d) - sum(size > 1)))
    spread <- min(pooled_sd, stats::IQR(d * sqrt(n_i / (n_i - 1))) / 1.34)
    if (!(spread > 0)) spread <- pooled_sd
    if (!(spread > 0)) {
      stop("the individuals of each population have identical values: give `bw`", call. = FALSE)
    }
    bw <- 0.9 * spread * mean(size)^(-1 / 5)
    check_spread_scale(z, pop)
  }
  stopifnot(is.numeric(bw), length(bw) == 1, bw > 0)
  grid <- seq(min(z) - 3 * bw, max(z) + 3 * bw, length.out = n)
  densities <- vapply(split(z, factor(pop, pop_levels)), function(v) {
    rowMeans(matrix(stats::dnorm(grid, rep(v, each = n), bw), n))
  }, numeric(n))
  densities <- matrix(densities, n)
  populations <- data.frame(site = rep(pop_site, each = n),
                            species = rep(species_id[first], each = n),
                            level = "population", weight = rep(weight, each = n),
                            trait = grid, density = as.vector(densities),
                            stringsAsFactors = FALSE)
  site_levels <- unique(site_id)
  community <- vapply(site_levels, function(s) {
    i <- pop_site == s
    as.vector(densities[, i, drop = FALSE] %*% weight[i])
  }, numeric(n))
  communities <- data.frame(site = rep(site_levels, each = n), species = NA_character_,
                            level = "community", weight = 1, trait = grid,
                            density = as.vector(community), stringsAsFactors = FALSE)
  out <- rbind(populations, communities)
  attr(out, "bw") <- bw
  out
}

# warns when the spread within populations rises with their mean, so that one bandwidth fits
# small and large populations badly
check_spread_scale <- function(z, pop) {
  if (any(z <= 0)) return(invisible(NULL))
  size <- tapply(z, pop, length)
  keep <- names(size)[size >= 3]
  if (length(keep) < 10) return(invisible(NULL))
  centre <- tapply(z, pop, mean)[keep]
  spread <- tapply(z, pop, stats::sd)[keep]
  if (stats::sd(spread) == 0) return(invisible(NULL))
  r <- stats::cor(centre, spread, method = "spearman")
  if (r > 0.7) {
    warning("the spread within populations rises with their mean (Spearman r = ",
            format(round(r, 2), nsmall = 2), "), so one bandwidth fits small and large ",
            "populations badly; consider the logarithm of the trait, or give `bw`",
            call. = FALSE)
  }
  invisible(r)
}

#' Plot the trait distributions of populations within communities
#'
#' Draws, in one panel per site, the trait probability density of each population from
#' [trait_densities()], filled with the colour of its species and with adjustable transparency,
#' the density of the community as a black line and the individuals as a rug. With
#' `weighted = TRUE` (the default) each population is scaled by its relative abundance, so that
#' the populations of a site add up to the community line; with `weighted = FALSE` every
#' population integrates to one, which shows the shape of rare populations.
#'
#' With many species, `highlight` draws a few of them in colour and the rest in grey. For traits
#' whose spread grows with their mean, such as body size, see the note on the bandwidth in
#' [trait_densities()].
#'
#' The ggplot object can be changed further, for instance with `+ ggplot2::facet_wrap(~ site,
#' ncol = 2)`, `+ ggplot2::scale_fill_viridis_d()` or `+ ggplot2::theme_classic()`.
#'
#' @inheritParams trait_densities
#' @param sites Sites to draw; `NULL` (the default) draws all.
#' @param weighted Scale the population densities by their relative abundances.
#' @param community Draw the density of the community.
#' @param rug Draw the individuals as a rug.
#' @param alpha Transparency of the population densities (0, transparent, to 1, opaque).
#' @param colours Colours of the species: unnamed, in the order of first appearance in
#'   `individuals`, or named by species. `NULL` (the default) uses the Okabe-Ito palette for
#'   up to seven species and [grDevices::hcl.colors()] (`"Dark 3"`) for more. With
#'   `highlight`, the colours of the highlighted species only.
#' @param highlight Species to draw in colour, on top of the others, which are drawn in grey
#'   and left out of the legend. Useful when a site holds many species, whose colours cannot
#'   all be told apart.
#' @param ncol,scales Passed to [ggplot2::facet_wrap()].
#' @return A ggplot object.
#' @seealso [trait_densities()], [tstats()].
#' @family trait variation among individuals
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' set.seed(1)
#' ind <- data.frame(site = rep(paste0("s", 1:4), each = 40),
#'                   species = sample(c("a", "b", "c"), 160, replace = TRUE))
#' ind$size <- c(a = 10, b = 14, c = 17)[ind$species] + rnorm(160)
#' plot_population_traits(ind, "size")
#' plot_population_traits(ind, "size", weighted = FALSE, alpha = 0.2, sites = c("s1", "s2"))
#' @export
plot_population_traits <- function(individuals, trait = "trait", site = "site",
                                   species = "species", abundance = NULL, sites = NULL,
                                   weighted = TRUE, community = TRUE, rug = TRUE,
                                   alpha = 0.35, colours = NULL, highlight = NULL,
                                   bw = NULL, n = 512, ncol = NULL, scales = "free_y") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_population_traits() needs the package ggplot2", call. = FALSE) # nocov
  }
  stopifnot(is.numeric(alpha), length(alpha) == 1, alpha >= 0, alpha <= 1)
  d <- trait_densities(individuals, trait, site, species, abundance, bw, n)
  points <- individuals[!is.na(individuals[[trait]]), , drop = FALSE]
  points <- data.frame(site = as.character(points[[site]]),
                       species = as.character(points[[species]]), trait = points[[trait]],
                       stringsAsFactors = FALSE)
  site_levels <- unique(points$site)
  if (!is.null(sites)) {
    unknown <- setdiff(sites, site_levels)
    if (length(unknown) > 0) {
      stop("`sites` has values not in the data: ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
    site_levels <- site_levels[site_levels %in% sites]
  }
  species_levels <- unique(points$species)
  if (!is.null(highlight)) {
    stopifnot(is.character(highlight), length(highlight) >= 1, !anyDuplicated(highlight))
    unknown <- setdiff(highlight, species_levels)
    if (length(unknown) > 0) {
      stop("`highlight` has species not in the data: ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
    if (is.null(colours)) {
      colours <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9",
                   "#F0E442")
      if (length(highlight) > length(colours)) {
        stop("give `colours` to highlight more than seven species", call. = FALSE)
      }
      colours <- colours[seq_along(highlight)]
    }
    if (is.null(names(colours))) {
      check_colours(colours, length(highlight))
      names(colours) <- highlight
    }
    if (anyNA(colours[highlight])) {
      stop("`colours` has no colour for some highlighted species", call. = FALSE)
    }
    others <- setdiff(species_levels, highlight)
    colours <- c(colours[highlight], stats::setNames(rep("grey75", length(others)), others))
  }
  if (is.null(colours)) {
    okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9", "#F0E442")
    colours <- if (length(species_levels) <= length(okabe_ito)) {
      okabe_ito[seq_along(species_levels)]
    } else {
      grDevices::hcl.colors(length(species_levels), "Dark 3")
    }
  }
  if (is.null(names(colours))) {
    check_colours(colours, length(species_levels))
    names(colours) <- species_levels
  } else {
    missing_species <- setdiff(species_levels, names(colours))
    if (length(missing_species) > 0) {
      stop("`colours` has no colour for ", paste(missing_species, collapse = ", "),
           call. = FALSE)
    }
    colours <- colours[species_levels]
    check_colours(unname(colours), length(species_levels))
  }
  d <- d[d$site %in% site_levels, , drop = FALSE]
  points <- points[points$site %in% site_levels, , drop = FALSE]
  d$site <- factor(d$site, site_levels)
  points$site <- factor(points$site, site_levels)
  pops <- d[d$level == "population", , drop = FALSE]
  if (weighted) pops$density <- pops$density * pops$weight
  pops$species <- factor(pops$species, species_levels)
  # highlighted species are drawn last, on top of the others
  draw_order <- if (is.null(highlight)) species_levels else
    c(setdiff(species_levels, highlight), highlight)
  pops$group <- factor(paste(pops$site, pops$species),
                       unique(paste(rep(site_levels, each = length(draw_order)),
                                    rep(draw_order, length(site_levels)))))
  points$species <- factor(points$species, species_levels)
  .data <- ggplot2::.data
  p <- ggplot2::ggplot(pops, ggplot2::aes(x = .data$trait, y = .data$density)) +
    ggplot2::geom_area(ggplot2::aes(fill = .data$species, colour = .data$species,
                                    group = .data$group),
                       position = "identity", alpha = alpha, linewidth = 0.3)
  if (community) {
    p <- p + ggplot2::geom_line(data = d[d$level == "community", , drop = FALSE],
                                colour = "black", linewidth = 0.6)
  }
  if (rug) {
    p <- p + ggplot2::geom_rug(data = points, ggplot2::aes(x = .data$trait,
                                                           colour = .data$species),
                               inherit.aes = FALSE, sides = "b", alpha = 0.8,
                               length = ggplot2::unit(0.04, "npc"))
  }
  breaks <- if (is.null(highlight)) ggplot2::waiver() else highlight
  p +
    ggplot2::scale_fill_manual(values = colours, breaks = breaks, name = "Species") +
    ggplot2::scale_colour_manual(values = colours, breaks = breaks, name = "Species") +
    ggplot2::facet_wrap(ggplot2::vars(.data$site), ncol = ncol, scales = scales) +
    ggplot2::labs(x = trait, y = if (weighted) "Density (weighted by abundance)" else
      "Density") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
