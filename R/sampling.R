# Sampling schemes for individuals of real data, and rarefaction of statistics over the number
# of individuals measured per population.

#' Subsamples of individuals
#'
#' Draws a subsample of the individuals of each population of real data, with the sampling
#' schemes of [sample_individuals()]: a fixed number of individuals per population (`n`), a
#' share of the individuals (`fraction`), only the most abundant species of each site
#' (`cover`), and a minimum number of individuals per population (`min_n`). It shows how much a
#' result depends on sampling effort, alone or through [trait_rarefaction()].
#'
#' The schemes are applied in this order: species outside `cover` are removed, individuals are
#' drawn without replacement (`n` from each population, all of them when the population is
#' smaller, or each one with probability `fraction`), and populations left with fewer than
#' `min_n` individuals are removed.
#'
#' @param individuals Data frame with one row per individual.
#' @param n Number of individuals drawn from each population, or `NULL` to keep all.
#' @param fraction Probability that each individual is kept; used when `n` is `NULL`.
#' @param cover Share (0 to 1) of the abundance of each site covered by the species kept, from
#'   the most abundant down (Pakeman & Quested 2007), or `NULL` to keep every species.
#'   Abundances are those of the `abundance` column, or the numbers of individuals in the data.
#'   Species tied in abundance are taken in the order of the data, and populations with an
#'   abundance of zero are never kept.
#' @param min_n Minimum number of individuals of a population after sampling.
#' @param site,species Names of the columns with the site and the species.
#' @param abundance Name of a column with the abundance of each population (the same value for
#'   all its individuals), or `NULL` to use the numbers of individuals in the data.
#' @return The rows of `individuals` kept, in their original order.
#' @references Pakeman, R. J., & Quested, H. M. (2007). Sampling plant functional traits: what
#'   proportion of the species need to be measured? *Applied Vegetation Science*, 10, 91--96.
#'   \doi{10.1111/j.1654-109x.2007.tb00507.x}
#' @seealso [trait_rarefaction()], [sample_individuals()].
#' @examples
#' set.seed(1)
#' ind <- data.frame(site = rep(c("s1", "s2"), each = 50),
#'                   species = sample(c("a", "b", "c"), 100, TRUE, prob = c(6, 3, 1)),
#'                   trait = rnorm(100))
#' table(subsample_individuals(ind, n = 5)[c("site", "species")])
#' # only the species making up 80% of the individuals of each site
#' table(subsample_individuals(ind, cover = 0.8)[c("site", "species")])
#' @export
subsample_individuals <- function(individuals, n = NULL, fraction = NULL, cover = NULL,
                                  min_n = 1, site = "site", species = "species",
                                  abundance = NULL) {
  stopifnot(is.data.frame(individuals))
  missing_columns <- setdiff(c(site, species, abundance), names(individuals))
  if (length(missing_columns) > 0) {
    stop("`individuals` has no column ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  stopifnot(length(min_n) == 1, min_n >= 1, min_n == round(min_n))
  site_id <- as.character(individuals[[site]])
  pop <- paste(site_id, as.character(individuals[[species]]), sep = "\r")
  pop_levels <- unique(pop)
  pop_index <- match(pop, pop_levels)
  first <- match(seq_along(pop_levels), pop_index)
  keep <- rep(TRUE, nrow(individuals))
  if (!is.null(cover)) {
    size <- if (is.null(abundance)) {
      tabulate(pop_index, length(pop_levels))
    } else {
      a <- individuals[[abundance]]
      if (anyNA(a) || any(tapply(a, pop_index, function(v) length(unique(v))) > 1)) {
        stop("`abundance` must have one value, not missing, for each population",
             call. = FALSE)
      }
      as.numeric(a[first])
    }
    dominant <- dominant_species(size, match(site_id[first], unique(site_id)), cover)
    keep <- dominant[pop_index]
  }
  if (!is.null(n)) {
    stopifnot(length(n) == 1, n >= 1, n == round(n))
    chosen <- unlist(lapply(split(which(keep), pop_index[keep]), function(i) {
      if (length(i) <= n) i else i[sample.int(length(i), n)]
    }), use.names = FALSE)
    keep <- seq_along(keep) %in% chosen
  } else if (!is.null(fraction)) {
    stopifnot(length(fraction) == 1, fraction > 0, fraction <= 1)
    keep <- keep & stats::runif(length(keep)) < fraction
  }
  size_kept <- tabulate(pop_index[keep], length(pop_levels))
  keep <- keep & size_kept[pop_index] >= min_n
  individuals[keep, , drop = FALSE]
}

#' Rarefaction of a statistic over sampling effort
#'
#' Computes a statistic of the individuals, such as the mean T-statistics of [tstats()] or
#' the decomposition of [trait_flex()], on repeated subsamples of decreasing effort, with the
#' schemes of [subsample_individuals()]. By default effort is the number of individuals per
#' population: `n` individuals are drawn from each population (populations smaller than `n`
#' keep all their individuals), for a series of `n`. With `cover`, effort is instead the share
#' of the abundance of each site covered by the species measured, from the most abundant down;
#' all individuals of those species are kept, or `n` of them when a single `n` is given. The
#' curve shows how much effort is needed before the statistic stops changing. This rarefies
#' individuals within populations or species within sites; the rarefaction of species richness
#' (for instance `vegan::rarefy()` or the package iNEXT) answers another question.
#'
#' @param individuals Data frame with one row per individual.
#' @param statistic Function of a data frame of individuals that returns a named numeric
#'   vector. The default gives the mean over sites of the three T-statistics of [tstats()] for
#'   the column `trait`.
#' @param n Numbers of individuals per population; the default runs from 2 to the median size
#'   of the populations. With `cover`, a single number of individuals per population, or
#'   `NULL` (the default) to keep all of them.
#' @param cover Shares (0 to 1) of abundance covered by the species measured, for a rarefaction
#'   over species instead of individuals; `NULL` (the default) for a rarefaction over `n`.
#' @param abundance Name of a column with the abundance of each population, used with `cover`
#'   (passed to [subsample_individuals()]).
#' @param reps Number of subsamples for each level of effort. With `cover` and no `n` the
#'   subsamples are all equal, and a single one is computed.
#' @param level Coverage of the interval of the statistic over subsamples.
#' @param min_n,site,species Passed to [subsample_individuals()]. With `min_n = 2` (the
#'   default) populations of one individual are left out, as in [tstats()] analyses.
#' @return A data frame of class `trait_rarefaction` with the columns `effort` (`n` or
#'   `cover`, as the attribute `by` says), `statistic` (the names returned by `statistic`),
#'   `mean`, `low` and `high` (quantiles over subsamples) and `reps` (subsamples with a finite
#'   value), with a [plot()][plot.trait_rarefaction] method. The values of all subsamples are
#'   in the attribute `values`.
#' @seealso [subsample_individuals()], [tstats()].
#' @examples
#' set.seed(1)
#' ind <- data.frame(site = rep(paste0("s", 1:5), each = 40),
#'                   species = sample(c("a", "b", "c"), 200, replace = TRUE))
#' ind$trait <- c(a = 10, b = 12, c = 15)[ind$species] + rnorm(200)
#' rare <- trait_rarefaction(ind, n = c(2, 4, 8), reps = 10)
#' rare
#' # over the share of abundance covered by the species measured
#' trait_rarefaction(ind, cover = c(0.5, 0.8, 1))
#' @export
trait_rarefaction <- function(individuals, statistic = NULL, n = NULL, cover = NULL,
                              reps = 20, level = 0.95, min_n = 2, site = "site",
                              species = "species", abundance = NULL) {
  stopifnot(is.data.frame(individuals), length(reps) == 1, reps >= 1,
            is.numeric(level), length(level) == 1, level > 0, level < 1)
  if (is.null(statistic)) {
    statistic <- function(d) {
      t <- tstats(d, site = site, species = species, nperm = 0)
      colMeans(t[c("T_IP.IC", "T_IC.IR", "T_PC.PR")], na.rm = TRUE)
    }
  }
  stopifnot(is.function(statistic))
  by <- if (is.null(cover)) "n" else "cover"
  if (by == "n") {
    if (is.null(n)) {
      size <- table(paste(individuals[[site]], individuals[[species]], sep = "\r"))
      n <- seq(2, max(2, stats::median(size)))
    }
    stopifnot(is.numeric(n), all(n >= 1), all(n == round(n)))
    effort <- sort(unique(n))
  } else {
    stopifnot(is.numeric(cover), all(cover > 0 & cover <= 1))
    if (!is.null(n)) stopifnot(length(n) == 1, n >= 1, n == round(n))
    if (is.null(n)) reps <- 1
    effort <- sort(unique(cover))
  }
  values <- lapply(effort, function(k) {
    v <- lapply(seq_len(reps), function(r) {
      d <- if (by == "n") {
        subsample_individuals(individuals, n = k, min_n = min_n, site = site,
                              species = species)
      } else {
        subsample_individuals(individuals, n = n, cover = k, min_n = min_n, site = site,
                              species = species, abundance = abundance)
      }
      out <- statistic(d)
      if (!is.numeric(out) || is.null(names(out)) || any(names(out) == "")) {
        stop("`statistic` must return a named numeric vector", call. = FALSE)
      }
      out
    })
    do.call(rbind, v)
  })
  names_out <- colnames(values[[1]])
  rows <- lapply(seq_along(effort), function(i) {
    m <- values[[i]]
    m[!is.finite(m)] <- NA
    data.frame(effort = effort[i], statistic = names_out, mean = colMeans(m, na.rm = TRUE),
               low = apply(m, 2, stats::quantile, (1 - level) / 2, na.rm = TRUE, names = FALSE),
               high = apply(m, 2, stats::quantile, 1 - (1 - level) / 2, na.rm = TRUE,
                            names = FALSE),
               reps = colSums(!is.na(m)), row.names = NULL, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$mean[out$reps == 0] <- NA_real_
  attr(out, "values") <- stats::setNames(values, effort)
  attr(out, "by") <- by
  class(out) <- c("trait_rarefaction", class(out))
  out
}

#' Plot a rarefaction of a statistic
#'
#' Draws the mean of each statistic from [trait_rarefaction()] against sampling effort, with
#' the interval over subsamples as a band, one panel per statistic. Dashed vertical lines mark
#' reference levels of effort: for a rarefaction over `cover`, by default the published rule of
#' measuring at least the species that make up 80% of the vegetation (Pakeman & Quested 2007).
#' No such rule is drawn by default for the number of individuals per population; give it in
#' `reference` when the field or study has one.
#'
#' @param x Result of [trait_rarefaction()].
#' @param reference Levels of effort marked by dashed lines; named values are written next to
#'   the lines. `NULL` (the default) marks the 80% rule for `cover` and nothing for `n`;
#'   `FALSE` draws no line.
#' @param colour Colour of the line and band; the default is the blue of the Okabe-Ito palette.
#' @param alpha Transparency of the band.
#' @param scales Passed to [ggplot2::facet_wrap()].
#' @param ... Not used.
#' @return A ggplot object.
#' @references Pakeman, R. J., & Quested, H. M. (2007). Sampling plant functional traits: what
#'   proportion of the species need to be measured? *Applied Vegetation Science*, 10, 91--96.
#'   \doi{10.1111/j.1654-109x.2007.tb00507.x}
#' @seealso [trait_rarefaction()].
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' set.seed(1)
#' ind <- data.frame(site = rep(paste0("s", 1:5), each = 40),
#'                   species = sample(c("a", "b", "c", "d"), 200, TRUE, prob = c(5, 3, 1, 1)))
#' ind$trait <- c(a = 10, b = 12, c = 15, d = 18)[ind$species] + rnorm(200)
#' plot(trait_rarefaction(ind, n = 2:10, reps = 10))
#' plot(trait_rarefaction(ind, cover = seq(0.3, 1, by = 0.1)))
#' # a rule of the study, for the number of individuals
#' plot(trait_rarefaction(ind, n = 2:10, reps = 10), reference = c("protocol" = 5))
#' @export
plot.trait_rarefaction <- function(x, reference = NULL, colour = "#0072B2", alpha = 0.25,
                                   scales = "free_y", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot() of trait_rarefaction() needs the package ggplot2", call. = FALSE) # nocov
  }
  check_colours(colour, 1)
  by <- attr(x, "by")
  if (is.null(by)) by <- "n"
  if (is.null(reference)) {
    reference <- if (by == "cover") c("80% rule (Pakeman & Quested 2007)" = 0.8) else numeric(0)
  } else if (isFALSE(reference)) {
    reference <- numeric(0)
  }
  stopifnot(is.numeric(reference))
  d <- as.data.frame(unclass(x))
  d$statistic <- factor(d$statistic, unique(d$statistic))
  .data <- ggplot2::.data
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$effort, y = .data$mean))
  if (length(reference) > 0) {
    lines <- data.frame(effort = unname(reference),
                        label = if (is.null(names(reference))) "" else names(reference))
    p <- p +
      ggplot2::geom_vline(data = lines, ggplot2::aes(xintercept = .data$effort),
                          linetype = "dashed", colour = "grey35") +
      ggplot2::geom_text(data = lines, ggplot2::aes(x = .data$effort, y = Inf,
                                                    label = .data$label),
                         inherit.aes = FALSE, angle = 90, hjust = 1.05, vjust = -0.4,
                         size = 3, colour = "grey35")
  }
  p +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$low, ymax = .data$high), fill = colour,
                         alpha = alpha, na.rm = TRUE) +
    ggplot2::geom_line(colour = colour, na.rm = TRUE) +
    ggplot2::geom_point(colour = colour, size = 1.5, na.rm = TRUE) +
    ggplot2::facet_wrap(ggplot2::vars(.data$statistic), scales = scales) +
    ggplot2::scale_x_continuous(labels = if (by == "cover") function(v) paste0(100 * v, "%")
                                else ggplot2::waiver()) +
    ggplot2::labs(x = if (by == "cover") "Share of abundance covered by the species measured"
                  else "Individuals per population", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

# TRUE for the populations of the most abundant species of each site that together reach
# `cover` of its abundance (the species that crosses the threshold is included)
dominant_species <- function(abundance, site, cover) {
  stopifnot(is.numeric(cover), length(cover) == 1, cover > 0, cover <= 1)
  keep <- logical(length(abundance))
  for (i in split(seq_along(abundance), site)) {
    i <- i[abundance[i] > 0]
    if (length(i) == 0) next
    o <- i[order(-abundance[i])]
    share <- cumsum(abundance[o]) / sum(abundance[o])
    before <- c(0, share[-length(share)])
    keep[o[before < cover - 1e-12]] <- TRUE
  }
  keep
}
