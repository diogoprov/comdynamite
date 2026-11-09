#' Trait-gradient analysis
#'
#' Splits the trait value of each species into a position along the gradient of community
#' trait means (\eqn{\beta}) and a difference from the species it co-occurs with
#' (\eqn{\alpha}), as in the trait-gradient analysis of Ackerly & Cornwell (2007), with three
#' changes: the mean of each site leaves out the focal species, the weights of the site means
#' and of the species means are chosen separately, and [trait_gradient_test()] and
#' [trait_gradient_phylo()] give a null model that keeps the richness of sites and the
#' frequency of species and a phylogenetic model that carries the error of \eqn{\beta} and
#' \eqn{\alpha}.
#'
#' @details
#' With \eqn{t_{ij}} the trait value of species \eqn{i} in site \eqn{j} and \eqn{a_{ij}} its
#' weight, the mean of site \eqn{j} is
#' \deqn{p_j = \frac{\sum_i a_{ij} t_{ij}}{\sum_i a_{ij}},}
#' the mean of species \eqn{i} is \eqn{\bar t_i = \sum_j a_{ij} t_{ij} / \sum_j a_{ij}}, its
#' position along the gradient is
#' \deqn{\beta_i = \frac{\sum_j a_{ij} p_j}{\sum_j a_{ij}}}
#' and its difference from the co-occurring species is \eqn{\alpha_i = \bar t_i - \beta_i}.
#' The niche breadth \eqn{R_i} is the range of \eqn{p_j} over the sites of the species and
#' the slope \eqn{b_i} is the weighted regression of \eqn{t_{ij}} on \eqn{p_j}, the response
#' of the species to the gradient within its range.
#'
#' @section The three changes:
#' **Leave-one-out site means.** In the original analysis \eqn{p_j} includes the focal species,
#' so \eqn{\beta_i} contains \eqn{\bar t_i}: in sites with few species, \eqn{\beta} follows the
#' trait of the species itself and the correlation between \eqn{\beta} and \eqn{\bar t}, and
#' between \eqn{\beta} and \eqn{\alpha}, is inflated. With `leave_one_out = TRUE` (the default)
#' \eqn{\beta_i}, \eqn{R_i} and \eqn{b_i} use
#' \eqn{p_j^{(-i)} = (\sum_{k \ne i} a_{kj} t_{kj}) / \sum_{k \ne i} a_{kj}}, and sites where
#' species \eqn{i} is alone are left out of its values (counted in `n_alone`). The jackknife of
#' the original script does the same, as an option; `leave_one_out = FALSE` gives the original
#' values.
#'
#' **Weights.** The weights of the site means (`site_weights`) and of the species means and
#' \eqn{\beta} (`species_weights`) are chosen apart, as abundance or occurrence; the original
#' uses one choice for both. Compare the two with [trait_gradient_weights()].
#'
#' **Null model and error.** See [trait_gradient_test()] and [trait_gradient_phylo()].
#'
#' Trait values missing in a site where the species occurs (`NA` in `trait`, with the row kept)
#' are replaced by the mean of the species, as in the original; the slope \eqn{b_i} uses only
#' measured values.
#'
#' @param data Data frame with one row per individual or per population (species in a site),
#'   such as from [sample_individuals()]. Rows of one population are averaged.
#' @param trait Name of the column with trait values.
#' @param site,species Names of the columns with the site and the species.
#' @param abundance Name of the column with the abundance of the species in the site (one value
#'   per population), or `NULL` to weight each occurrence equally.
#' @param site_weights,species_weights `"abundance"` or `"occurrence"`: weights of the site
#'   means and of the species means and \eqn{\beta}. Both default to `"abundance"` when
#'   `abundance` is given.
#' @param leave_one_out Leave the focal species out of the site means (the default).
#' @param min_slope Smallest number of sites with measured values to fit the slope \eqn{b_i}.
#' @return An object of class `mc_trait_gradient`, a list with
#'   \describe{
#'     \item{species}{one row per species: `species`, `n_sites`, `n_alone` (sites where the
#'       species is alone), `mean` (\eqn{\bar t_i}), `beta`, `alpha`, `breadth` (\eqn{R_i}),
#'       `slope` (\eqn{b_i}), `slope_se`, `intercept` and `min_site`, `max_site` (range of the
#'       site means of its sites).}
#'     \item{sites}{one row per site: `site`, `n_species` and `mean` (\eqn{p_j}, with all
#'       species).}
#'     \item{populations}{one row per population: `site`, `species`, `weight`, `trait`
#'       (with species means filled in), `measured` and `site_mean` (\eqn{p_j}, or
#'       \eqn{p_j^{(-i)}}).}
#'     \item{settings}{the arguments.}
#'   }
#' @references Ackerly, D. D., & Cornwell, W. K. (2007). A trait-based approach to community
#'   assembly: partitioning of species trait values into within- and among-community
#'   components. *Ecology Letters*, 10, 135--145. \doi{10.1111/j.1461-0248.2006.01006.x}
#' @seealso [trait_gradient_test()], [trait_gradient_phylo()], [trait_gradient_weights()],
#'   [plot.mc_trait_gradient()].
#' @family trait variation among individuals
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(20)
#' sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
#'                               competition_matrix(tree, "none"), n_steps = 100,
#'                               seed_steps = 1:5)
#' ind <- sample_individuals(sim, n = 5, local_adaptation = 0.5, intraspecific_sd = 1)
#' tg <- trait_gradient(ind, abundance = "abundance")
#' head(tg$species)
#' # the original, with the focal species in the site means
#' original <- trait_gradient(ind, abundance = "abundance", leave_one_out = FALSE)
#' c(leave_one_out = cor(tg$species$beta, tg$species$alpha),
#'   original = cor(original$species$beta, original$species$alpha))
#' @export
trait_gradient <- function(data, trait = "trait", site = "site", species = "species",
                           abundance = NULL,
                           site_weights = if (is.null(abundance)) "occurrence" else "abundance",
                           species_weights = site_weights, leave_one_out = TRUE,
                           min_slope = 3) {
  site_weights <- match.arg(site_weights, c("abundance", "occurrence"))
  species_weights <- match.arg(species_weights, c("abundance", "occurrence"))
  pops <- tg_populations(data, trait, site, species, abundance)
  if ("abundance" %in% c(site_weights, species_weights) && is.null(abundance)) {
    stop("weights by abundance need the column `abundance`", call. = FALSE)
  }
  stopifnot(is.logical(leave_one_out), length(leave_one_out) == 1, !is.na(leave_one_out))
  check_count(min_slope, "min_slope", min = 2)
  settings <- list(trait = trait, site_weights = site_weights,
                   species_weights = species_weights, leave_one_out = leave_one_out,
                   min_slope = min_slope)
  out <- tg_compute(pops, settings)
  out$settings <- settings
  class(out) <- "mc_trait_gradient"
  out
}

# one row per population (site x species): site, species, abundance, trait (NA if not
# measured); several rows of a population are averaged
tg_populations <- function(data, trait, site, species, abundance) {
  if (!is.data.frame(data)) stop("`data` must be a data frame", call. = FALSE)
  needed <- c(trait, site, species, abundance)
  missing_columns <- setdiff(needed, names(data))
  if (length(missing_columns) > 0) {
    stop("`data` has no column ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(data[[trait]])) stop("`", trait, "` must be numeric", call. = FALSE)
  key <- paste(data[[site]], data[[species]], sep = "\r")
  first <- !duplicated(key)
  index <- match(key, key[first])
  values <- tapply(data[[trait]], index, function(v) if (all(is.na(v))) NA_real_ else
    mean(v, na.rm = TRUE))
  pops <- data.frame(site = as.character(data[[site]][first]),
                     species = as.character(data[[species]][first]),
                     abundance = 1, trait = as.vector(values), stringsAsFactors = FALSE)
  if (!is.null(abundance)) {
    a <- data[[abundance]]
    if (!is.numeric(a) || anyNA(a) || any(a < 0)) {
      stop("`", abundance, "` must be numbers, 0 or more, without missing values",
           call. = FALSE)
    }
    differs <- tapply(a, index, function(v) diff(range(v)) > 0)
    if (any(differs)) {
      stop("`", abundance, "` must have one value per population (species in a site)",
           call. = FALSE)
    }
    pops$abundance <- a[first]
    pops <- pops[pops$abundance > 0, , drop = FALSE]
  }
  if (all(is.na(pops$trait))) stop("`", trait, "` has no values", call. = FALSE)
  unmeasured <- tapply(!is.na(pops$trait), pops$species, any)
  if (!all(unmeasured)) {
    stop("species without any trait value: ",
         paste(names(unmeasured)[!unmeasured], collapse = ", "), call. = FALSE)
  }
  if (length(unique(pops$site)) < 2) stop("`data` needs two sites or more", call. = FALSE)
  rownames(pops) <- NULL
  pops
}

# the analysis of a table of populations (site, species, abundance, trait)
tg_compute <- function(pops, settings, slopes = TRUE) {
  measured <- !is.na(pops$trait)
  w_site <- if (settings$site_weights == "abundance") pops$abundance else rep(1, nrow(pops))
  w_species <- if (settings$species_weights == "abundance") pops$abundance else
    rep(1, nrow(pops))
  # species means with weights of species, used to fill unmeasured populations
  fill <- tapply(w_species[measured] * pops$trait[measured], pops$species[measured], sum) /
    tapply(w_species[measured], pops$species[measured], sum)
  t_full <- ifelse(measured, pops$trait, fill[pops$species])
  site_sum <- tapply(w_site * t_full, pops$site, sum)
  site_w <- tapply(w_site, pops$site, sum)
  site_mean <- site_sum / site_w
  if (settings$leave_one_out) {
    rest <- site_w[pops$site] - w_site
    p <- ifelse(rest > 0, (site_sum[pops$site] - w_site * t_full) / rest, NA_real_)
  } else {
    p <- site_mean[pops$site]
  }
  p <- as.vector(p)
  species_ids <- sort(unique(pops$species))
  rows <- split(seq_len(nrow(pops)), factor(pops$species, species_ids))
  stats <- t(vapply(rows, function(r) {
    w <- w_species[r]
    tbar <- sum(w * t_full[r]) / sum(w)
    ok <- !is.na(p[r])
    beta <- if (any(ok)) sum(w[ok] * p[r][ok]) / sum(w[ok]) else NA_real_
    range_p <- if (any(ok)) range(p[r][ok]) else c(NA_real_, NA_real_)
    fit <- c(NA_real_, NA_real_, NA_real_)
    use <- ok & measured[r]
    if (slopes && sum(use) >= settings$min_slope && diff(range(p[r][use])) > 0) {
      x <- p[r][use]
      y <- pops$trait[r][use]
      ww <- w[use]
      fit <- weighted_slope(x, y, ww)
    }
    c(length(r), sum(!ok), tbar, beta, tbar - beta, range_p[2] - range_p[1], fit,
      range_p)
  }, numeric(11)))
  species_table <- data.frame(species = species_ids, stats, row.names = NULL)
  names(species_table) <- c("species", "n_sites", "n_alone", "mean", "beta", "alpha",
                            "breadth", "intercept", "slope", "slope_se", "min_site",
                            "max_site")
  species_table <- species_table[c("species", "n_sites", "n_alone", "mean", "beta", "alpha",
                                   "breadth", "slope", "slope_se", "intercept", "min_site",
                                   "max_site")]
  species_table[c("n_sites", "n_alone")] <- lapply(species_table[c("n_sites", "n_alone")],
                                                   as.integer)
  site_ids <- names(site_mean)
  sites <- data.frame(site = site_ids, n_species = as.integer(table(pops$site)[site_ids]),
                      mean = as.vector(site_mean), row.names = NULL)
  populations <- data.frame(site = pops$site, species = pops$species, weight = w_species,
                            trait = t_full, measured = measured, site_mean = p,
                            stringsAsFactors = FALSE)
  list(species = species_table, sites = sites, populations = populations)
}

# intercept, slope and standard error of the slope of a weighted least-squares line
weighted_slope <- function(x, y, w) {
  sw <- sum(w)
  mx <- sum(w * x) / sw
  my <- sum(w * y) / sw
  sxx <- sum(w * (x - mx)^2)
  slope <- sum(w * (x - mx) * (y - my)) / sxx
  intercept <- my - slope * mx
  n <- length(x)
  se <- if (n > 2) sqrt(sum(w * (y - intercept - slope * x)^2) / (n - 2) / sxx) else NA_real_
  c(intercept, slope, se)
}

#' @export
print.mc_trait_gradient <- function(x, ...) {
  s <- x$species
  cat("Trait-gradient analysis of", nrow(s), "species in", nrow(x$sites), "sites\n")
  cat("Site means:", if (x$settings$leave_one_out) "without the focal species" else
    "with the focal species", "| weights of sites:", x$settings$site_weights,
    "| of species:", x$settings$species_weights, "\n")
  ok <- stats::complete.cases(s$beta, s$alpha)
  cat("Correlation of beta and alpha:",
      format(signif(stats::cor(s$beta[ok], s$alpha[ok]), 3)), "\n")
  cat("Variance of species means: ", format(signif(stats::var(s$mean), 3)),
      "\n  of beta: ", format(signif(stats::var(s$beta, na.rm = TRUE), 3)),
      "\n  of alpha: ", format(signif(stats::var(s$alpha, na.rm = TRUE), 3)), "\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.mc_trait_gradient <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$species
}

#' Trait-gradient analysis with the two choices of weights
#'
#' Runs [trait_gradient()] with each combination of weights of the site means and of the
#' species means (abundance or occurrence) and puts the values of \eqn{\beta} and
#' \eqn{\alpha} of each species side by side, with the correlations among the combinations.
#' Large differences mean that the result depends on the dominant species.
#'
#' @inheritParams trait_gradient
#' @param abundance Name of the column with abundances (needed here).
#' @param ... Further arguments to [trait_gradient()].
#' @return A list with `beta` and `alpha`, data frames with one row per species and one column
#'   per combination (`abundance/abundance`, `abundance/occurrence`, `occurrence/abundance`,
#'   `occurrence/occurrence`, sites first), and `correlation`, the correlations of
#'   \eqn{\beta} and of \eqn{\alpha} between each combination and the first.
#' @seealso [trait_gradient()]
#' @family trait variation among individuals
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(20)
#' sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
#'                               competition_matrix(tree, "none"), n_steps = 100,
#'                               seed_steps = 1:5)
#' ind <- sample_individuals(sim, n = 5, intraspecific_sd = 1)
#' trait_gradient_weights(ind, abundance = "abundance")$correlation
#' @export
trait_gradient_weights <- function(data, abundance, ...) {
  if (missing(abundance) || is.null(abundance)) {
    stop("trait_gradient_weights() needs the column `abundance`", call. = FALSE)
  }
  choices <- expand.grid(site = c("abundance", "occurrence"),
                         species = c("abundance", "occurrence"), stringsAsFactors = FALSE)
  choices <- choices[order(choices$site), ]
  labels <- paste(choices$site, choices$species, sep = "/")
  fits <- Map(function(s, p) trait_gradient(data, abundance = abundance, site_weights = s,
                                            species_weights = p, ...)$species,
              choices$site, choices$species)
  pick <- function(column) {
    out <- data.frame(fits[[1]]$species, lapply(fits, `[[`, column), row.names = NULL)
    names(out) <- c("species", labels)
    out
  }
  beta <- pick("beta")
  alpha <- pick("alpha")
  correlation <- data.frame(weights = labels,
                            beta = vapply(labels, function(l)
                              stats::cor(beta[[labels[1]]], beta[[l]], use = "complete.obs"),
                              numeric(1)),
                            alpha = vapply(labels, function(l)
                              stats::cor(alpha[[labels[1]]], alpha[[l]],
                                         use = "complete.obs"), numeric(1)),
                            row.names = NULL)
  list(beta = beta, alpha = alpha, correlation = correlation)
}

#' Null model of the trait-gradient analysis
#'
#' Compares the values of [trait_gradient()] with those of communities in which the species
#' are moved among sites while the richness of each site and the number of sites of each
#' species stay as observed (the curveball algorithm of Strona et al. 2014). Each population
#' keeps its trait value and abundance when it moves, so the null breaks only the association
#' between the traits of species and the sites they share.
#'
#' @details
#' The null model of Ackerly & Cornwell (2007) permutes the site of each population
#' (`null = "permute"`), which changes the richness of sites and can place a species twice in
#' one site. Under that null, species of few sites have extreme \eqn{\beta} by chance and
#' rich sites pull every \eqn{\beta} toward the overall mean. The curveball null keeps both
#' margins, so the remaining signal comes from which species co-occur.
#'
#' The chain starts at the observed communities, runs `burn_in` trades and then `thin` trades
#' between null communities. Statistics: the correlation of \eqn{\beta} and \eqn{\alpha}, the
#' variance of \eqn{\beta} and the variance of \eqn{\alpha}, and the \eqn{\beta} and
#' \eqn{\alpha} of each species. Tests are two-sided.
#'
#' @param x Result of [trait_gradient()].
#' @param n_null Number of null communities.
#' @param null `"curveball"` (the default) or `"permute"`, the null of the original analysis.
#' @param burn_in,thin Trades of the curveball chain before the first null community and
#'   between null communities.
#' @return An object of class `mc_trait_gradient_test`, a list with `statistics` (one row per
#'   statistic: `observed`, `null_mean`, `null_sd`, `ses` and `p`), `species` (one row per
#'   species: observed \eqn{\beta} and \eqn{\alpha} with their `ses_` and `p_`), `null` (a
#'   matrix with one row per null community and one column per statistic) and the settings.
#' @references Ackerly, D. D., & Cornwell, W. K. (2007). A trait-based approach to community
#'   assembly: partitioning of species trait values into within- and among-community
#'   components. *Ecology Letters*, 10, 135--145. \doi{10.1111/j.1461-0248.2006.01006.x}
#'
#'   Strona, G., Nappo, D., Boccacci, F., Fattorini, S., & San-Miguel-Ayanz, J. (2014). A fast
#'   and unbiased procedure to randomize ecological binary matrices with fixed row and column
#'   totals. *Nature Communications*, 5, 4114. \doi{10.1038/ncomms5114}
#' @seealso [trait_gradient()]
#' @family trait variation among individuals
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(20)
#' sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
#'                               competition_matrix(tree, "none"), n_steps = 100,
#'                               seed_steps = 1:5)
#' ind <- sample_individuals(sim, n = 5, intraspecific_sd = 1)
#' tg <- trait_gradient(ind, abundance = "abundance")
#' trait_gradient_test(tg, n_null = 99)
#' @export
trait_gradient_test <- function(x, n_null = 999, null = c("curveball", "permute"),
                                burn_in = 1000, thin = 100) {
  if (!inherits(x, "mc_trait_gradient")) {
    stop("`x` must come from trait_gradient()", call. = FALSE)
  }
  null <- match.arg(null)
  check_count(n_null, "n_null")
  check_count(burn_in, "burn_in", min = 0)
  check_count(thin, "thin")
  pops <- tg_records(x)
  settings <- x$settings
  observed <- tg_null_stats(x)
  species_ids <- x$species$species
  k <- length(species_ids)
  site_of <- pops$site
  if (null == "curveball") site_of <- curveball_sites(site_of, pops$species, burn_in)
  sims <- matrix(NA_real_, n_null, 3 + 2 * k)
  for (b in seq_len(n_null)) {
    site_of <- if (null == "curveball") curveball_sites(site_of, pops$species, thin) else
      sample(pops$site)
    moved <- pops
    moved$site <- site_of
    sims[b, ] <- tg_null_stats(tg_compute(moved, settings, slopes = FALSE))
  }
  labels <- c("cor_beta_alpha", "var_beta", "var_alpha", paste0("beta_", species_ids),
              paste0("alpha_", species_ids))
  colnames(sims) <- labels
  compared <- null_compare(observed, sims, "two.sided")
  sd_null <- apply(sims, 2, stats::sd, na.rm = TRUE)
  statistics <- data.frame(statistic = labels[1:3], observed = observed[1:3],
                           null_mean = colMeans(sims[, 1:3], na.rm = TRUE),
                           null_sd = sd_null[1:3], ses = compared$ses[1:3],
                           p = compared$p[1:3], row.names = NULL)
  b_cols <- 3 + seq_len(k)
  a_cols <- 3 + k + seq_len(k)
  species <- data.frame(species = species_ids, beta = observed[b_cols],
                        ses_beta = compared$ses[b_cols], p_beta = compared$p[b_cols],
                        alpha = observed[a_cols], ses_alpha = compared$ses[a_cols],
                        p_alpha = compared$p[a_cols], row.names = NULL)
  structure(list(statistics = statistics, species = species, null = sims,
                 settings = c(settings, list(null = null, n_null = n_null, burn_in = burn_in,
                                             thin = thin))),
            class = "mc_trait_gradient_test")
}

# the populations of a trait_gradient() result as the input of tg_compute(), with
# unmeasured trait values back to NA and the abundance as weight
tg_records <- function(x) {
  p <- x$populations
  data.frame(site = p$site, species = p$species, abundance = p$weight,
             trait = ifelse(p$measured, p$trait, NA_real_), stringsAsFactors = FALSE)
}

# correlation of beta and alpha, their variances, and beta and alpha of each species
tg_null_stats <- function(x) {
  s <- x$species
  ok <- stats::complete.cases(s$beta, s$alpha)
  r <- if (sum(ok) > 2) stats::cor(s$beta[ok], s$alpha[ok]) else NA_real_
  c(r, stats::var(s$beta, na.rm = TRUE), stats::var(s$alpha, na.rm = TRUE), s$beta, s$alpha)
}

# `n_trades` curveball trades on the sites of populations; species never repeat in a site
curveball_sites <- function(site_of, species, n_trades) {
  if (n_trades == 0) return(site_of)
  site_ids <- unique(site_of)
  if (length(site_ids) < 2) return(site_of)
  for (i in seq_len(n_trades)) {
    pair <- sample(site_ids, 2)
    in_a <- which(site_of == pair[1])
    in_b <- which(site_of == pair[2])
    only_a <- in_a[!species[in_a] %in% species[in_b]]
    only_b <- in_b[!species[in_b] %in% species[in_a]]
    pool <- c(only_a, only_b)
    if (length(only_a) == 0 || length(only_b) == 0) next
    to_a <- pool[sample.int(length(pool), length(only_a))]
    site_of[pool] <- pair[2]
    site_of[to_a] <- pair[1]
  }
  site_of
}

#' @export
print.mc_trait_gradient_test <- function(x, ...) {
  s <- x$settings
  cat("Null model of the trait-gradient analysis:", s$n_null, "communities,", s$null, "\n")
  st <- x$statistics
  st[-1] <- lapply(st[-1], function(v) signif(v, 3))
  print(st, row.names = FALSE)
  invisible(x)
}

#' Phylogenetic models of beta and alpha, with their sampling error
#'
#' Fits models of trait evolution to the \eqn{\beta}, \eqn{\alpha} and mean trait of the
#' species from [trait_gradient()], with the sampling error of each value added to the
#' phylogenetic covariance, and compares a model without phylogeny, Brownian motion and an
#' Ornstein-Uhlenbeck process by AIC. It tells whether the position of species along the
#' gradient (\eqn{\beta}) and their difference from co-occurring species (\eqn{\alpha}) are
#' conserved or labile.
#'
#' @details
#' The standard error of each value comes from a bootstrap of the sites of each species: the
#' sites (with their trait values and abundances) are resampled with replacement, and
#' \eqn{\bar t_i}, \eqn{\beta_i} and \eqn{\alpha_i} recomputed with the site means of the
#' observed data (without the focal species, if so chosen). The site means are held fixed
#' because the question is how well the sites of a species measure its position, not how the
#' gradient itself would change.
#'
#' For each response \eqn{y} the model is \eqn{y \sim N(\mu 1, \sigma^2 C + S)}, with
#' \eqn{S} the diagonal matrix of squared standard errors and \eqn{C} the identity
#' (`"white"`), the shared branch lengths (`"BM"`), or, for `"OU"`, the covariance of an
#' Ornstein-Uhlenbeck process with a root at \eqn{\mu} and strength \eqn{a},
#' \deqn{C_{ij} = e^{-a(d_i + d_j - 2 s_{ij})} \frac{1 - e^{-2 a s_{ij}}}{2a},}
#' where \eqn{d_i} is the distance from the root to tip \eqn{i} and \eqn{s_{ij}} the length
#' shared by \eqn{i} and \eqn{j}. The parameters are estimated by restricted maximum
#' likelihood, so AIC compares the covariance models with the same mean. Without the error
#' term, species with few sites, whose values are noisy, would look like fast evolution.
#'
#' @param x Result of [trait_gradient()].
#' @param tree Phylogeny (`phylo`) with branch lengths; its tips are matched to species by
#'   name.
#' @param responses Which of `"beta"`, `"alpha"` and `"mean"` to model.
#' @param models Which of `"white"`, `"BM"` and `"OU"` to fit.
#' @param n_boot Number of bootstrap samples for the standard errors.
#' @param min_sites Smallest number of sites of a species with a value (sites where it is
#'   alone do not count); species with fewer sites, or not in `tree`, are left out.
#' @return An object of class `mc_trait_gradient_phylo`, a list with `models` (one row per
#'   response and model: `response`, `model`, `mu`, `sigma2`, `a`, `half_life` (\eqn{\ln 2 / a},
#'   in units of branch length), `loglik`, `n_par`, `aic` and `delta_aic` within each
#'   response), `data` (one row per species used: values and standard errors) and `tree`
#'   (the tree pruned to those species).
#' @references Ackerly, D. D., & Cornwell, W. K. (2007). A trait-based approach to community
#'   assembly: partitioning of species trait values into within- and among-community
#'   components. *Ecology Letters*, 10, 135--145. \doi{10.1111/j.1461-0248.2006.01006.x}
#' @seealso [trait_gradient()], [trait_ou()] to simulate traits under these models.
#' @family trait variation among individuals
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(20)
#' sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree, lambda = 1),
#'                               competition_matrix(tree, "none"), n_steps = 100,
#'                               seed_steps = 1:5)
#' ind <- sample_individuals(sim, n = 5, local_adaptation = 0.5, intraspecific_sd = 1)
#' tg <- trait_gradient(ind, abundance = "abundance")
#' trait_gradient_phylo(tg, tree, n_boot = 99)
#' @export
trait_gradient_phylo <- function(x, tree, responses = c("beta", "alpha", "mean"),
                                 models = c("white", "BM", "OU"), n_boot = 199,
                                 min_sites = 3) {
  if (!inherits(x, "mc_trait_gradient")) {
    stop("`x` must come from trait_gradient()", call. = FALSE)
  }
  check_tree(tree)
  responses <- unique(match.arg(responses, several.ok = TRUE))
  models <- unique(match.arg(models, several.ok = TRUE))
  check_count(n_boot, "n_boot", min = 2)
  check_count(min_sites, "min_sites", min = 2)
  boot <- tg_boot_se(x, n_boot)
  s <- x$species
  used <- (s$n_sites - s$n_alone) >= min_sites & s$species %in% tree$tip.label &
    stats::complete.cases(s$beta, s$alpha)
  if (sum(used) < 4) {
    stop("fewer than 4 species have ", min_sites, " sites or more and are in `tree`",
         call. = FALSE)
  }
  data <- data.frame(species = s$species[used], n_sites = s$n_sites[used] - s$n_alone[used],
                     mean = s$mean[used], se_mean = boot[used, "mean"],
                     beta = s$beta[used], se_beta = boot[used, "beta"],
                     alpha = s$alpha[used], se_alpha = boot[used, "alpha"],
                     row.names = NULL)
  pruned <- ape::keep.tip(tree, data$species)
  covs <- tg_tree_parts(pruned, data$species)
  fits <- lapply(responses, function(r) {
    y <- data[[r]]
    v <- data[[paste0("se_", r)]]^2
    do.call(rbind, lapply(models, function(m) {
      f <- gls_known_error(y, v, covs, m)
      data.frame(response = r, model = m, f, row.names = NULL)
    }))
  })
  table <- do.call(rbind, fits)
  table$delta_aic <- table$aic - stats::ave(table$aic, table$response, FUN = min)
  structure(list(models = table, data = data, tree = pruned,
                 settings = c(x$settings, list(n_boot = n_boot, min_sites = min_sites))),
            class = "mc_trait_gradient_phylo")
}

# bootstrap standard errors of mean, beta and alpha of each species, resampling its sites with
# the site means of the observed data held fixed
tg_boot_se <- function(x, n_boot) {
  p <- x$populations
  rows <- split(seq_len(nrow(p)), factor(p$species, x$species$species))
  out <- t(vapply(rows, function(r) {
    ok <- r[!is.na(p$site_mean[r])]
    if (length(ok) < 2) return(c(NA_real_, NA_real_, NA_real_))
    draws <- replicate(n_boot, {
      pick <- r[sample.int(length(r), length(r), replace = TRUE)]
      w <- p$weight[pick]
      tbar <- sum(w * p$trait[pick]) / sum(w)
      keep <- !is.na(p$site_mean[pick])
      beta <- if (any(keep)) sum(w[keep] * p$site_mean[pick][keep]) / sum(w[keep]) else NA
      c(tbar, beta, tbar - beta)
    })
    apply(draws, 1, stats::sd, na.rm = TRUE)
  }, numeric(3)))
  colnames(out) <- c("mean", "beta", "alpha")
  out
}

# shared lengths and root-to-tip distances of the tips, in the order of `species`
tg_tree_parts <- function(tree, species) {
  shared <- ape::vcv(tree)[species, species]
  list(shared = shared, depth = diag(shared), height = max(diag(shared)))
}

# REML fit of y ~ N(mu, sigma2 C + diag(v)) for a covariance model
gls_known_error <- function(y, v, covs, model) {
  n <- length(y)
  ones <- rep(1, n)
  cov_of <- function(a) {
    switch(model,
           white = diag(n),
           BM = covs$shared,
           OU = {
             d <- covs$depth
             s <- covs$shared
             exp(-a * (outer(d, d, "+") - 2 * s)) * -expm1(-2 * a * s) / (2 * a)
           })
  }
  reml <- function(log_sigma2, a = NULL) {
    vmat <- exp(log_sigma2) * cov_of(a) + diag(v, n)
    ch <- tryCatch(chol(vmat), error = function(e) NULL)
    if (is.null(ch)) return(-1e10) # nocov: a singular covariance (tips at distance 0)
    vi_1 <- backsolve(ch, forwardsolve(t(ch), ones))
    vi_y <- backsolve(ch, forwardsolve(t(ch), y))
    xvx <- sum(ones * vi_1)
    mu <- sum(ones * vi_y) / xvx
    r <- y - mu
    vi_r <- vi_y - mu * vi_1
    ll <- -0.5 * (2 * sum(log(diag(ch))) + log(xvx) + sum(r * vi_r) + (n - 1) * log(2 * pi))
    attr(ll, "mu") <- mu
    ll
  }
  scale <- stats::var(y) / if (model == "white") 1 else covs$height
  lower <- log(scale) - 12
  upper <- log(scale) + 8
  if (model == "OU") {
    h <- covs$height
    objective <- function(par) -reml(par[1], exp(par[2]) / h)
    starts <- expand.grid(log(scale), log(c(0.1, 1, 5)))
    fits <- lapply(seq_len(nrow(starts)), function(i) {
      stats::optim(unlist(starts[i, ]), objective, method = "L-BFGS-B",
                   lower = c(lower, log(1e-3)), upper = c(upper, log(50)))
    })
    best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
    a <- exp(best$par[2]) / h
    log_s2 <- best$par[1]
  } else {
    best <- stats::optimize(function(z) -reml(z), c(lower, upper))
    log_s2 <- best$minimum
    a <- NA_real_
  }
  ll <- reml(log_s2, a)
  k <- 2 + (model == "OU")
  data.frame(mu = attr(ll, "mu"), sigma2 = exp(log_s2), a = a,
             half_life = if (model == "OU") log(2) / a else NA_real_,
             loglik = as.numeric(ll), n_par = k, aic = -2 * as.numeric(ll) + 2 * k)
}

#' @export
print.mc_trait_gradient_phylo <- function(x, ...) {
  cat("Phylogenetic models of the trait-gradient analysis:", nrow(x$data), "species\n")
  m <- x$models
  num <- vapply(m, is.numeric, logical(1))
  m[num] <- lapply(m[num], function(v) signif(v, 3))
  print(m, row.names = FALSE)
  invisible(x)
}

#' Plot a trait-gradient analysis
#'
#' `type = "gradient"` draws the trait value of each population against the mean of its site,
#' with the line of each species over the range of its sites (its slope \eqn{b_i}) and the
#' species means at \eqn{(\beta_i, \bar t_i)}; the dashed line is \eqn{y = x}, where
#' \eqn{\alpha = 0}. `type = "beta_alpha"` draws \eqn{\alpha} against \eqn{\beta}. For
#' [trait_gradient_test()], the null distributions of the three statistics with the observed
#' value.
#'
#' @param x Result of [trait_gradient()] or [trait_gradient_test()].
#' @param type `"gradient"` or `"beta_alpha"`.
#' @param species_labels Write the names of species in `type = "beta_alpha"`.
#' @param colour Colour of the species lines and means.
#' @param point_size Size of the points.
#' @param ... Not used.
#' @return A ggplot object.
#' @seealso [trait_gradient()], [trait_gradient_test()]
#' @family trait variation among individuals
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
#' set.seed(1)
#' tree <- treeshapesim::sim_tree_shape(20)
#' sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
#'                               competition_matrix(tree, "none"), n_steps = 100,
#'                               seed_steps = 1:5)
#' ind <- sample_individuals(sim, n = 5, local_adaptation = 0.5, intraspecific_sd = 1)
#' tg <- trait_gradient(ind, abundance = "abundance")
#' plot(tg)
#' plot(tg, type = "beta_alpha", species_labels = TRUE)
#' plot(trait_gradient_test(tg, n_null = 99))
#' @export
plot.mc_trait_gradient <- function(x, type = c("gradient", "beta_alpha"),
                                   species_labels = FALSE, colour = "#0072B2",
                                   point_size = 1.5, ...) {
  need_ggplot("plot() of trait_gradient()")
  .data <- ggplot2::.data
  type <- match.arg(type)
  check_colours(colour, 1)
  s <- x$species
  if (type == "beta_alpha") {
    p <- ggplot2::ggplot(s, ggplot2::aes(x = .data$beta, y = .data$alpha)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::geom_point(colour = colour, size = 2 * point_size) +
      ggplot2::labs(x = expression(beta ~ "(position along the gradient)"),
                    y = expression(alpha ~ "(difference from co-occurring species)")) +
      ggplot2::theme_minimal()
    if (species_labels) {
      p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$species), vjust = -0.8, size = 3)
    }
    return(p)
  }
  pops <- x$populations[!is.na(x$populations$site_mean) & x$populations$measured, ]
  lines <- s[!is.na(s$slope), ]
  lines$y_min <- lines$intercept + lines$slope * lines$min_site
  lines$y_max <- lines$intercept + lines$slope * lines$max_site
  ggplot2::ggplot() +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(data = pops, ggplot2::aes(x = .data$site_mean, y = .data$trait),
                        colour = "grey60", size = point_size) +
    ggplot2::geom_segment(data = lines, ggplot2::aes(x = .data$min_site, xend = .data$max_site,
                                                     y = .data$y_min, yend = .data$y_max),
                          colour = colour, linewidth = 0.6) +
    ggplot2::geom_point(data = s, ggplot2::aes(x = .data$beta, y = .data$mean),
                        colour = colour, size = 2 * point_size) +
    ggplot2::labs(x = if (x$settings$leave_one_out) "Site mean (without the species)" else
      "Site mean", y = paste("Trait value:", x$settings$trait)) +
    ggplot2::theme_minimal()
}

#' @rdname plot.mc_trait_gradient
#' @export
plot.mc_trait_gradient_test <- function(x, ...) {
  need_ggplot("plot() of trait_gradient_test()")
  .data <- ggplot2::.data
  labels <- c(cor_beta_alpha = "Correlation of beta and alpha", var_beta = "Variance of beta",
              var_alpha = "Variance of alpha")
  null <- data.frame(statistic = rep(labels, each = nrow(x$null)),
                     value = as.vector(x$null[, names(labels)]))
  observed <- data.frame(statistic = labels, value = x$statistics$observed)
  null$statistic <- factor(null$statistic, labels)
  observed$statistic <- factor(observed$statistic, labels)
  ggplot2::ggplot(null, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = 30, fill = "grey70", colour = "white") +
    ggplot2::geom_vline(data = observed, ggplot2::aes(xintercept = .data$value),
                        colour = "#D55E00", linewidth = 0.9) +
    ggplot2::facet_wrap(~statistic, scales = "free") +
    ggplot2::labs(x = NULL, y = "Null communities") +
    ggplot2::theme_minimal()
}
