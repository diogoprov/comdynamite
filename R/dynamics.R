#' Simulate a metacommunity
#'
#' Runs the process-based model of Thompson et al. (2020) on a landscape. At
#' each time step, abundances are drawn as
#' \eqn{N' \sim \mathrm{Poisson}(r N / (1 + N \alpha'))}{N' ~ Poisson(r N / (1 + N alpha'))},
#' where \eqn{r} holds the growth rates of each species at each site ([growth_rates()]) and
#' \eqn{\alpha}{alpha} the competition coefficients ([competition_matrix()]). Each individual
#' then emigrates with
#' probability `dispersal`, and emigrants are redistributed among sites:
#'
#' * spatially implicit landscapes ([landscape_implicit()]): the emigrants of
#'   each species are distributed multinomially, each site receiving them with
#'   probability proportional to the number of emigrants from the other sites
#'   (a uniform kernel; an individual may return to its own site);
#' * spatially explicit landscapes ([landscape_grid()]): the number of
#'   immigrants of each species at each site is drawn from a Poisson
#'   distribution whose mean is the expected flow of emigrants through the
#'   dispersal kernel. This approximation does not conserve the number of
#'   dispersing individuals exactly, but has the right expectation. The code of
#'   Thompson et al. (2020) does the same. Few emigrants may then produce no
#'   immigrant, so with `extirpation` and low dispersal, recolonisation is rarer
#'   than in mcomsimr, which conserves dispersers. With `immigration = "multinomial"`,
#'   the emigrants of each species at each site are instead distributed among the other
#'   sites multinomially, with the probabilities of the kernel, which conserves them. It
#'   costs about as much as the Poisson draw at low dispersal and up to 20 times more at high
#'   dispersal, where the two give nearly the same results.
#'
#' After dispersal, each population is lost with probability `extirpation`,
#' independently of its size. Thompson et al. (2020) used 0.002 in their
#' competition-colonisation scenario and no extirpation in the others.
#'
#' During colonisation, Poisson(`seed_mean`) individuals of every species are
#' added to every site at each of the `seed_steps`. After the last of these,
#' species absent from all sites can never return and are no longer computed.
#'
#' Distributions from elsewhere, such as ranges, biogeographic reconstructions or the output of
#' gen3sis, enter as `pool` and `initial`. With `pool`, a species can only occur where its value
#' is `TRUE` (or 1): seeds go only there, and individuals that disperse to other sites are lost,
#' as in a hostile matrix. A species without any site in `pool` dies out. `initial` replaces the
#' first draw of Poisson(`seed_mean`) individuals with given abundances, or with one individual
#' where a matrix of presences has 1; seeding still follows `seed_steps`, so use
#' `seed_steps = numeric(0)` to start from `initial` alone. Both take one row per site and one
#' column per species, named by species. Rows are matched by name when they have names (`site1`,
#' `site2` and so on in [simulate_metacommunity()]) and by position otherwise. Without them,
#' results with the same seed are the same as before they existed. A tree from elsewhere, such as
#' one simulated with gen3sis or a model of state-dependent diversification, enters through
#' [species_traits()]; drop extinct lineages first, since they have no distribution.
#'
#' @section Article and archived code:
#' The article of Thompson et al. (2020) and the code archived with it
#' (\doi{10.5281/zenodo.3708350}) differ in points that change results. Beta is
#' gamma / alpha in the code and in their figures, but gamma - alpha in their
#' Table 1. Their environment is converted to ranks, so it is uniform between 0
#' and 1, and the lowest dispersal is 0.00001, not 0.0001. Emigrants are Poisson in
#' the code and binomial in the article, which comdynamite follows. In their
#' competition-colonisation scenario, the code also halves the immigrants. See
#' the article "Validating the model of Thompson et al. (2020)" on the package
#' website for the complete list.
#'
#' @param landscape An `mc_landscape` object.
#' @param traits Species traits from [species_traits()], or, for a niche with several
#'   environmental axes, a list of them, one per column of `env`, with the same species; see
#'   [growth_rates()].
#' @param env Optional environmental axes of each site, a matrix or data frame with one row per
#'   site and one column per element of `traits` (such as the `env` of several landscapes on
#'   the same grid). Needed when `traits` is a list; replaces the `env` of `landscape`. Not
#'   available with `env_time`.
#' @param alpha Competition matrix from [competition_matrix()]. Species are matched to the
#'   rows of `traits` by name, or by position if `alpha` has no names.
#' @param dispersal Emigration probability, one value or one per species.
#' @param extirpation Probability that a population is lost at each time step,
#'   one value or one per species.
#' @param n_steps Number of time steps.
#' @param env_time Optional matrix (sites x `n_steps`) with the environment at
#'   each step, e.g. from [env_over_time()]. If `NULL`, the environment of the
#'   landscape is constant.
#' @param seed_steps Steps at which individuals are added during colonisation.
#' @param seed_mean Mean number of individuals of each species added to each
#'   site at each seeding step (and at the start).
#' @param record_steps Optional steps at which abundances are recorded.
#' @param immigration How emigrants reach other sites on spatially explicit landscapes:
#'   `"poisson"` (the default, as in the code of Thompson et al. 2020) or `"multinomial"`,
#'   which conserves them. Spatially implicit landscapes are always multinomial.
#' @param pool Optional matrix (sites x species), logical or 0/1, with the sites where each
#'   species can occur.
#' @param initial Optional matrix (sites x species) of abundances or 0/1 presences at the start,
#'   only where `pool` allows.
#' @return An `mc_sim` object: a list with `abundance` (sites x species matrix
#'   at the last step), `history` (sites x species x recorded steps array, or
#'   `NULL`), and the `landscape`, `traits`, `alpha`, `dispersal`, `extirpation` and `pool`
#'   used, `env` when given, and `tree` when the traits came from [species_traits()] with a
#'   tree.
#' @references Thompson, P. L., Guzman, L. M., De Meester, L., Horváth, Z.,
#'   Ptacnik, R., Vanschoenwinkel, B., Viana, D. S., & Chase, J. M. (2020). A
#'   process-based metacommunity framework linking local and regional scale
#'   community ecology. *Ecology Letters*, 23, 1314--1329.
#'   \doi{10.1111/ele.13568}
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' ls <- landscape_implicit(30)
#' tr <- treeshapesim::sim_tree_shape(20)
#' sim <- simulate_metacommunity(ls, species_traits(tr),
#'                               competition_matrix(tr, "phylo"),
#'                               dispersal = 0.01, n_steps = 200)
#' sim
#' # a niche with two environmental axes
#' ls2 <- landscape_implicit(30)
#' niche2 <- species_traits(tr, rmax_range = c(1, 1))
#' sim2 <- simulate_metacommunity(ls, list(species_traits(tr), niche2),
#'                                competition_matrix(tr, "phylo"), n_steps = 200,
#'                                env = cbind(ls$env, ls2$env))
#' @export
simulate_metacommunity <- function(landscape, traits, alpha, dispersal = 0.01,
                                   n_steps = 1000, env_time = NULL,
                                   seed_steps = seq(10, 100, 10), seed_mean = 0.5,
                                   record_steps = NULL, extirpation = 0, pool = NULL,
                                   initial = NULL, immigration = c("poisson", "multinomial"),
                                   env = NULL) {
  stopifnot(inherits(landscape, "mc_landscape"))
  immigration <- match.arg(immigration)
  several <- is.list(traits) && !is.data.frame(traits)
  if (several || !is.null(env)) {
    if (!is.null(env_time)) stop("`env_time` is not available with `env`", call. = FALSE)
    trait_list <- if (several) lapply(traits, as_traits) else list(as_traits(traits))
    env <- as.data.frame(env %||% stop("a list of traits needs `env`, one column per element",
                                       call. = FALSE))
    if (nrow(env) != length(landscape$env) || ncol(env) != length(trait_list)) {
      stop("`env` needs one row per site (", length(landscape$env), ") and one column per ",
           "set of traits (", length(trait_list), ")", call. = FALSE)
    }
    traits <- trait_list[[1]]
  }
  if (!identical(dim(alpha), rep(nrow(traits), 2L))) {
    stop("`alpha` needs one row and one column per species (", nrow(traits), ")", call. = FALSE)
  }
  species <- rownames(traits)
  if (!is.null(rownames(alpha)) && !is.null(colnames(alpha))) {
    if (!setequal(rownames(alpha), species) || !setequal(colnames(alpha), species)) {
      stop("names of species in `alpha` and `traits` differ", call. = FALSE)
    }
    alpha <- alpha[species, species]
  }
  sites <- paste0("site", seq_along(landscape$env))
  r <- if (!is.null(env)) {
    rownames(env) <- sites
    growth_rates(env, if (several) trait_list else traits)
  } else if (is.null(env_time)) {
    growth_rates(stats::setNames(landscape$env, sites), traits)
  } else {
    stopifnot(nrow(env_time) == length(landscape$env), ncol(env_time) >= n_steps)
    function(step) growth_rates(stats::setNames(env_time[, step], sites), traits)
  }
  out <- mc_dynamics(r, alpha, dispersal, n_steps = n_steps, seed_steps = seed_steps,
                     seed_mean = seed_mean, disp_matrix = landscape$dispersal,
                     record_steps = record_steps, sites = sites, species = rownames(traits),
                     extirpation = extirpation, pool = pool, initial = initial,
                     immigration = immigration)
  if (several) traits <- trait_list
  structure(c(out, list(landscape = landscape, env = env, traits = traits, alpha = alpha,
                        dispersal = dispersal, extirpation = extirpation,
                        pool = out$pool, tree = attr(traits, "tree"))), class = "mc_sim")
}

#' Metacommunity dynamics from growth rates
#'
#' The engine of [simulate_metacommunity()], for users who want to supply
#' their own growth rates (e.g. responses to several environmental variables).
#'
#' @inheritParams simulate_metacommunity
#' @param r Sites x species matrix of density-independent growth rates, or a
#'   function of the time step that returns such a matrix.
#' @param disp_matrix Row-stochastic dispersal matrix (sites x sites), or
#'   `NULL` for uniform dispersal.
#' @param sites,species Optional names of sites and species.
#' @return List with `abundance` and `history`, as in
#'   [simulate_metacommunity()], and `pool` as a logical matrix, or `NULL` without it.
#' @examples
#' set.seed(1)
#' # one species alone converges to (r - 1) / alpha = (5 - 1) / 0.05 = 80
#' out <- mc_dynamics(matrix(5, 20, 1), matrix(0.05), dispersal = 0, n_steps = 200)
#' mean(out$abundance)
#' @export
mc_dynamics <- function(r, alpha, dispersal, n_steps = 1000,
                        seed_steps = seq(10, 100, 10), seed_mean = 0.5,
                        disp_matrix = NULL, record_steps = NULL,
                        sites = NULL, species = NULL, extirpation = 0, pool = NULL,
                        initial = NULL, immigration = c("poisson", "multinomial")) {
  immigration <- match.arg(immigration)
  r_fun <- if (is.function(r)) r else function(step) r
  r1 <- r_fun(1)
  m <- nrow(r1)
  k_all <- ncol(r1)
  stopifnot(all(dim(alpha) == k_all), n_steps >= 1)
  if (!is.null(disp_matrix)) stopifnot(all(dim(disp_matrix) == m))
  dispersal <- rep_len(dispersal, k_all)
  extirpation <- rep_len(extirpation, k_all)
  stopifnot(all(dispersal >= 0 & dispersal <= 1), all(extirpation >= 0 & extirpation <= 1))
  if (m == 1 && any(dispersal > 0)) {
    stop("dispersal needs at least two sites: use `dispersal = 0` with one site", call. = FALSE)
  }
  sites <- sites %||% rownames(r1) %||% paste0("site", seq_len(m))
  species <- species %||% colnames(r1) %||% paste0("sp", seq_len(k_all))
  if (!is.null(pool)) {
    pool <- site_species(pool, sites, species, "`pool`")
    if (!all(pool %in% c(0, 1))) stop("`pool` must be logical or 0/1", call. = FALSE)
    pool <- pool == 1
  }
  if (!is.null(initial)) {
    initial <- site_species(initial, sites, species, "`initial`")
    if (any(initial < 0) || any(initial != round(initial))) {
      stop("`initial` must hold whole numbers of individuals, zero or more", call. = FALSE)
    }
    if (!is.null(pool) && any(initial[!pool] > 0)) {
      stop("`initial` has individuals outside `pool`", call. = FALSE)
    }
    storage.mode(initial) <- "integer"
  }
  record_steps <- sort(unique(record_steps[record_steps >= 1 & record_steps <= n_steps]))
  history <- if (length(record_steps))
    array(0L, c(m, k_all, length(record_steps)),
          dimnames = list(sites, species, paste0("step", record_steps))) else NULL

  alive <- seq_len(k_all)
  k <- k_all
  p_em <- rep(dispersal, each = m)
  p_ex <- rep(extirpation, each = m)
  at <- t(alpha)
  inside <- pool
  n <- if (is.null(initial)) matrix(stats::rpois(m * k, seed_mean), m, k) else unname(initial)
  if (!is.null(pool)) n[!pool] <- 0L
  last_seed <- max(c(0, seed_steps))
  for (step in seq_len(n_steps)) {
    if (step %in% seed_steps) {
      seeds <- stats::rpois(m * k, seed_mean)
      n <- n + if (is.null(inside)) seeds else seeds * inside
    }
    keep <- colSums(n) > 0
    if (step > last_seed && !all(keep)) {
      alive <- alive[keep]
      n <- n[, keep, drop = FALSE]
      at <- at[keep, keep, drop = FALSE]
      k <- length(alive)
      p_em <- rep(dispersal[alive], each = m)
      p_ex <- rep(extirpation[alive], each = m)
      if (!is.null(inside)) inside <- pool[, alive, drop = FALSE]
    }
    if (k > 0) {
      r_now <- r_fun(step)[, alive, drop = FALSE]
      lam <- r_now * n / (1 + n %*% at)
      n <- matrix(stats::rpois(m * k, lam), m, k)
      em <- matrix(stats::rbinom(m * k, n, p_em), m, k)
      tot <- colSums(em)
      if (any(tot > 0)) {
        if (!is.null(disp_matrix) && immigration == "multinomial") {
          n <- n - em + place_emigrants(em, disp_matrix)
        } else if (!is.null(disp_matrix)) {
          n <- n - em + matrix(stats::rpois(m * k, crossprod(disp_matrix, em)), m, k)
        } else {
          w <- sweep(-em, 2, tot, "+")
          imm <- matrix(0L, m, k)
          for (j in which(tot > 0)) imm[, j] <- stats::rmultinom(1, tot[j], w[, j])
          n <- n - em + imm
        }
      }
      # no draws without extirpation, so that seeds give the same results as before
      if (any(p_ex > 0)) n[stats::runif(m * k) < p_ex] <- 0L
      if (!is.null(inside)) n[!inside] <- 0L
    }
    if (!is.null(history) && step %in% record_steps)
      history[, alive, match(step, record_steps)] <- n
  }
  abundance <- matrix(0L, m, k_all, dimnames = list(sites, species))
  abundance[, alive] <- n
  list(abundance = abundance, history = history, pool = pool)
}

# immigrants when the emigrants of each species at each site go to the other sites
# multinomially, with the probabilities of the rows of `disp_matrix`
place_emigrants <- function(em, disp_matrix) {
  m <- nrow(em)
  immigrants <- matrix(0L, m, ncol(em))
  for (cell in which(em > 0)) {
    from <- (cell - 1L) %% m + 1L
    species <- (cell - 1L) %/% m + 1L
    immigrants[, species] <- immigrants[, species] +
      drop(stats::rmultinom(1, em[cell], disp_matrix[from, ]))
  }
  immigrants
}

# a sites x species matrix in the order of `sites` and `species`: rows by name when they have
# names, by position otherwise; columns by species name
site_species <- function(x, sites, species, name) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !(is.numeric(x) || is.logical(x))) {
    stop(name, " must be a matrix of sites x species", call. = FALSE)
  }
  if (is.null(colnames(x))) stop(name, " needs species names in its columns", call. = FALSE)
  missing <- setdiff(species, colnames(x))
  if (length(missing)) {
    stop(name, " has no column for the species ", species_list(missing), call. = FALSE)
  }
  extra <- setdiff(colnames(x), species)
  if (length(extra)) {
    warning(name, " has species that are not simulated, dropped: ", species_list(extra),
            call. = FALSE)
  }
  rows <- if (is.null(rownames(x))) {
    if (nrow(x) != length(sites)) {
      stop(name, " has ", nrow(x), " rows and the landscape ", length(sites), " sites",
           call. = FALSE)
    }
    seq_along(sites)
  } else {
    match(sites, rownames(x))
  }
  if (anyNA(rows)) {
    stop(name, " has no row for the sites ", species_list(sites[is.na(rows)]), call. = FALSE)
  }
  if (anyNA(x)) stop(name, " has missing values", call. = FALSE)
  out <- x[rows, species, drop = FALSE]
  dimnames(out) <- list(sites, species)
  out
}

#' @export
print.mc_sim <- function(x, ...) {
  a <- x$abundance
  cat("<mc_sim>", nrow(a), "sites x", ncol(a), "species\n")
  cat("species persisting:", sum(colSums(a) > 0), "of", ncol(a), "\n")
  cat("mean local richness:", format(mean(rowSums(a > 0)), digits = 3), "\n")
  if (!is.null(x$history)) cat("recorded steps:", dim(x$history)[3], "\n")
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
