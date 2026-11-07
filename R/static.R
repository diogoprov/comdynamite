#' Simulate a static metacommunity
#'
#' Draws species abundances directly from unimodal responses to an
#' environmental gradient, without population dynamics, competition or
#' dispersal. The expected abundance of species j at site i is
#' `r_max[j] * exp(-(env[i] - optimum[j])^2 / (2 * sigma[j]^2))`, the Gaussian
#' response of Dray & Legendre (2008), also used by Peres-Neto et al. (2012). By
#' default, observed abundances are Poisson counts with that mean, as in
#' Peres-Neto & Kembel (2015) and Duarte et al. (2016), where `r_max` lies between
#' 0 and 30. In the two earlier studies `r_max` lies between 0.5 and 1, abundances
#' are the expected values (`sampling = "expected"`), and Dray & Legendre (2008)
#' draw presences from a binomial distribution with the expected value as the
#' probability (`sampling = "binomial"`). Because nothing but the
#' gradient structures the data, this is a baseline for
#' [simulate_metacommunity()] and a fast generator for evaluating methods.
#'
#' With several gradients, for example an environmental and a spatial one, each
#' species has its own response to each gradient and the expected abundance is
#' the product of the responses, as in Peres-Neto et al. (2012).
#'
#' @param traits Data frame from [species_traits()], or any data frame with
#'   columns `optimum`, `sigma` and `r_max`, or a list of them, one per
#'   gradient, with the same species. Here `r_max` is the maximum expected
#'   abundance, not a growth rate, so set `rmax_range` accordingly. With several
#'   gradients the maxima multiply, so a range such as `c(1, 1)` in all but one
#'   of them keeps abundances on the scale of the first.
#' @param env Value of the gradient at each site, or a list or data frame with
#'   one gradient per element, in the order of `traits`.
#' @param incidence Return presences and absences instead of Poisson counts? A
#'   species is present where its count is positive. Only for
#'   `sampling = "poisson"`.
#' @param sampling How to obtain the data from the expected values: `"poisson"`
#'   draws counts, `"binomial"` draws presences and absences, which needs
#'   expected values of at most 1, and `"expected"` returns the expected values.
#' @return List of class `mc_sim` with `abundance`, a matrix with one row per
#'   site and one column per species, plus `env`, `traits` and, when the traits came from
#'   [species_traits()] with a tree, `tree`. Sites and species with no individuals are kept,
#'   so remove them if a method requires it.
#' @references
#' Dray, S., & Legendre, P. (2008). Testing the species traits-environment
#' relationships: the fourth-corner problem revisited. *Ecology*, 89,
#' 3400--3412. \doi{10.1890/08-0349.1}
#'
#' Duarte, L. D. S., et al. (2016). Dissecting phylogenetic fuzzy weighting:
#' theory and application in metacommunity phylogenetics. *Methods in Ecology
#' and Evolution*, 7, 937--946. \doi{10.1111/2041-210X.12547}
#'
#' Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient analysis:
#' environmental drivers of phylogenetic variation across ecological
#' communities. *Plant Ecology*, 216, 709--724. \doi{10.1007/s11258-014-0405-0}
#'
#' Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects
#' of spatial contingency and environmental filtering on metacommunity
#' phylogenetics. *Ecology*, 93, S14--S30. \doi{10.1890/11-0494.1}
#' @examplesIf requireNamespace("treeshapesim", quietly = TRUE)
#' set.seed(1)
#' tr <- treeshapesim::sim_tree_shape(20)
#' traits <- species_traits(tr, lambda = 1, optima_range = c(-1, 101),
#'                          rmax_range = c(0, 30))
#' sim <- simulate_static(traits, env = runif(100, 0, 100))
#' sim
#'
#' # an environmental and a spatial gradient
#' space <- species_traits(tr, lambda = 1, optima_range = c(-1, 101),
#'                         rmax_range = c(1, 1))
#' gradients <- data.frame(env = runif(100, 0, 100), space = 1:100)
#' simulate_static(list(traits, space), gradients)
#' @export
simulate_static <- function(traits, env, incidence = FALSE,
                            sampling = c("poisson", "binomial", "expected")) {
  sampling <- match.arg(sampling)
  trait_list <- if (is.list(traits) && !is.data.frame(traits)) traits else list(traits)
  trait_list <- lapply(trait_list, as_traits)
  env_list <- if (is.list(env)) as.list(env) else list(env)
  stopifnot(sampling == "poisson" || !incidence, length(trait_list) == length(env_list),
            all(vapply(env_list, is.numeric, logical(1))),
            length(unique(lengths(env_list))) == 1, length(env_list[[1]]) >= 1,
            all(vapply(trait_list, function(x) all(x$r_max >= 0), logical(1))))
  species <- rownames(trait_list[[1]])
  if (!all(vapply(trait_list, function(x) setequal(rownames(x), species), logical(1)))) {
    stop("sets of traits have different species", call. = FALSE)
  }
  trait_list <- lapply(trait_list, function(x) x[species, , drop = FALSE])
  expected <- Reduce(`*`, Map(growth_rates, unname(env_list), trait_list))
  if (sampling == "binomial" && any(expected > 1)) {
    stop("binomial sampling needs expected values of at most 1, so set r_max accordingly",
         call. = FALSE)
  }
  n <- length(expected)
  drawn <- switch(sampling, poisson = stats::rpois(n, expected),
                  binomial = stats::rbinom(n, 1, expected), expected = expected)
  abundance <- matrix(drawn, nrow(expected), ncol(expected), dimnames = list(
                        rownames(expected) %||% paste0("site", seq_along(env_list[[1]])),
                        rownames(trait_list[[1]])))
  if (incidence) abundance[] <- (abundance > 0) * 1L
  structure(list(abundance = abundance, history = NULL, env = env, traits = traits,
                 tree = attr(traits, "tree")), class = "mc_sim")
}

#' Simulate presences and absences from logistic responses
#'
#' The two-gradient simulation of Appendix C of Peres-Neto et al. (2012), after
#' Ives & Helmus (2011). The probability that species j occurs at site t along
#' gradient x is `plogis(a[j] + b[j] * x[t])`, where `a` is a normal deviate and
#' the slope `b` evolves on the tree. With several gradients, each with its own
#' `a` and `b`, the probabilities multiply. A value drawn for each site from
#' `reduction` then reduces them, which adds noise.
#'
#' Appendix C says that a species is present when a normally generated value is
#' smaller than its probability. Taken literally, a species with a probability
#' near zero is present half of the time, so the default is a uniform value
#' (a Bernoulli draw), and `draw = "normal"` follows the text.
#'
#' @param tree A `phylo` tree with branch lengths. Appendix C uses a fully
#'   balanced tree with 32 species.
#' @param gradients Value of the gradient at each site, or a list or data frame
#'   with one gradient per element, such as [logistic_gradients()]. Gradients are
#'   centred.
#' @param intercept_mean,intercept_sd Mean and standard deviation of `a`. Appendix C does
#'   not give these values, and the defaults are those of the code of Ives & Helmus (2011).
#' @param slope_sd Standard deviation of the slopes `b` among the species of the tree, as in
#'   the code of Ives & Helmus (2011). Appendix C does not give this value.
#' @param lambda Pagel's lambda of the slopes.
#' @param reduction Range of the uniform distribution of the site reduction.
#' @param draw Distribution of the value compared with each probability.
#' @return List of class `mc_sim` with `abundance`, a 0/1 matrix with one row per
#'   site and one column per species; `probability`, the matrix of probabilities;
#'   `env`, the gradients; and `traits`, a list with the matrices `intercept`
#'   and `slope`, with one row per species and one column per gradient; and `tree`.
#' @references
#' Ives, A. R., & Helmus, M. R. (2011). Generalized linear mixed models for
#' phylogenetic analyses of community structure. *Ecological Monographs*, 81,
#' 511--525. \doi{10.1890/10-1264.1}
#'
#' Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects
#' of spatial contingency and environmental filtering on metacommunity
#' phylogenetics. *Ecology*, 93, S14--S30. \doi{10.1890/11-0494.1}
#' @examples
#' set.seed(1)
#' tr <- ape::compute.brlen(ape::stree(32, "balanced"), method = "Grafen")
#' simulate_logistic(tr, logistic_gradients(1)[c("x", "y")])
#' @export
simulate_logistic <- function(tree, gradients, intercept_mean = 0, intercept_sd = 0.5,
                              slope_sd = 1, lambda = 1, reduction = c(0.5, 1),
                              draw = c("uniform", "normal")) {
  draw <- match.arg(draw)
  stopifnot(inherits(tree, "phylo"), intercept_sd >= 0, slope_sd >= 0,
            length(reduction) == 2)
  grads <- if (is.list(gradients)) as.list(gradients) else list(gradients)
  stopifnot(all(vapply(grads, is.numeric, logical(1))),
            length(unique(lengths(grads))) == 1)
  sp <- tree$tip.label
  k <- length(sp)
  n <- length(grads[[1]])
  a <- vapply(grads, function(g) stats::rnorm(k, intercept_mean, intercept_sd), numeric(k))
  b <- vapply(grads, function(g) slope_sd * drop(scale(trait_bm(tree, lambda))), numeric(k))
  a <- matrix(a, k, dimnames = list(sp, names(grads)))
  b <- matrix(b, k, dimnames = list(sp, names(grads)))
  probability <- matrix(1, n, k, dimnames = list(paste0("site", seq_len(n)), sp))
  for (g in seq_along(grads)) {
    x <- grads[[g]] - mean(grads[[g]])
    probability <- probability * stats::plogis(sweep(outer(x, b[, g]), 2, a[, g], "+"))
  }
  probability <- probability * stats::runif(n, reduction[1], reduction[2])
  value <- if (draw == "uniform") stats::runif(n * k) else stats::rnorm(n * k)
  abundance <- (value < probability) * 1L
  structure(list(abundance = abundance, history = NULL, probability = probability,
                 env = gradients, traits = list(intercept = a, slope = b), tree = tree),
            class = "mc_sim")
}

#' Gradients of the seven scenarios of Peres-Neto et al. (2012)
#'
#' Two gradients on a transect, as in Appendix C of Peres-Neto et al. (2012).
#' Unless stated, `x` is evenly spaced on \[-15, 15\] and `y` is `x` in random
#' order, so that the two are independent.
#'
#' 1. The two gradients are equally important and independent.
#' 2. They are equally important and correlated: `y` is `x` plus a N(0, 1) error.
#' 3. They are independent, and `y` is less important: data come from `y`, and
#'    `y_test` is `y` with unit variance plus a N(0, 1) error. On the scale of `y` this
#'    error would have no effect, and with unit variance the results agree with Table 5.
#' 4. Neither is important: data come from scenario 1, and `x_test` and `y_test`
#'    are N(0, 1) deviates, to estimate the type I error.
#' 5. `y` has a broad spatial scale, the first Moran's eigenvector map (MEM), and `x` is
#'    in random order along the transect.
#' 6. `y` has an intermediate spatial scale, the third MEM, and `x` is as in scenario 5.
#' 7. `x` is the first MEM and `y` the third.
#'
#' MEMs are rescaled to \[-15, 15\].
#'
#' @param scenario Number of the scenario, from 1 to 7.
#' @param n_sites Number of sites.
#' @return Data frame with one row per site: `x` and `y`, the gradients that
#'   generate the data, and `x_test` and `y_test`, the ones to analyse.
#' @references
#' Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects
#' of spatial contingency and environmental filtering on metacommunity
#' phylogenetics. *Ecology*, 93, S14--S30. \doi{10.1890/11-0494.1}
#' @examples
#' set.seed(1)
#' head(logistic_gradients(6))
#' @export
logistic_gradients <- function(scenario, n_sites = 31) {
  stopifnot(length(scenario) == 1, scenario %in% 1:7)
  minimum <- if (scenario >= 5) 8 else 4
  if (n_sites < minimum) {
    stop("scenario ", scenario, " needs at least ", minimum, " sites", call. = FALSE)
  }
  x <- seq(-15, 15, length.out = n_sites)
  y <- sample(x)
  to_range <- function(v) rescale(v, c(-15, 15))
  mem <- if (scenario >= 5) transect_mem(n_sites)
  switch(as.character(scenario),
    "2" = y <- x + stats::rnorm(n_sites),
    "5" = {
      x <- sample(x)
      y <- to_range(mem[, 1])
    },
    "6" = {
      x <- sample(x)
      y <- to_range(mem[, 3])
    },
    "7" = {
      x <- to_range(mem[, 1])
      y <- to_range(mem[, 3])
    }
  )
  test <- switch(as.character(scenario),
    "3" = list(x, drop(scale(y)) + stats::rnorm(n_sites)),
    "4" = list(stats::rnorm(n_sites), stats::rnorm(n_sites)),
    list(x, y)
  )
  data.frame(x = x, y = y, x_test = test[[1]], y_test = test[[2]])
}

# Moran's eigenvector maps of a regular transect, from the broadest scale: the
# eigenvectors with positive eigenvalues of the doubly centred similarity
# matrix, with distances above one step truncated to four steps
transect_mem <- function(n) {
  d <- as.matrix(stats::dist(seq_len(n)))
  d[d > 1] <- 4
  s <- 1 - d^2 / max(d)^2
  diag(s) <- 0
  centring <- diag(n) - 1 / n
  e <- eigen(centring %*% s %*% centring, symmetric = TRUE)
  e$vectors[, e$values > 1e-8, drop = FALSE]
}
