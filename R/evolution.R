#' Evolution of local populations
#'
#' Settings for the evolution of the niche optimum within each population (a species at a site)
#' during [simulate_metacommunity()] or [mc_dynamics()]. Without them, every population of a
#' species keeps the optimum of the species. With them, each population has its own genetic
#' mean \eqn{\bar z}{zbar}, which responds to selection by the local environment, drifts, mixes
#' with immigrants and starts from a sample of its founders.
#'
#' At each step, in the order of the model:
#'
#' 1. **Seeds** carry the optimum of the species in `traits`, the mean of a regional pool that
#'    does not evolve. They mix with residents in proportion to their numbers; a population
#'    founded by seeds alone starts as below.
#' 2. **Phenotype and growth.** The phenotype is \eqn{p = \bar z + b (E - E_0)}{p = zbar + b
#'    (E - E0)}, a linear reaction norm with slope `plasticity` (b) around the reference
#'    environment \eqn{E_0}{E0}, the mean environment of the sites at the first step. The
#'    plastic part is immediate and reversible: every individual, migrants included, grows
#'    with the phenotype of the site where it is. Growth is the mean fitness of a population
#'    whose phenotypes are normal around \eqn{p}{p} with variance \eqn{P = G + V_E}{P = G +
#'    VE}, under the Gaussian response of [growth_rates()]: \eqn{r = r_{max} \sqrt{\sigma^2 /
#'    (\sigma^2 + P)} \exp(-(E - p)^2 / (2(\sigma^2 + P)))}{r = r_max sqrt(sigma^2 / (sigma^2
#'    + P)) exp(-(E - p)^2 / (2 (sigma^2 + P)))}. The variance widens the response and lowers its
#'    peak (the load); with \eqn{P = 0}{P = 0} it is the response of [growth_rates()].
#' 3. **Selection and drift**, in each population present after reproduction (\eqn{N > 0}):
#'    \deqn{\Delta \bar z = \frac{1}{g} \frac{G}{\sigma^2 + P} (E - p) + \epsilon, \quad
#'    \epsilon \sim N\left(0, \frac{G}{g N_e}\right), \quad N_e = \nu N,}{Delta zbar = (1/g)
#'    G / (sigma^2 + P) (E - p) + e, e ~ N(0, G / (g Ne)), Ne = nu N,}
#'    the response to Gaussian stabilising selection of Lande (1976), with additive genetic
#'    variance \eqn{G}{G} (`genetic_var`), phenotypic variance \eqn{P = G + V_E}{P = G + VE}
#'    (`env_var`), the niche breadth \eqn{\sigma}{sigma} of the species as the width of the
#'    fitness function, \eqn{g}{g} steps per generation (`steps_per_generation`) and the ratio
#'    \eqn{\nu}{nu} of effective to census size (`ne_ratio`), the same fitness function as in
#'    growth. \eqn{G}{G} is constant: it does not erode with drift nor grow when divergent
#'    populations mix.
#' 4. **Gene flow.** After dispersal, the genetic mean of a population is the mean of residents
#'    and immigrants weighted by their numbers. Immigrants carry the mean of the populations
#'    they come from, weighted by the emigrants each sends through the dispersal kernel (or, on
#'    spatially implicit landscapes, by the emigrants of all other sites). A population with no
#'    residents after reproduction and emigration (absent, or lost all its members) is founded
#'    by its \eqn{k}{k} immigrants: with `founder = TRUE`, their mean is drawn from
#'    \eqn{N(\bar z_{in}, (G + V_B)/k)}{N(zbar_in, (G + VB) / k)}, the sampling of a few
#'    individuals, where \eqn{V_B}{VB} is the variance among the means of their sources,
#'    weighted as above. Populations present at the start and those founded by seeds alone are
#'    drawn in the same way around the optimum of the species. Immigrants into a population
#'    with residents carry the expected mean of their sources, without sampling noise. The
#'    phenotype recorded after dispersal (`phenotype` in the output) mixes residents and
#'    immigrants still carrying the plastic response to their site of origin, as when
#'    individuals are sampled right after dispersal; it does not enter the dynamics.
#' 5. Populations lost to extirpation, competition or the `pool` lose their mean.
#'
#' Species with \eqn{G = 0}{G = 0} keep the optimum of the species exactly and draw no random
#' numbers, so with `genetic_var = 0`, `env_var = 0` and `plasticity = 0` the abundances are the
#' same as without evolution for the same seed. With [mc_dynamics()], growth rates come from
#' `env` and `traits` of `evolution`, and `r` only sets the numbers of sites and species.
#'
#' Two cases have known answers, checked in the tests of the package. A population alone,
#' without selection (very wide niche) and of constant size \eqn{N}{N}, accumulates a variance
#' of \eqn{G t / (\nu N)}{G t / (nu N)} after \eqn{t}{t} generations (Lande 1976). Two
#' populations of the same size in environments \eqn{E_1}{E1} and \eqn{E_2}{E2}, exchanging a
#' share \eqn{m}{m} of individuals each generation, with selection \eqn{a = G / (\sigma^2 +
#' P)}{a = G / (sigma^2 + P)} and no drift, reach the divergence \eqn{d^* = (1 - 2m) a (E_1 -
#' E_2) / (1 - (1 - 2m)(1 - a))}{d* = (1 - 2m) a (E1 - E2) / (1 - (1 - 2m)(1 - a))}, which
#' follows from the recursion of the steps above.
#'
#' @section Several traits:
#' With a niche of several environmental axes (a list of `traits` and an `env` with one column
#' per axis in [simulate_metacommunity()]), each population has a vector of genetic means, one
#' trait per axis: the optimum on that axis. `genetic_var` and `env_var` are then covariance
#' matrices \eqn{G}{G} and \eqn{E}{E}, with \eqn{P = G + E}{P = G + E}. The niche breadths
#' make the diagonal matrix \eqn{\Omega}{Omega} of Gaussian selection, and growth and the
#' response to selection are the multivariate versions of the steps above (Lande 1979):
#' \deqn{r = r_{max} \sqrt{\frac{|\Omega|}{|\Omega + P|}} \exp\left(-\frac{1}{2} d^\top
#' (\Omega + P)^{-1} d\right), \quad \Delta \bar z = \frac{1}{g} G (\Omega + P)^{-1} d +
#' \epsilon, \quad \epsilon \sim N\left(0, \frac{G}{g N_e}\right),}{r = r_max sqrt(|Omega| /
#' |Omega + P|) exp(-d' (Omega + P)^-1 d / 2), Delta zbar = G (Omega + P)^-1 d / g + e, e ~
#' N(0, G / (g Ne)),}
#' with \eqn{d}{d} the difference between the environment and the phenotype on each axis and
#' \eqn{r_{max}}{r_max} the product over axes. An axis with an infinite niche breadth
#' (`sigma = Inf`) exerts no selection: its trait, such as a shape that does not affect
#' growth, changes only through its genetic covariance with the selected traits and through
#' drift, the correlated response of Lande (1979). Founders draw their means from
#' \eqn{N(\bar z_{in}, G / k)}{N(zbar_in, G / k)} plus the variance among their sources trait by
#' trait. The environment can change on every axis (`env_time` as a list of matrices), and
#' [simulate_evolving_metacommunity()] takes several traits in the same way. With \eqn{G = 0}{G = 0} the
#' abundances are those without evolution for the same seed.
#'
#' @param genetic_var Additive genetic variance of the optimum \eqn{G}{G}, in squared units of
#'   the environment, one value or one per species. For several traits (see the section
#'   Several traits), a covariance matrix with one row and column per trait, or a list with
#'   one per species, such as from [genetic_matrices()].
#' @param ne_ratio Ratio of effective to census population size \eqn{\nu}{nu}; `Inf` turns drift
#'   off.
#' @param plasticity Slope of the reaction norm b, one value or one per species: 0 for a fixed
#'   phenotype, 1 for a phenotype that tracks the environment fully around \eqn{E_0}{E0}. For
#'   several traits, one value, one per trait or a matrix of species x traits; the plastic
#'   response of each trait is to its own axis.
#' @param env_var Environmental (non-heritable) variance of the phenotype \eqn{V_E}{VE}, one
#'   value or one per species; for several traits, a covariance matrix (or one per species),
#'   which need not be proportional to `genetic_var`.
#' @param steps_per_generation Time steps per generation \eqn{g}{g}; selection and drift per
#'   step are divided by it (at least 1); dispersal and demography are per step.
#' @param founder Should populations founded at a step start from the sampled mean of their
#'   founders (`TRUE`) or from the mean of their source (`FALSE`)?
#' @param immigrant_selection Should immigrants face viability selection in the site they reach
#'   before they mix with residents? With `FALSE` (the default), the immigrants from all sources
#'   join the population and the mean of the merged population is selected. With `TRUE`, the
#'   immigrants from each source arrive in proportion to the mean fitness of their source
#'   phenotype in the new site, relative to a perfectly adapted population (the same Gaussian
#'   as the growth rate, with the variance of the phenotypes), so that maladapted immigrants
#'   are filtered source by source, as individuals are in an individual-based model, and the
#'   mean and variance of the immigrants are those of the survivors. One trait only, on a
#'   spatially explicit landscape with `immigration = "poisson"`; used in the article
#'   "Replicating Leibold et al. (2019)".
#' @param mutation_var `NULL` (the default) for a constant genetic variance, or the mutational
#'   variance \eqn{V_m}{Vm} added per generation, one value or one per species. Each population
#'   then carries its own genetic variance, which loses a fraction \eqn{1 / (2 N_e)}{1 / (2 Ne)}
#'   and gains \eqn{V_m}{Vm} each generation (neutral mutation and drift, Lynch & Hill 1986),
#'   so that it approaches \eqn{2 N_e V_m}{2 Ne Vm} in a population of constant size and stays
#'   low in small ones. Founders keep the variance of their sources times
#'   \eqn{1 - 1 / (2k)}{1 - 1 / (2k)} for \eqn{k}{k} founders, and gene flow averages the
#'   variances of residents and immigrants by their numbers. The variance between the means of
#'   residents and immigrants (or of sources) is not added, as if selection removed it, since
#'   selection does not otherwise act on the variance; without this choice, gene flow between
#'   diverged populations inflates the variance without bound. `genetic_var` is then the
#'   variance of the populations at the start and of seeds. Needs a finite `ne_ratio`; one trait only.
#' @param variance_selection With `mutation_var`, should stabilising selection also remove
#'   genetic variance? With `TRUE`, each generation loses \eqn{G^2 / (\sigma^2 + G + V_E)}{G^2 /
#'   (sigma^2 + G + VE)}, the reduction of the genetic variance by Gaussian selection of width
#'   \eqn{\sigma}{sigma} on a phenotype of variance \eqn{G + V_E}{G + VE} (Bulmer 1971). Under
#'   sexual reproduction most of this reduction is linkage disequilibrium, which recombination
#'   restores, so the full loss holds for clonal inheritance, as in the model of Urban & De
#'   Meester (2009); with mutation, the variance then settles near the balance between
#'   mutation and selection (Lande 1975; Turelli 1984) rather than at \eqn{2 N_e V_m}{2 Ne Vm},
#'   and standing variation is eroded instead of carrying a load. `FALSE` (the default) keeps
#'   the mutation-drift dynamics alone.
#' @param env,traits The environment of each site (a vector, or a function of the step that
#'   returns one) and the traits of the species (from [species_traits()]). Filled in by
#'   [simulate_metacommunity()]; needed only with [mc_dynamics()].
#' @return An `mc_evolution` object, to pass as `evolution`.
#' @references Lande, R. (1976). Natural selection and random genetic drift in phenotypic
#'   evolution. *Evolution*, 30, 314--334. \doi{10.1111/j.1558-5646.1976.tb00911.x}
#'
#'   Lande, R. (1979). Quantitative genetic analysis of multivariate evolution, applied to
#'   brain:body size allometry. *Evolution*, 33, 402--416.
#'   \doi{10.1111/j.1558-5646.1979.tb04694.x}
#'
#'   Bulmer, M. G. (1971). The effect of selection on genetic variability. *The American
#'   Naturalist*, 105, 201--211. \doi{10.1086/282718}
#'
#'   Lande, R. (1975). The maintenance of genetic variability by mutation in a polygenic
#'   character with linked loci. *Genetical Research*, 26, 221--235.
#'   \doi{10.1017/S0016672300016037}
#'
#'   Turelli, M. (1984). Heritable genetic variation via mutation-selection balance: Lerch's
#'   zeta meets the abdominal bristle. *Theoretical Population Biology*, 25, 138--193.
#'   \doi{10.1016/0040-5809(84)90017-0}
#'
#'   Lynch, M. & Hill, W. G. (1986). Phenotypic evolution by neutral mutation. *Evolution*, 40,
#'   915--935. \doi{10.1111/j.1558-5646.1986.tb00561.x}
#' @family simulate metacommunities
#' @examples
#' set.seed(1)
#' ls <- landscape_implicit(20)
#' tr <- treeshapesim::sim_tree_shape(10)
#' traits <- species_traits(tr)
#' sim <- simulate_metacommunity(ls, traits, competition_matrix(tr, "phylo"),
#'                               dispersal = 0.001, n_steps = 300,
#'                               evolution = evolution_control(genetic_var = 4, ne_ratio = 0.5))
#' # local genetic means follow the environment of each site
#' head(sim$evolution$genetic[, 1:4])
#' @export
evolution_control <- function(genetic_var, ne_ratio = 1, plasticity = 0, env_var = 0,
                              steps_per_generation = 1, founder = TRUE, env = NULL,
                              traits = NULL, immigrant_selection = FALSE, mutation_var = NULL,
                              variance_selection = FALSE) {
  if (missing(genetic_var)) stop("`genetic_var` is needed", call. = FALSE)
  # several traits: covariance matrices (checked when the simulation starts)
  covariance <- function(v) is.matrix(v) || is.list(v)
  stopifnot(covariance(genetic_var) || (is.numeric(genetic_var) && all(genetic_var >= 0)),
            covariance(env_var) || (is.numeric(env_var) && all(env_var >= 0)),
            is.numeric(ne_ratio), length(ne_ratio) == 1, ne_ratio > 0,
            is.numeric(plasticity),
            is.numeric(steps_per_generation), length(steps_per_generation) == 1,
            steps_per_generation >= 1, is.logical(founder), length(founder) == 1,
            is.logical(immigrant_selection), length(immigrant_selection) == 1,
            !is.na(immigrant_selection))
  if (!is.null(mutation_var)) {
    if (!is.numeric(mutation_var) || any(is.na(mutation_var)) || any(mutation_var < 0)) {
      stop("`mutation_var` must be one or more non-negative numbers", call. = FALSE)
    }
    if (!is.finite(ne_ratio)) {
      stop("`mutation_var` needs drift: a finite `ne_ratio`", call. = FALSE)
    }
  }
  if (!is.logical(variance_selection) || length(variance_selection) != 1 ||
      is.na(variance_selection)) {
    stop("`variance_selection` must be TRUE or FALSE", call. = FALSE)
  }
  if (variance_selection && is.null(mutation_var)) {
    stop("`variance_selection` needs `mutation_var`", call. = FALSE)
  }
  structure(list(genetic_var = genetic_var, ne_ratio = ne_ratio, plasticity = plasticity,
                 env_var = env_var, steps_per_generation = steps_per_generation,
                 founder = founder, env = env, traits = traits,
                 immigrant_selection = immigrant_selection, mutation_var = mutation_var,
                 variance_selection = variance_selection),
            class = "mc_evolution")
}

#' @export
print.mc_evolution <- function(x, ...) {
  show <- function(v) {
    if (is.list(v)) return(paste(length(v), "matrices"))
    if (is.matrix(v)) return(paste0(nrow(v), " x ", ncol(v), " matrix"))
    format(v)
  }
  cat("<mc_evolution> G =", show(x$genetic_var), "| Ne/N =", format(x$ne_ratio),
      "| plasticity =", show(x$plasticity), "| VE =", show(x$env_var),
      "| steps per generation =", format(x$steps_per_generation),
      "| founder sampling:", x$founder, "\n")
  invisible(x)
}

# the state of evolution inside mc_dynamics(): checks the settings against the sites and
# species and returns the functions that the loop calls
evolution_state <- function(evolution, m, k_all, species) {
  if (!inherits(evolution, "mc_evolution")) {
    stop("`evolution` must come from `evolution_control()`", call. = FALSE)
  }
  if (is.null(evolution$env) || is.null(evolution$traits)) {
    stop("`evolution` needs `env` and `traits` with `mc_dynamics()`", call. = FALSE)
  }
  traits <- as_traits(evolution$traits)
  if (nrow(traits) != k_all) {
    stop("`traits` of `evolution` need one row per species (", k_all, ")", call. = FALSE)
  }
  if (!is.null(species) && all(species %in% rownames(traits))) traits <- traits[species, ]
  if (is.matrix(evolution$genetic_var) || is.list(evolution$genetic_var) ||
      is.matrix(evolution$env_var) || is.list(evolution$env_var)) {
    stop("covariance matrices in `evolution` need a niche of several axes (a list of ",
         "`traits`), one trait per axis", call. = FALSE)
  }
  env_fun <- if (is.function(evolution$env)) evolution$env else {
    function(step) evolution$env
  }
  e1 <- env_fun(1)
  if (length(e1) != m) stop("`env` of `evolution` needs one value per site (", m, ")",
                            call. = FALSE)
  each <- function(x) {
    if (!length(x) %in% c(1, k_all)) stop("settings of `evolution` take one value or one per ",
                                          "species (", k_all, ")", call. = FALSE)
    rep_len(x, k_all)
  }
  list(env = function(step) unname(as.numeric(env_fun(step))), ref = mean(e1),
       optimum = traits$optimum, sigma = traits$sigma, r_max = traits$r_max,
       g = each(evolution$genetic_var), b = each(evolution$plasticity),
       ve = each(evolution$env_var), nu = evolution$ne_ratio,
       per = evolution$steps_per_generation, founder = evolution$founder,
       imm_sel = isTRUE(evolution$immigrant_selection),
       dyn = !is.null(evolution$mutation_var),
       vm = if (is.null(evolution$mutation_var)) NULL else each(evolution$mutation_var),
       vsel = isTRUE(evolution$variance_selection))
}

# immigrants filtered by viability selection in the site they reach (immigrant_selection =
# TRUE): the expected number from each source is weighted by the fitness of the source
# phenotype in the destination; returns the expected numbers and the mean and variance of the
# genetic means, and the mean phenotype, of the immigrants that arrive (sites x species)
immigrants_selected <- function(em, za, pa, e, width, disp_matrix, vg = NULL) {
  m <- nrow(em)
  k <- ncol(em)
  if (!is.matrix(width)) width <- matrix(width, m, k, byrow = TRUE)
  if (is.null(vg)) vg <- matrix(0, m, k)
  lam <- zi <- z2 <- pi <- vi <- matrix(0, m, k)
  for (s in seq_len(k)) {
    src <- em[, s] > 0
    if (!any(src)) next
    # w[j, i]: fitness of the phenotype of source j in site i, relative to the optimum
    w <- exp(-outer(pa[src, s], e, "-")^2 / (2 * width[src, s]))
    flow <- disp_matrix[src, , drop = FALSE] * w * em[src, s]
    lam[, s] <- colSums(flow)
    zi[, s] <- colSums(flow * za[src, s])
    z2[, s] <- colSums(flow * za[src, s]^2)
    pi[, s] <- colSums(flow * pa[src, s])
    vi[, s] <- colSums(flow * vg[src, s])
  }
  has <- lam > 0
  zi[has] <- zi[has] / lam[has]
  pi[has] <- pi[has] / lam[has]
  vi[has] <- vi[has] / lam[has]
  z2[has] <- pmax(z2[has] / lam[has] - zi[has]^2, 0)
  zi[!has] <- za[!has]
  pi[!has] <- pa[!has]
  vi[!has] <- vg[!has]
  list(lambda = lam, mean = zi, between = z2, phenotype = pi, variance = vi)
}

# growth rates of the phenotypes p (sites x species) in environment e: the mean fitness under
# the response of growth_rates() of phenotypes with variance v around p (v = 0 gives it exactly)
# the same with the variance v as a matrix (sites x species), one value per population
evolution_growth_cells <- function(p, e, sigma, r_max, v) {
  width <- sweep(v, 2, sigma^2, "+")
  r <- exp(-(e - p)^2 / (2 * width))
  shrink <- sqrt(sweep(1 / width, 2, sigma^2, "*"))
  shrink[, !is.finite(sigma)] <- 1
  sweep(r * shrink, 2, r_max, "*")
}

evolution_growth <- function(p, e, sigma, r_max, v = 0) {
  width <- sigma^2 + v
  r <- exp(-sweep((e - p)^2, 2, 2 * width, "/"))
  # a niche of infinite breadth: no selection and no cost of the variance
  shrink <- ifelse(is.finite(sigma), sqrt(sigma^2 / width), 1)
  sweep(r, 2, r_max * shrink, "*")
}

# mean trait of the immigrants into each site, from the emigrants em and their trait x
# (sites x species), through the dispersal matrix or, without it, from all other sites
immigrant_mean <- function(em, x, disp_matrix, fallback) {
  x[is.na(x)] <- 0
  if (!is.null(disp_matrix)) {
    num <- crossprod(disp_matrix, em * x)
    den <- crossprod(disp_matrix, em)
  } else {
    num <- sweep(-(em * x), 2, colSums(em * x), "+")
    den <- sweep(-em, 2, colSums(em), "+")
  }
  out <- num / den
  out[!is.finite(out)] <- fallback[!is.finite(out)]
  out
}

# the genetic variance g (sites x species) after a fraction f of a generation (a number or a
# matrix) in populations of size n: drift removes 1 / (2 Ne) and mutation adds Vm per generation
# (Lynch & Hill 1986) and, with `variance_selection`, Gaussian selection of width `width`
# (= g + sigma^2 + VE) removes g^2 / width (Bulmer 1971)
variance_update <- function(g, n, width, ev, alive, f) {
  m <- nrow(g)
  k <- ncol(g)
  vm <- matrix(ev$vm[alive], m, k, byrow = TRUE)
  keep <- pmax(1 - 1 / (2 * ev$nu * n), 0)^f
  sel <- if (isTRUE(ev$vsel)) g^2 / width else 0
  pmax(g * keep + f * (vm - sel), 0)
}
