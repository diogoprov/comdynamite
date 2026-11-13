# Metacommunity dynamics from growth rates

The engine of
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
for users who want to supply their own growth rates (e.g. responses to
several environmental variables).

## Usage

``` r
mc_dynamics(
  r,
  alpha,
  dispersal,
  n_steps = 1000,
  seed_steps = seq(10, 100, 10),
  seed_mean = 0.5,
  disp_matrix = NULL,
  record_steps = NULL,
  sites = NULL,
  species = NULL,
  extirpation = 0,
  pool = NULL,
  initial = NULL,
  immigration = c("poisson", "multinomial"),
  evolution = NULL,
  disturbance = 0,
  pool_from = NULL,
  recruitment = c("beverton_holt", "lottery"),
  microsites = NULL,
  mortality = 1
)
```

## Arguments

- r:

  Sites x species matrix of density-independent growth rates, or a
  function of the time step that returns such a matrix.

- alpha:

  Competition matrix from
  [`competition_matrix()`](https://provetelab.org/comdynamite/reference/competition_matrix.md).
  Species are matched to the rows of `traits` by name, or by position if
  `alpha` has no names.

- dispersal:

  Emigration probability, one value or one per species.

- n_steps:

  Number of time steps.

- seed_steps:

  Steps at which individuals are added during colonisation.

- seed_mean:

  Mean number of individuals of each species added to each site at each
  seeding step (and at the start).

- disp_matrix:

  Row-stochastic dispersal matrix (sites x sites), or `NULL` for uniform
  dispersal.

- record_steps:

  Optional steps at which abundances are recorded.

- sites, species:

  Optional names of sites and species.

- extirpation:

  Probability that a population is lost at each time step, one value or
  one per species.

- pool:

  Optional matrix (sites x species), logical or 0/1, with the sites
  where each species can occur.

- initial:

  Optional matrix (sites x species) of abundances or 0/1 presences at
  the start, only where `pool` allows.

- immigration:

  How emigrants reach other sites on spatially explicit landscapes:
  `"poisson"` (the default, as in the code of Thompson et al. 2020) or
  `"multinomial"`, which conserves them. Spatially implicit landscapes
  are always multinomial.

- evolution:

  Optional settings from
  [`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md)
  for the evolution of the niche optimum within each population:
  selection by the local environment, drift, gene flow, founder sampling
  and plasticity. `NULL` (the default) keeps the optimum of each species
  in all its populations, as before. With a niche of several axes, one
  trait evolves per axis, with covariance matrices \\G\\ and \\E\\ (see
  the section Several traits of
  [`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md)
  and
  [`genetic_matrices()`](https://provetelab.org/comdynamite/reference/genetic_matrices.md)).

- disturbance:

  Probability that a site loses all its populations at each time step
  (one value; 0, the default, for none).

- pool_from:

  Optional matrix (sites x species) with the first step at which each
  species can occur in each site (1 for the start), for example to open
  a new site to two species at different times, as in the
  colonisation-time experiments of Urban & De Meester (2009). Combined
  with `pool`. Until then, seeds of that species do not establish there
  and the site is not a destination of its emigrants: those that would
  go there stay in their site, so that a closed site is not a sink
  (sites outside `pool` are).

- recruitment:

  How populations grow at each step: `"beverton_holt"` (the default),
  the Beverton-Holt competition of Thompson et al. (2020), or
  `"lottery"`, recruitment into a fixed number of microsites per site,
  as in the individual-based models of Urban & De Meester (2009) and
  Leibold et al. (2019). With a lottery, each adult survives the step
  with probability `1 - mortality`; adults produce Poisson offspring
  with mean equal to their growth rate (fecundity times the fitness of
  their phenotype, so that selection acts on recruits); offspring
  emigrate with probability `dispersal`; and the microsites left vacant
  by the survivors in each site are filled by a draw without replacement
  among the offspring present there (local and immigrant), so that
  emigration costs the site of origin little while offspring outnumber
  the vacancies. Competition is through the shared microsites, so
  `alpha` is not used. One trait at most with `evolution`; not available
  in
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md).

- microsites:

  With `recruitment = "lottery"`, the number of microsites of each site
  (one value or one per site). `initial` must fit in them; seeds take
  vacant microsites by lottery with the residents. In the colonisation
  history, `immigrants` and the migrant load count immigrant recruits,
  and `propagule` the individuals of a new population after recruitment.

- mortality:

  With `recruitment = "lottery"`, the probability that an adult dies at
  each step, whatever its phenotype: 1 (the default) for discrete
  generations, lower values for overlapping generations, with a
  generation time of about `1 / mortality` steps. With `evolution`,
  selection, mutation and the loss of variance to selection act on
  recruits, so that their effect per step scales with the share of
  recruits in the population; drift comes from the sampling of both
  survivors and recruits, so that it approaches that of a Moran model at
  low mortality; founders are sampled once, when they recruit;
  `steps_per_generation` must stay 1.

## Value

List with `abundance`, `history`, the colonisation history (`arrival`,
`propagule`, `immigrants`, `migrant_load`) and `n_steps`, as in
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
and `pool` as a logical matrix, or `NULL` without it.

## See also

Other simulate metacommunities:
[`as_metacommunity()`](https://provetelab.org/comdynamite/reference/as_metacommunity.md),
[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md),
[`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md),
[`genetic_matrices()`](https://provetelab.org/comdynamite/reference/genetic_matrices.md),
[`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md),
[`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md),
[`maladaptation()`](https://provetelab.org/comdynamite/reference/maladaptation.md),
[`persistence()`](https://provetelab.org/comdynamite/reference/persistence.md),
[`range_history()`](https://provetelab.org/comdynamite/reference/range_history.md),
[`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md),
[`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md),
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
[`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md),
[`speciation_control()`](https://provetelab.org/comdynamite/reference/speciation_control.md)

## Examples

``` r
set.seed(1)
# one species alone converges to (r - 1) / alpha = (5 - 1) / 0.05 = 80
out <- mc_dynamics(matrix(5, 20, 1), matrix(0.05), dispersal = 0, n_steps = 200)
mean(out$abundance)
#> [1] 78.05
```
