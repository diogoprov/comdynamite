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
  immigration = c("poisson", "multinomial")
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

## Value

List with `abundance`, `history`, the colonisation history (`arrival`,
`propagule`, `immigrants`, `migrant_load`) and `n_steps`, as in
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
and `pool` as a logical matrix, or `NULL` without it.

## See also

Other simulate metacommunities:
[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md),
[`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md),
[`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md),
[`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md),
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
[`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md)

## Examples

``` r
set.seed(1)
# one species alone converges to (r - 1) / alpha = (5 - 1) / 0.05 = 80
out <- mc_dynamics(matrix(5, 20, 1), matrix(0.05), dispersal = 0, n_steps = 200)
mean(out$abundance)
#> [1] 78.05
```
