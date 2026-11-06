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
  species = NULL
)
```

## Arguments

- r:

  Sites x species matrix of density-independent growth rates, or a
  function of the time step that returns such a matrix.

- alpha:

  Competition matrix from
  [`competition_matrix()`](https://provetelab.org/comdynamite/reference/competition_matrix.md),
  with species in the same order as the rows of `traits`.

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

## Value

List with `abundance` and `history`, as in
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).

## Examples

``` r
set.seed(1)
# one species alone converges to (r - 1) / alpha = (5 - 1) / 0.05 = 80
out <- mc_dynamics(matrix(5, 20, 1), matrix(0.05), dispersal = 0, n_steps = 200)
mean(out$abundance)
#> [1] 78.05
```
