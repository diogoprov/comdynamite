# Simulate a metacommunity

Runs the process-based model of Thompson et al. (2020) on a landscape.
At each time step, abundances are drawn as
`N' ~ Poisson(r N / (1 + N alpha'))`, where `r` holds the growth rates
of each species at each site
([`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md))
and `alpha` the competition coefficients
([`competition_matrix()`](https://provetelab.org/comdynamite/reference/competition_matrix.md)).
Each individual then emigrates with probability `dispersal`, and
emigrants are redistributed among sites:

## Usage

``` r
simulate_metacommunity(
  landscape,
  traits,
  alpha,
  dispersal = 0.01,
  n_steps = 1000,
  env_time = NULL,
  seed_steps = seq(10, 100, 10),
  seed_mean = 0.5,
  record_steps = NULL,
  extirpation = 0
)
```

## Arguments

- landscape:

  An `mc_landscape` object.

- traits:

  Species traits from
  [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md).

- alpha:

  Competition matrix from
  [`competition_matrix()`](https://provetelab.org/comdynamite/reference/competition_matrix.md),
  with species in the same order as the rows of `traits`.

- dispersal:

  Emigration probability, one value or one per species.

- n_steps:

  Number of time steps.

- env_time:

  Optional matrix (sites x `n_steps`) with the environment at each step,
  e.g. from
  [`env_over_time()`](https://provetelab.org/comdynamite/reference/env_over_time.md).
  If `NULL`, the environment of the landscape is constant.

- seed_steps:

  Steps at which individuals are added during colonisation.

- seed_mean:

  Mean number of individuals of each species added to each site at each
  seeding step (and at the start).

- record_steps:

  Optional steps at which abundances are recorded.

- extirpation:

  Probability that a population is lost at each time step, one value or
  one per species.

## Value

An `mc_sim` object: a list with `abundance` (sites x species matrix at
the last step), `history` (sites x species x recorded steps array, or
`NULL`), and the `landscape`, `traits`, `alpha`, `dispersal` and
`extirpation` used.

## Details

- spatially implicit landscapes
  ([`landscape_implicit()`](https://provetelab.org/comdynamite/reference/landscape_implicit.md)):
  the emigrants of each species are distributed multinomially, each site
  receiving them with probability proportional to the number of
  emigrants from the other sites (a uniform kernel; an individual may
  return to its own site);

- spatially explicit landscapes
  ([`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md)):
  the number of immigrants of each species at each site is drawn from a
  Poisson distribution whose mean is the expected flow of emigrants
  through the dispersal kernel. This approximation does not conserve the
  number of dispersing individuals exactly, but has the right
  expectation. The code of Thompson et al. (2020) does the same. Few
  emigrants may then produce no immigrant, so with `extirpation` and low
  dispersal, recolonisation is rarer than in mcomsimr, which conserves
  dispersers.

After dispersal, each population is lost with probability `extirpation`,
independently of its size. Thompson et al. (2020) used 0.002 in their
competition-colonisation scenario and no extirpation in the others.

During colonisation, Poisson(`seed_mean`) individuals of every species
are added to every site at each of the `seed_steps`. After the last of
these, species absent from all sites can never return and are no longer
computed.

## References

Thompson, P. L., Guzman, L. M., De Meester, L., Horváth, Z., Ptacnik,
R., Vanschoenwinkel, B., Viana, D. S., & Chase, J. M. (2020). A
process-based metacommunity framework linking local and regional scale
community ecology. *Ecology Letters*, 23, 1314–1329.
[doi:10.1111/ele.13568](https://doi.org/10.1111/ele.13568)

## Examples

``` r
set.seed(1)
ls <- landscape_implicit(30)
tr <- ape::rcoal(20)
sim <- simulate_metacommunity(ls, species_traits(tr),
                              competition_matrix(tr, "phylo"),
                              dispersal = 0.01, n_steps = 200)
sim
#> <mc_sim> 30 sites x 20 species
#> species persisting: 17 of 20 
#> mean local richness: 8 
```
