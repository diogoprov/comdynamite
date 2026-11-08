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
  extirpation = 0,
  pool = NULL,
  initial = NULL,
  immigration = c("poisson", "multinomial"),
  env = NULL
)
```

## Arguments

- landscape:

  An `mc_landscape` object.

- traits:

  Species traits from
  [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md),
  or, for a niche with several environmental axes, a list of them, one
  per column of `env`, with the same species; see
  [`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md).

- alpha:

  Competition matrix from
  [`competition_matrix()`](https://provetelab.org/comdynamite/reference/competition_matrix.md).
  Species are matched to the rows of `traits` by name, or by position if
  `alpha` has no names.

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

- env:

  Optional environmental axes of each site, a matrix or data frame with
  one row per site and one column per element of `traits` (such as the
  `env` of several landscapes on the same grid). Needed when `traits` is
  a list; replaces the `env` of `landscape`. Not available with
  `env_time`.

## Value

An `mc_sim` object: a list with `abundance` (sites x species matrix at
the last step), `history` (sites x species x recorded steps array, or
`NULL`), and the `landscape`, `traits`, `alpha`, `dispersal`,
`extirpation` and `pool` used, `env` when given, and `tree` when the
traits came from
[`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md)
with a tree.

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
  dispersers. With `immigration = "multinomial"`, the emigrants of each
  species at each site are instead distributed among the other sites
  multinomially, with the probabilities of the kernel, which conserves
  them. It costs about as much as the Poisson draw at low dispersal and
  up to 20 times more at high dispersal, where the two give nearly the
  same results.

After dispersal, each population is lost with probability `extirpation`,
independently of its size. Thompson et al. (2020) used 0.002 in their
competition-colonisation scenario and no extirpation in the others.

During colonisation, Poisson(`seed_mean`) individuals of every species
are added to every site at each of the `seed_steps`. After the last of
these, species absent from all sites can never return and are no longer
computed.

Distributions from elsewhere, such as ranges, biogeographic
reconstructions or the output of gen3sis, enter as `pool` and `initial`.
With `pool`, a species can only occur where its value is `TRUE` (or 1):
seeds go only there, and individuals that disperse to other sites are
lost, as in a hostile matrix. A species without any site in `pool` dies
out. `initial` replaces the first draw of Poisson(`seed_mean`)
individuals with given abundances, or with one individual where a matrix
of presences has 1; seeding still follows `seed_steps`, so use
`seed_steps = numeric(0)` to start from `initial` alone. Both take one
row per site and one column per species, named by species. Rows are
matched by name when they have names (`site1`, `site2` and so on in
`simulate_metacommunity()`) and by position otherwise. Without them,
results with the same seed are the same as before they existed. A tree
from elsewhere, such as one simulated with gen3sis or a model of
state-dependent diversification, enters through
[`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md);
drop extinct lineages first, since they have no distribution.

## Article and archived code

The article of Thompson et al. (2020) and the code archived with it
([doi:10.5281/zenodo.3708350](https://doi.org/10.5281/zenodo.3708350) )
differ in points that change results. Beta is gamma / alpha in the code
and in their figures, but gamma - alpha in their Table 1. Their
environment is converted to ranks, so it is uniform between 0 and 1, and
the lowest dispersal is 0.00001, not 0.0001. Emigrants are Poisson in
the code and binomial in the article, which comdynamite follows. In
their competition-colonisation scenario, the code also halves the
immigrants. See the article "Validating the model of Thompson et al.
(2020)" on the package website for the complete list.

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
tr <- treeshapesim::sim_tree_shape(20)
sim <- simulate_metacommunity(ls, species_traits(tr),
                              competition_matrix(tr, "phylo"),
                              dispersal = 0.01, n_steps = 200)
sim
#> <mc_sim> 30 sites x 20 species
#> species persisting: 18 of 20 
#> mean local richness: 7.03 
# a niche with two environmental axes
ls2 <- landscape_implicit(30)
niche2 <- species_traits(tr, rmax_range = c(1, 1))
sim2 <- simulate_metacommunity(ls, list(species_traits(tr), niche2),
                               competition_matrix(tr, "phylo"), n_steps = 200,
                               env = cbind(ls$env, ls2$env))
```
