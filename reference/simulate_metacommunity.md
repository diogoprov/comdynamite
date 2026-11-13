# Simulate a metacommunity

Runs the process-based model of Thompson et al. (2020) on a landscape.
At each time step, abundances are drawn as \\N' \sim \mathrm{Poisson}(r
N / (1 + N \alpha'))\\, where \\r\\ holds the growth rates of each
species at each site
([`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md))
and \\\alpha\\ the competition coefficients
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
  env = NULL,
  evolution = NULL,
  disturbance = 0,
  pool_from = NULL,
  recruitment = c("beverton_holt", "lottery"),
  microsites = NULL,
  mortality = 1
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
  If `NULL`, the environment of the landscape is constant. With `env`
  (several axes), a list with one such matrix per column of `env`
  (matched by name when named), or an array of sites x steps x axes.

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
  a list; replaces the `env` of `landscape`. With `env_time`, it gives
  the axes (names and number), and the environment of each step comes
  from `env_time`. Columns without names take the names of `traits`, or
  `axis1`, `axis2`, ...

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

An `mc_sim` object: a list with `abundance` (sites x species matrix at
the last step), `history` (sites x species x recorded steps array, or
`NULL`), and the `landscape`, `traits`, `alpha`, `dispersal`,
`extirpation` and `pool` used, `env` when given, and `tree` when the
traits came from
[`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md)
with a tree, the colonisation history (`arrival`, `propagule`,
`immigrants`, `migrant_load`; see the section Colonisation history) and
`n_steps`. With `evolution`, also `evolution`, a list with the genetic
mean of the optimum in each population (`genetic`), the plastic part of
its phenotype at the last step (`plastic`), the phenotype after
dispersal (`phenotype`, which mixes residents and immigrants), the
genetic means and the environment at `record_steps` (`history`,
`environment_history`), the reference environment and the settings, and,
with `mutation_var`, the genetic variance of each population at the last
step (`genetic_variance`); `NA` where a population is absent. With
several evolving traits, `genetic`, `plastic` and `phenotype` are arrays
of sites x species x traits, `history` has the recorded steps as a
fourth dimension, and `genetic_var`, `env_var` (one matrix per species)
and `plasticity` (species x traits) are kept. With a niche of several
axes,
[`as_mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md)
keeps every axis (`env` and one `optimum_` column of traits per axis),
while
[`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md)
uses the first axis unless given `trait` and `env`, and
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md)
needs a single axis unless the traits evolved (`source = "evolution"`,
one column per trait).

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
competition-colonisation scenario and no extirpation in the others. With
`disturbance`, each site then loses all its populations at once with
that probability, as the patch extinctions of Leibold et al. (2019); the
site is recolonised by dispersal.

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

## Colonisation history

The simulation records, for every population (a species at a site)
present at the last step, when and how it was founded, for
arrival-dependent trait deviations in
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md)
(founder effects and monopolization):

- `arrival`, the step at which the population was last founded: the
  first step at whose end it is present after being absent at the end of
  the previous one. Populations present at the end of step 1, from
  `initial` or the first draw of seeds, have arrival 1; a population
  lost to extirpation or competition and founded again gets the step of
  the new arrival.

- `propagule`, the number of founders: the seeds (or `initial`
  individuals) and immigrants that reached the site at the step of
  arrival.

- `immigrants`, the number of immigrants and seeds that reached the
  population after its arrival.

- `migrant_load`, the sum over the steps since arrival of the share of
  the population made up of newcomers: immigrants over the abundance
  after dispersal, plus seeds over the abundance before reproduction.
  Each newcomer carries the species mean, so under complete mixing a
  deviation of the population shrinks by a factor \\1 - m_t \approx
  e^{-m_t}\\ at a step where newcomers make up a share \\m_t\\, and by
  about \\e^{-L}\\ over the steps with migrant load \\L\\.

Each is a sites x species matrix, `NA` for populations absent at the
last step. Recording them takes no random draws, so abundances are the
same as without them for the same seed.
[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md)
returns them as a table.

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

Urban, M. C., & De Meester, L. (2009). Community monopolization: local
adaptation enhances priority effects in an evolving metacommunity.
*Proceedings of the Royal Society B*, 276, 4129–4138.
[doi:10.1098/rspb.2009.1382](https://doi.org/10.1098/rspb.2009.1382)

Leibold, M. A., Urban, M. C., De Meester, L., Klausmeier, C. A., &
Vanoverbeke, J. (2019). Regional neutrality evolves through local
adaptive niche evolution. *Proceedings of the National Academy of
Sciences*, 116, 2612–2617.
[doi:10.1073/pnas.1808615116](https://doi.org/10.1073/pnas.1808615116)

## See also

Other simulate metacommunities:
[`as_metacommunity()`](https://provetelab.org/comdynamite/reference/as_metacommunity.md),
[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md),
[`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md),
[`genetic_matrices()`](https://provetelab.org/comdynamite/reference/genetic_matrices.md),
[`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md),
[`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md),
[`maladaptation()`](https://provetelab.org/comdynamite/reference/maladaptation.md),
[`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md),
[`persistence()`](https://provetelab.org/comdynamite/reference/persistence.md),
[`range_history()`](https://provetelab.org/comdynamite/reference/range_history.md),
[`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md),
[`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md),
[`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md),
[`speciation_control()`](https://provetelab.org/comdynamite/reference/speciation_control.md)

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
