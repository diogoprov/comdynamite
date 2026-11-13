# Simulate a metacommunity with speciation

Runs the model of
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
with the evolution of the niche optimum within populations of
[`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md),
on species that bud as their populations diverge, by point mutation or
by fission
([`speciation_control()`](https://provetelab.org/comdynamite/reference/speciation_control.md)),
and die out when their last population is lost. The tree is not given:
it grows out of the dynamics. The community starts from the species of
`traits` (usually one ancestor) and is closed: there is no immigration
from outside the landscape.

## Usage

``` r
simulate_evolving_metacommunity(
  landscape,
  traits,
  speciation,
  evolution = NULL,
  dispersal = 0.01,
  extirpation = 0,
  n_steps = 1000,
  competition = c("phylo", "none"),
  inter_max = 1,
  inter_min = 0.25,
  scale = 0.05,
  distance_scale = 100,
  niche_overlap = TRUE,
  initial = NULL,
  seed_mean = 0.5,
  record_steps = NULL,
  immigration = c("poisson", "multinomial"),
  max_species = 500
)
```

## Arguments

- landscape:

  An `mc_landscape` object.

- traits:

  Traits of the species at the start, from
  [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md)
  (one row per species).

- speciation:

  Settings from
  [`speciation_control()`](https://provetelab.org/comdynamite/reference/speciation_control.md).

- evolution:

  Settings from
  [`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md);
  `NULL` keeps the optimum of each species fixed (`genetic_var = 0`),
  which rules out ecological speciation.

- dispersal:

  Emigration probability, one value or one per species.

- extirpation:

  Probability that a population is lost at each time step, one value or
  one per species.

- n_steps:

  Number of time steps.

- competition:

  `"phylo"` or `"none"`.

- inter_max, inter_min, scale, distance_scale:

  Competition between species, relative to competition within a species:
  between sisters just split (`inter_max`) and between distant species
  (`inter_min`); the overall scale (competition within a species); and
  the distance in generations over which the excess over `inter_min`
  decays by a factor of \\e\\.

- niche_overlap:

  Should competition between species also fall with the difference of
  their phenotypes in each site?

- initial:

  Optional matrix (sites x species) of abundances at the start; by
  default Poisson(`seed_mean`) individuals of each species in every
  site.

- seed_mean:

  Mean number of individuals of each species added to each site at each
  seeding step (and at the start).

- record_steps:

  Optional steps at which abundances are recorded.

- immigration:

  How emigrants reach other sites on spatially explicit landscapes:
  `"poisson"` (the default, as in the code of Thompson et al. 2020) or
  `"multinomial"`, which conserves them. Spatially implicit landscapes
  are always multinomial.

- max_species:

  Maximum number of species ever, including extinct ones; the simulation
  stops with an error beyond it.

## Value

An `mc_evolving` object, a list with:

- `abundance`: sites x extant species. With protracted speciation
  (`completion` of
  [`speciation_control()`](https://provetelab.org/comdynamite/reference/speciation_control.md)),
  the species are the good species, each with the individuals of its
  extant lineages; otherwise every lineage is a species.

- `tree`: the reconstructed tree of those species, with branch lengths
  in generations, nodes at times of initiation and a root edge from the
  start to the crown, or `NULL` with fewer than two.

- `tree_full`: the tree of all lineages, extinct ones included.

- `species`: a table with the parent, birth (initiation), end, fate,
  order of budding, time of completion and status (`"good"` or
  `"incipient"` at its end) of every lineage.

- `events`: one row per budding, with the new lineages.

- `lineages` and `lineage_species`: abundances of the extant lineages
  and the species each belongs to.

- `richness` and `richness_good`: extant lineages, and extant lineages
  that have completed, at each generation (the second can be below the
  number of good species, some of which survive only in incipient
  lineages).

- `genetic`: genetic means of the populations of the extant lineages.

- `traits`: of the extant species, with the mean genetic optimum over
  the populations of their lineages as `optimum`.

- `alpha`: among the extant species at the end, without the niche
  overlap, which varies among sites.

- The records and colonisation history of
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  for the extant species: records summed over their extant lineages
  (lineages extinct by the end are left out); colonisation, in each
  site, of the most abundant lineage.

- `series`: one row per lineage alive in each generation, with its
  `occupancy` (sites), regional `abundance`, `genetic` mean (weighted by
  abundance), and the sites it `colonised` and lost (`extirpated`) since
  the previous generation, out of those where it was absent
  (`absent_before`) and present (`present_before`). Sites that a lineage
  passed to a new one are not counted as lost, and the transitions are
  `NA` in the generation of birth.
  [`persistence()`](https://provetelab.org/comdynamite/reference/persistence.md)
  turns them into rates of colonisation and extirpation.

- `ranges`: the changes of presence of every lineage in every site, one
  row per change (`generation`, `lineage`, `site` and `present`), with
  the starting sites at generation 0 and the sites a new lineage took
  from its mother at its birth.
  [`range_history()`](https://provetelab.org/comdynamite/reference/range_history.md)
  maps them on the tree as ranges of areas.

- `evolution`: with `evolution`, the state of the extant species at the
  end, for
  [`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md)
  with `source = "evolution"`: genetic means in each site (`genetic`,
  weighted over their lineages by abundance; `NA` where absent), the
  plastic part of the phenotype (`plastic`), the variance among the
  lineages of a species in each site (`within_var`, 0 for species of one
  lineage), the environment of each site at the end (`environment`), the
  same genetic means, variance among lineages and environment at
  `record_steps` for the lineages extant at the end (`history`,
  `within_history`, `environment_history`), the reference environment of
  plasticity (`reference`), the genetic and environmental variances and
  the slope of plasticity of each species (`genetic_var`, `env_var`,
  `plasticity`) and the `settings`; otherwise `NULL`.

- The `landscape`, `dispersal`, `extirpation`, `speciation` and
  `n_steps` used.

[`as_metacommunity()`](https://provetelab.org/comdynamite/reference/as_metacommunity.md)
turns it into an `mc_sim` for the other functions of the package.

## Details

Competition follows the tree as it grows and, with
`niche_overlap = TRUE`, the phenotypes in each site: \\\alpha\_{ii}\\ =
`scale` and, between species in site \\s\\, \$\$\alpha\_{ij}(s) = c \\
(q_0 + (q_1 - q_0) e^{-d\_{ij} / \delta}) \\ \exp\left(-\frac{(p\_{is} -
p\_{js})^2}{2(\sigma_i^2 + \sigma_j^2)}\right),\$\$ with \\c\\ =
`scale`, \\q_1\\ = `inter_max`, \\q_0\\ = `inter_min`, \\d\_{ij}\\ the
patristic distance in generations, \\\delta\\ = `distance_scale`, \\p\\
the phenotypes (genetic mean and plastic part) and \\\sigma\\ the niche
breadths. The last factor is the overlap of two Gaussian niches, 1 for
equal phenotypes, so that species that have diverged in niche compete
less. With `inter_max = 1`, sisters just split with the same phenotype
compete as members of one species. With `inter_min = 0` and no niche
divergence, distant species do not compete at all, and richness is
limited only by speciation and extinction. With `competition = "none"`,
species do not compete with each other.

New species inherit the traits of their mother (`sigma`, `r_max`,
`dispersal`, `extirpation`, the settings of evolution) and, in each of
their sites, its abundance and genetic mean. Several species at the
start need a tree in `traits` (from
[`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md))
with branch lengths in generations, so that their distances are known.

Without a split, the simulation is the same as
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
with the same `evolution` and the same seed (with one species, or
several with `competition = "none"`).

## References

Hagen, O., Flück, B., Fopp, F., Cabral, J. S., Hartig, F., Pontarp, M.,
Rangel, T. F., & Pellissier, L. (2021). gen3sis: A general engine for
eco-evolutionary simulations of the processes that shape Earth's
biodiversity. *PLOS Biology*, 19, e3001340.
[doi:10.1371/journal.pbio.3001340](https://doi.org/10.1371/journal.pbio.3001340)

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
[`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md),
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
[`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md),
[`speciation_control()`](https://provetelab.org/comdynamite/reference/speciation_control.md)

## Examples

``` r
set.seed(1)
ls <- landscape_grid(6, disp_scale = 0.5)
ancestor <- data.frame(optimum = 50, sigma = 20, r_max = 5, row.names = "sp1")
sim <- simulate_evolving_metacommunity(
  ls, ancestor, dispersal = 0.002, n_steps = 400,
  speciation = speciation_control("geographic", threshold = 100),
  evolution = evolution_control(genetic_var = 20, ne_ratio = 0.5))
sim
#> <mc_evolving> 36 sites, 400 generations
#>   lineages ever: 247 | extant: 233 | extinct: 14 | buddings: 237 
#>   speciation: geographic 
```
