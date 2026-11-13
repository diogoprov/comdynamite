# Speciation settings

Settings for speciation in
[`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md).
Speciation is checked once per generation (`steps_per_generation` of
[`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md)).
A species buds when its populations fall into two or more groups: the
largest group (most sites, then most individuals) keeps the identity and
age of the species, and each other group becomes a new species with its
populations, abundances and genetic means. Several new species at once
branch off one after another, the most divergent first, at the same time
(internal branches of length zero), so that the tree stays binary.

## Usage

``` r
speciation_control(
  mode = c("geographic", "ecological", "point", "fission"),
  threshold,
  rate,
  divergence = c("cluster", "gradual"),
  migration_min = 0.005,
  reset_generations = 10,
  homogenization = 1,
  min_generations = 20,
  migration_max = 0.001,
  partition = c("spatial", "random"),
  jump_var = 0,
  completion = 0,
  completion_dist = c("fixed", "exponential")
)
```

## Arguments

- mode:

  `"geographic"`, `"ecological"`, `"point"` or `"fission"`.

- threshold:

  For `"geographic"` and `"ecological"`: divergence in generations, or
  difference of genetic means in units of the environment, that
  separates groups; `Inf` turns speciation off.

- rate:

  For `"point"`: probability per individual per generation; for
  `"fission"`: rate per species per generation.

- divergence:

  For `"geographic"`: `"cluster"` or `"gradual"`.

- migration_min:

  For `divergence = "cluster"`: migration rate per generation that links
  two populations. Where clusters form depends on dispersal, the kernel
  and abundances: the rate of the strongest link of each population
  should straddle it. In the design of
  `data-raw/evolving_migration_min.R` (strongest links with median
  0.006), species did not split at all at 0.002 or below and split
  almost as with no links at 0.02 or above.

- reset_generations:

  For `divergence = "cluster"`: generations in a row in the same cluster
  after which divergence returns to 0.

- homogenization:

  For `divergence = "gradual"`: reduction of divergence per migrant.

- min_generations:

  For `"ecological"`: generations in a row of difference and reduced
  gene flow before groups split.

- migration_max:

  For `"ecological"`: migration rate per generation below which gene
  flow counts as reduced.

- partition:

  For `"fission"`: `"spatial"` or `"random"`.

- jump_var:

  Variance of the shift of the genetic mean of each new species.

- completion:

  Time to complete speciation, in generations (fixed or mean); 0 for
  instantaneous speciation.

- completion_dist:

  `"fixed"` or `"exponential"`.

## Value

An `mc_speciation` object.

## Details

Four modes start new species (the initiation of Hubert et al. 2015):

- `"geographic"`: each pair of populations of a species has a divergence
  \\d\\, in generations. With `divergence = "cluster"` (the default,
  close to gen3sis; Hagen et al. 2021), populations linked by a
  migration rate of at least `migration_min` form clusters (chains of
  links). \\d\\ grows by 1 per generation between populations of
  different clusters; within a cluster it stays as it is, and returns to
  0 only after `reset_generations` generations in a row in the same
  cluster, so that a single generation of contact does not erase a long
  isolation. With `divergence = "gradual"`, \\d\\ grows by 1 per
  generation and gene flow lowers it in proportion to the migrants
  exchanged: \\d \leftarrow \max(0, d + 1 - h M)\\, with \\M\\ the
  expected migrants between the two populations in the generation and
  \\h\\ = `homogenization`. Populations are grouped by single linkage:
  two are in the same group when a chain of pairs with \\d \<\\
  `threshold` links them.

- `"ecological"`: each pair of populations counts the generations in a
  row in which their genetic means (`evolution`, with genetic variance)
  differ by at least `threshold` and their migration rate is below
  `migration_max`. Populations are grouped by single linkage on pairs
  whose count is below `min_generations`. Requiring the difference to
  last, with reduced gene flow, keeps drift and founder sampling in
  small populations from splitting species in a uniform environment.

- `"point"`: point mutation, the mode of neutral theory and of MESS
  (Overcast et al. 2021). Each individual becomes, with probability
  `rate` per generation, the first individual of a new species in its
  site. The mother keeps its populations, except where a mutant takes
  its last individual. Mutants arrive in the generation of their birth.

- `"fission"`: random fission. Each species with two or more populations
  splits with probability \\1 - e^{-\lambda}\\ per generation
  (\\\lambda\\ = `rate`). With `partition = "spatial"`, the populations
  of the new group are those nearest (by the dispersal kernel) to a
  random occupied site, in a random number from 1 to all but one; with
  `"random"` (or on a spatially implicit landscape), each population
  joins it with probability 1/2, with at least one population on each
  side.

Gene flow between two populations \\a\\ and \\b\\ is measured by the
migration rate, the larger of the expected migrants from one to the
other in the generation over the size of the receiving population, so
that it does not depend on abundance. A new population takes the
divergences (or counts) of the population that sent it most migrants,
and a population lost at any step of a generation founds anew.

New species take the genetic means of their populations; with
`jump_var > 0` their genetic means (and the optimum of their regional
pool) shift by one Gaussian draw with that variance per new species, as
in the punctuated change of MESS (Overcast et al. 2021).

**Protracted speciation.** With `completion > 0`, each new lineage is
incipient until it completes, after `completion` generations
(`completion_dist = "fixed"`) or after an exponential time with that
mean (`"exponential"`, the rate \\\lambda_2\\ of the protracted
birth-death model; Etienne & Rosindell 2012). Completion only changes
what is recognised as a species: incipient lineages are separate
lineages in the dynamics. At the end, a lineage that has completed is a
good species, and an incipient lineage belongs to the species its parent
belonged to at its birth, as in the labels of the PBD. The reconstructed
tree has one tip per good species with extant lineages, represented by
the good lineage itself or, if it is extinct, by its earliest extant
member. Its nodes are times of initiation of lineages, which for a
species whose lineages are not monophyletic may predate the initiation
of the species itself.

## References

Etienne, R. S., & Rosindell, J. (2012). Prolonging the past counteracts
the pull of the present: protracted speciation can explain observed
slowdowns in diversification. *Systematic Biology*, 61, 204–213.
[doi:10.1093/sysbio/syr091](https://doi.org/10.1093/sysbio/syr091)

Hagen, O., Flück, B., Fopp, F., Cabral, J. S., Hartig, F., Pontarp, M.,
Rangel, T. F., & Pellissier, L. (2021). gen3sis: A general engine for
eco-evolutionary simulations of the processes that shape Earth's
biodiversity. *PLOS Biology*, 19, e3001340.
[doi:10.1371/journal.pbio.3001340](https://doi.org/10.1371/journal.pbio.3001340)

Hubert, N., Calcagno, V., Etienne, R. S., & Mouquet, N. (2015).
Metacommunity speciation models and their implications for
diversification theory. *Ecology Letters*, 18, 864–881.
[doi:10.1111/ele.12458](https://doi.org/10.1111/ele.12458)

Overcast, I., Ruffley, M., Rosindell, J., Harmon, L., Borges, P. A. V.,
Emerson, B. C., Etienne, R. S., Gillespie, R., Krehenwinkel, H., Mahler,
D. L., Massol, F., Parent, C. E., Patiño, J., Peter, B., Week, B.,
Wagner, C., Hickerson, M. J., & Rominger, A. (2021). A unified model of
species abundance, genetic diversity, and functional diversity reveals
the mechanisms structuring ecological communities. *Molecular Ecology
Resources*, 21, 2782–2800.
[doi:10.1111/1755-0998.13514](https://doi.org/10.1111/1755-0998.13514)

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
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
[`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md)

## Examples

``` r
speciation_control("geographic", threshold = 50)
#> <mc_speciation> geographic | threshold = 50 | divergence: cluster | migration_min = 0.005 | reset_generations = 10
speciation_control("geographic", threshold = 50, divergence = "gradual", homogenization = 0.5)
#> <mc_speciation> geographic | threshold = 50 | divergence: gradual | homogenization = 0.5
speciation_control("ecological", threshold = 20, min_generations = 30)
#> <mc_speciation> ecological | threshold = 20 | min_generations = 30 | migration_max = 0.001
speciation_control("point", rate = 1e-4, completion = 50)
#> <mc_speciation> point | rate = 1e-04 | completion: fixed 50
speciation_control("fission", rate = 0.01, partition = "random")
#> <mc_speciation> fission | rate = 0.01 | partition: random
```
