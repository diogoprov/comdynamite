# Evolution of local populations

Settings for the evolution of the niche optimum within each population
(a species at a site) during
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
or
[`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md).
Without them, every population of a species keeps the optimum of the
species. With them, each population has its own genetic mean \\\bar z\\,
which responds to selection by the local environment, drifts, mixes with
immigrants and starts from a sample of its founders.

## Usage

``` r
evolution_control(
  genetic_var,
  ne_ratio = 1,
  plasticity = 0,
  env_var = 0,
  steps_per_generation = 1,
  founder = TRUE,
  env = NULL,
  traits = NULL
)
```

## Arguments

- genetic_var:

  Additive genetic variance of the optimum \\G\\, in squared units of
  the environment, one value or one per species.

- ne_ratio:

  Ratio of effective to census population size \\\nu\\; `Inf` turns
  drift off.

- plasticity:

  Slope of the reaction norm b, one value or one per species: 0 for a
  fixed phenotype, 1 for a phenotype that tracks the environment fully
  around \\E_0\\.

- env_var:

  Environmental (non-heritable) variance of the phenotype \\V_E\\, one
  value or one per species.

- steps_per_generation:

  Time steps per generation \\g\\; selection and drift per step are
  divided by it (at least 1); dispersal and demography are per step.

- founder:

  Should populations founded at a step start from the sampled mean of
  their founders (`TRUE`) or from the mean of their source (`FALSE`)?

- env, traits:

  The environment of each site (a vector, or a function of the step that
  returns one) and the traits of the species (from
  [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md)).
  Filled in by
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md);
  needed only with
  [`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md).

## Value

An `mc_evolution` object, to pass as `evolution`.

## Details

At each step, in the order of the model:

1.  **Seeds** carry the optimum of the species in `traits`, the mean of
    a regional pool that does not evolve. They mix with residents in
    proportion to their numbers; a population founded by seeds alone
    starts as below.

2.  **Phenotype and growth.** The phenotype is \\p = \bar z + b (E -
    E_0)\\, a linear reaction norm with slope `plasticity` (b) around
    the reference environment \\E_0\\, the mean environment of the sites
    at the first step. The plastic part is immediate and reversible:
    every individual, migrants included, grows with the phenotype of the
    site where it is. Growth is the mean fitness of a population whose
    phenotypes are normal around \\p\\ with variance \\P = G + V_E\\,
    under the Gaussian response of
    [`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md):
    \\r = r\_{max} \sqrt{\sigma^2 / (\sigma^2 + P)} \exp(-(E - p)^2 /
    (2(\sigma^2 + P)))\\. The variance widens the response and lowers
    its peak (the load); with \\P = 0\\ it is the response of
    [`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md).

3.  **Selection and drift**, in each population present after
    reproduction (\\N \> 0\\): \$\$\Delta \bar z = \frac{1}{g}
    \frac{G}{\sigma^2 + P} (E - p) + \epsilon, \quad \epsilon \sim
    N\left(0, \frac{G}{g N_e}\right), \quad N_e = \nu N,\$\$ the
    response to Gaussian stabilising selection of Lande (1976), with
    additive genetic variance \\G\\ (`genetic_var`), phenotypic variance
    \\P = G + V_E\\ (`env_var`), the niche breadth \\\sigma\\ of the
    species as the width of the fitness function, \\g\\ steps per
    generation (`steps_per_generation`) and the ratio \\\nu\\ of
    effective to census size (`ne_ratio`), the same fitness function as
    in growth. \\G\\ is constant: it does not erode with drift nor grow
    when divergent populations mix.

4.  **Gene flow.** After dispersal, the genetic mean of a population is
    the mean of residents and immigrants weighted by their numbers.
    Immigrants carry the mean of the populations they come from,
    weighted by the emigrants each sends through the dispersal kernel
    (or, on spatially implicit landscapes, by the emigrants of all other
    sites). A population with no residents after reproduction and
    emigration (absent, or lost all its members) is founded by its \\k\\
    immigrants: with `founder = TRUE`, their mean is drawn from \\N(\bar
    z\_{in}, (G + V_B)/k)\\, the sampling of a few individuals, where
    \\V_B\\ is the variance among the means of their sources, weighted
    as above. Populations present at the start and those founded by
    seeds alone are drawn in the same way around the optimum of the
    species. Immigrants into a population with residents carry the
    expected mean of their sources, without sampling noise. The
    phenotype recorded after dispersal (`phenotype` in the output) mixes
    residents and immigrants still carrying the plastic response to
    their site of origin, as when individuals are sampled right after
    dispersal; it does not enter the dynamics.

5.  Populations lost to extirpation, competition or the `pool` lose
    their mean.

Species with \\G = 0\\ keep the optimum of the species exactly and draw
no random numbers, so with `genetic_var = 0`, `env_var = 0` and
`plasticity = 0` the abundances are the same as without evolution for
the same seed. With
[`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md),
growth rates come from `env` and `traits` of `evolution`, and `r` only
sets the numbers of sites and species.

Two cases have known answers, checked in the tests of the package. A
population alone, without selection (very wide niche) and of constant
size \\N\\, accumulates a variance of \\G t / (\nu N)\\ after \\t\\
generations (Lande 1976). Two populations of the same size in
environments \\E_1\\ and \\E_2\\, exchanging a share \\m\\ of
individuals each generation, with selection \\a = G / (\sigma^2 + P)\\
and no drift, reach the divergence \\d^\* = (1 - 2m) a (E_1 - E_2) /
(1 - (1 - 2m)(1 - a))\\, which follows from the recursion of the steps
above.

## References

Lande, R. (1976). Natural selection and random genetic drift in
phenotypic evolution. *Evolution*, 30, 314–334.
[doi:10.1111/j.1558-5646.1976.tb00911.x](https://doi.org/10.1111/j.1558-5646.1976.tb00911.x)

## See also

Other simulate metacommunities:
[`as_metacommunity()`](https://provetelab.org/comdynamite/reference/as_metacommunity.md),
[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md),
[`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md),
[`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md),
[`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md),
[`persistence()`](https://provetelab.org/comdynamite/reference/persistence.md),
[`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md),
[`simulate_logistic()`](https://provetelab.org/comdynamite/reference/simulate_logistic.md),
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
[`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md),
[`speciation_control()`](https://provetelab.org/comdynamite/reference/speciation_control.md)

## Examples

``` r
set.seed(1)
ls <- landscape_implicit(20)
tr <- treeshapesim::sim_tree_shape(10)
traits <- species_traits(tr)
sim <- simulate_metacommunity(ls, traits, competition_matrix(tr, "phylo"),
                              dispersal = 0.001, n_steps = 300,
                              evolution = evolution_control(genetic_var = 4, ne_ratio = 0.5))
# local genetic means follow the environment of each site
head(sim$evolution$genetic[, 1:4])
#>             t1       t2       t3       t4
#> site1 33.82137 35.20115       NA 19.62512
#> site2 41.21913       NA 37.06979 23.92967
#> site3 57.64955 65.25984 57.58140 48.50037
#> site4 88.43806 91.00418 92.88088       NA
#> site5 28.77965 67.32722       NA 13.15587
#> site6 85.50623 91.24391 89.18069       NA
```
