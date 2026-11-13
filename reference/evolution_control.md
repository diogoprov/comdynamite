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
  traits = NULL,
  immigrant_selection = FALSE,
  mutation_var = NULL,
  variance_selection = FALSE
)
```

## Arguments

- genetic_var:

  Additive genetic variance of the optimum \\G\\, in squared units of
  the environment, one value or one per species. For several traits (see
  the section Several traits), a covariance matrix with one row and
  column per trait, or a list with one per species, such as from
  [`genetic_matrices()`](https://provetelab.org/comdynamite/reference/genetic_matrices.md).

- ne_ratio:

  Ratio of effective to census population size \\\nu\\; `Inf` turns
  drift off.

- plasticity:

  Slope of the reaction norm b, one value or one per species: 0 for a
  fixed phenotype, 1 for a phenotype that tracks the environment fully
  around \\E_0\\. For several traits, one value, one per trait or a
  matrix of species x traits; the plastic response of each trait is to
  its own axis.

- env_var:

  Environmental (non-heritable) variance of the phenotype \\V_E\\, one
  value or one per species; for several traits, a covariance matrix (or
  one per species), which need not be proportional to `genetic_var`.

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

- immigrant_selection:

  Should immigrants face viability selection in the site they reach
  before they mix with residents? With `FALSE` (the default), the
  immigrants from all sources join the population and the mean of the
  merged population is selected. With `TRUE`, the immigrants from each
  source arrive in proportion to the mean fitness of their source
  phenotype in the new site, relative to a perfectly adapted population
  (the same Gaussian as the growth rate, with the variance of the
  phenotypes), so that maladapted immigrants are filtered source by
  source, as individuals are in an individual-based model, and the mean
  and variance of the immigrants are those of the survivors. One trait
  only, on a spatially explicit landscape with
  `immigration = "poisson"`; used in the article "Replicating Leibold et
  al. (2019)".

- mutation_var:

  `NULL` (the default) for a constant genetic variance, or the
  mutational variance \\V_m\\ added per generation, one value or one per
  species. Each population then carries its own genetic variance, which
  loses a fraction \\1 / (2 N_e)\\ and gains \\V_m\\ each generation
  (neutral mutation and drift, Lynch & Hill 1986), so that it approaches
  \\2 N_e V_m\\ in a population of constant size and stays low in small
  ones. Founders keep the variance of their sources times \\1 - 1 /
  (2k)\\ for \\k\\ founders, and gene flow averages the variances of
  residents and immigrants by their numbers. The variance between the
  means of residents and immigrants (or of sources) is not added, as if
  selection removed it, since selection does not otherwise act on the
  variance; without this choice, gene flow between diverged populations
  inflates the variance without bound. `genetic_var` is then the
  variance of the populations at the start and of seeds. Needs a finite
  `ne_ratio`; one trait only.

- variance_selection:

  With `mutation_var`, should stabilising selection also remove genetic
  variance? With `TRUE`, each generation loses \\G^2 / (\sigma^2 + G +
  V_E)\\, the reduction of the genetic variance by Gaussian selection of
  width \\\sigma\\ on a phenotype of variance \\G + V_E\\ (Bulmer 1971).
  Under sexual reproduction most of this reduction is linkage
  disequilibrium, which recombination restores, so the full loss holds
  for clonal inheritance, as in the model of Urban & De Meester (2009);
  with mutation, the variance then settles near the balance between
  mutation and selection (Lande 1975; Turelli 1984) rather than at \\2
  N_e V_m\\, and standing variation is eroded instead of carrying a
  load. `FALSE` (the default) keeps the mutation-drift dynamics alone.

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

## Several traits

With a niche of several environmental axes (a list of `traits` and an
`env` with one column per axis in
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)),
each population has a vector of genetic means, one trait per axis: the
optimum on that axis. `genetic_var` and `env_var` are then covariance
matrices \\G\\ and \\E\\, with \\P = G + E\\. The niche breadths make
the diagonal matrix \\\Omega\\ of Gaussian selection, and growth and the
response to selection are the multivariate versions of the steps above
(Lande 1979): \$\$r = r\_{max} \sqrt{\frac{\|\Omega\|}{\|\Omega + P\|}}
\exp\left(-\frac{1}{2} d^\top (\Omega + P)^{-1} d\right), \quad \Delta
\bar z = \frac{1}{g} G (\Omega + P)^{-1} d + \epsilon, \quad \epsilon
\sim N\left(0, \frac{G}{g N_e}\right),\$\$ with \\d\\ the difference
between the environment and the phenotype on each axis and \\r\_{max}\\
the product over axes. An axis with an infinite niche breadth
(`sigma = Inf`) exerts no selection: its trait, such as a shape that
does not affect growth, changes only through its genetic covariance with
the selected traits and through drift, the correlated response of Lande
(1979). Founders draw their means from \\N(\bar z\_{in}, G / k)\\ plus
the variance among their sources trait by trait. The environment can
change on every axis (`env_time` as a list of matrices), and
[`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)
takes several traits in the same way. With \\G = 0\\ the abundances are
those without evolution for the same seed.

## References

Lande, R. (1976). Natural selection and random genetic drift in
phenotypic evolution. *Evolution*, 30, 314–334.
[doi:10.1111/j.1558-5646.1976.tb00911.x](https://doi.org/10.1111/j.1558-5646.1976.tb00911.x)

Lande, R. (1979). Quantitative genetic analysis of multivariate
evolution, applied to brain:body size allometry. *Evolution*, 33,
402–416.
[doi:10.1111/j.1558-5646.1979.tb04694.x](https://doi.org/10.1111/j.1558-5646.1979.tb04694.x)

Bulmer, M. G. (1971). The effect of selection on genetic variability.
*The American Naturalist*, 105, 201–211.
[doi:10.1086/282718](https://doi.org/10.1086/282718)

Lande, R. (1975). The maintenance of genetic variability by mutation in
a polygenic character with linked loci. *Genetical Research*, 26,
221–235.
[doi:10.1017/S0016672300016037](https://doi.org/10.1017/S0016672300016037)

Turelli, M. (1984). Heritable genetic variation via mutation-selection
balance: Lerch's zeta meets the abdominal bristle. *Theoretical
Population Biology*, 25, 138–193.
[doi:10.1016/0040-5809(84)90017-0](https://doi.org/10.1016/0040-5809%2884%2990017-0)

Lynch, M. & Hill, W. G. (1986). Phenotypic evolution by neutral
mutation. *Evolution*, 40, 915–935.
[doi:10.1111/j.1558-5646.1986.tb00561.x](https://doi.org/10.1111/j.1558-5646.1986.tb00561.x)

## See also

Other simulate metacommunities:
[`as_metacommunity()`](https://provetelab.org/comdynamite/reference/as_metacommunity.md),
[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md),
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
