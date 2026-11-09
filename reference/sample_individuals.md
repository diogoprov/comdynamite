# Individuals from a simulated metacommunity

Draws individuals from the populations of a simulated metacommunity and
gives each one a trait value, so that intraspecific variation and local
adaptation can be studied with the T-statistics of Violle et al. (2012)
([`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)).

## Usage

``` r
sample_individuals(
  x,
  n = 10,
  fraction = NULL,
  local_adaptation = 0,
  population_sd = 0,
  intraspecific_sd = 1,
  displacement = 0,
  displacement_width = NULL,
  step = NULL,
  env = NULL,
  cover = NULL,
  min_n = 1
)
```

## Arguments

- x:

  An `mc_sim` object from
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).

- n:

  Number of individuals measured in each population, or all of them when
  the population is smaller. Ignored when `fraction` is given.

- fraction:

  Optional probability that each individual is measured, so that the
  number measured in a population is Binomial(abundance, `fraction`),
  proportional to abundance on average.

- local_adaptation:

  Fraction (0 to 1) of the gap between the optimum of the species and
  the environment of the site that each population closes; one value or
  one per species (named by species, or in their order). The same holds
  for the two standard deviations.

- population_sd:

  Standard deviation of the deviations of populations from the species
  mean that do not track the environment; one value or one per species.

- intraspecific_sd:

  Standard deviation of trait values among the individuals of a
  population; one value or one per species.

- displacement:

  Strength of character displacement among the populations of a site
  (`c` in Details), relative to the selection that keeps each population
  near its mean without competitors; 0 (the default) for none.

- displacement_width:

  Width `w` of the Gaussian competition kernel of character
  displacement, in units of the trait. The default is the mean niche
  breadth (`sigma`) of the species.

- step:

  Recorded step to sample, as in `record_steps`; `NULL` (the default)
  uses the abundances at the last step.

- env:

  Optional environment of each site at the sampled step. The default is
  the environment of the landscape; give it for simulations with
  `env_time`.

- cover:

  Optional share (0 to 1) of the abundance of each site covered by the
  species measured: species are taken from the most abundant down until
  they reach `cover`, and the rest are not measured, as in field
  protocols that measure the species making up 80% of the abundance
  (Pakeman & Quested 2007). `NULL` (the default) measures every species.

- min_n:

  Populations with fewer individuals measured than `min_n` are left out,
  as when populations with a single individual are removed before
  computing T-statistics.

## Value

Data frame with one row per individual and the columns `site`,
`species`, `individual` (its number within the population), `abundance`
(of the population), `env` (of the site), `optimum` (of the species),
`population_mean` (the expected value of the population, \\m\_{ij}\\
above, after character displacement) and `trait`. Populations with no
individual measured are absent.

## Details

The trait is the environmental optimum of
[`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md),
the one that sets growth in
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).
The value of individual \\k\\ of species \\j\\ at site \\i\\ is

\$\$z\_{ijk} = \theta_j + h_j (E_i - \theta_j) + d\_{ij} + e\_{ijk},\$\$

where \\\theta_j\\ is the optimum of the species, \\E_i\\ the
environment of the site, \\h_j\\ (`local_adaptation`) the fraction of
the gap between them that the population closes, \\d\_{ij} \sim N(0,
s_p^2)\\ (with \\s_p\\ = `population_sd`) a deviation of the population
that does not track the environment (drift, founder effects), and
\\e\_{ijk} \sim N(0, s_e^2)\\ (with \\s_e\\ = `intraspecific_sd`) the
variation among individuals of a population. With `local_adaptation = 0`
and `population_sd = 0` all populations of a species share its mean;
with `local_adaptation = 1` every population sits on the environment of
its site and species no longer differ within a site.

With `displacement > 0`, the populations of each site also diverge from
one another, as in ecological character displacement. Their means are
the equilibrium of

\$\$\mu\_{ij} = m\_{ij} + c \sum_k p\_{ik} (\mu\_{ij} - \mu\_{ik})
\exp\left(-\frac{(\mu\_{ij} - \mu\_{ik})^2}{2 w^2}\right),\$\$

where \\m\_{ij} = \theta_j + h_j (E_i - \theta_j) + d\_{ij}\\ is the
mean above, \\c\\ is `displacement`, \\p\_{ik}\\ the share of species
\\k\\ among the individuals of the site and \\w\\ (`displacement_width`)
the width of a Gaussian competition kernel. Each population is pushed
away from the populations close to it in trait, more so from abundant
ones, and held near \\m\_{ij}\\, so \\c\\ sets the strength of
competition relative to the selection that keeps the population near
\\m\_{ij}\\. The push from one competitor is largest at a distance \\w\\
and at most \\c w e^{-1/2}\\, about \\0.6 c w\\. The equilibrium is
found by damped fixed-point iteration from \\m\_{ij}\\, with steps of
\\1 / (1 + c)\\. Populations with the same \\m\_{ij}\\, as with
`local_adaptation = 1` and no `population_sd`, are not pushed apart,
since the push is symmetric. With \\c \> 1\\, however, the equilibrium
with equal means is unstable: populations whose \\m\_{ij}\\ differ even
slightly end up about \\w \sqrt{2 \log c}\\ apart, and several
equilibria can exist, of which the iteration finds the one reached from
\\m\_{ij}\\. The size of the displacement therefore grows with \\c\\
only slowly once \\c\\ exceeds 1.

The traits of individuals are drawn after the simulation and do not feed
back on the dynamics: abundances are those of the species means. The
function describes the pattern that local adaptation would leave in
samples of individuals, not its demographic consequences.

## References

Violle, C., Enquist, B. J., McGill, B. J., Jiang, L., Albert, C. H.,
Hulshof, C., Jung, V., & Messier, J. (2012). The return of the variance:
intraspecific variability in community ecology. *Trends in Ecology &
Evolution*, 27, 244–252.
[doi:10.1016/j.tree.2011.11.014](https://doi.org/10.1016/j.tree.2011.11.014)

Pakeman, R. J., & Quested, H. M. (2007). Sampling plant functional
traits: what proportion of the species need to be measured? *Applied
Vegetation Science*, 10, 91–96.
[doi:10.1111/j.1654-109x.2007.tb00507.x](https://doi.org/10.1111/j.1654-109x.2007.tb00507.x)

## See also

[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md) for
the T-statistics;
[`subsample_individuals()`](https://provetelab.org/comdynamite/reference/subsample_individuals.md)
and
[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md)
for the same sampling schemes on real data.

Other trait variation among individuals:
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
[`neighbour_metrics`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
[`nested_variance()`](https://provetelab.org/comdynamite/reference/nested_variance.md),
[`null_intervals()`](https://provetelab.org/comdynamite/reference/null_intervals.md),
[`plot.mc_trait_gradient()`](https://provetelab.org/comdynamite/reference/plot.mc_trait_gradient.md),
[`plot.nested_variance()`](https://provetelab.org/comdynamite/reference/plot.nested_variance.md),
[`plot.trait_flex()`](https://provetelab.org/comdynamite/reference/plot.trait_flex.md),
[`plot.trait_rarefaction()`](https://provetelab.org/comdynamite/reference/plot.trait_rarefaction.md),
[`plot.tstats()`](https://provetelab.org/comdynamite/reference/plot.tstats.md),
[`plot_population_traits()`](https://provetelab.org/comdynamite/reference/plot_population_traits.md),
[`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md),
[`subsample_individuals()`](https://provetelab.org/comdynamite/reference/subsample_individuals.md),
[`trait_densities()`](https://provetelab.org/comdynamite/reference/trait_densities.md),
[`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md),
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md),
[`trait_gradient_brms()`](https://provetelab.org/comdynamite/reference/trait_gradient_brms.md),
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md),
[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md),
[`trait_gradient_weights()`](https://provetelab.org/comdynamite/reference/trait_gradient_weights.md),
[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
[`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md)

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(15)
sim <- simulate_metacommunity(landscape_implicit(10), species_traits(tree),
                              competition_matrix(tree, "none"), n_steps = 100,
                              seed_steps = 1:5)
ind <- sample_individuals(sim, n = 10, local_adaptation = 0.5, intraspecific_sd = 2)
head(ind)
#>    site species individual abundance      env  optimum population_mean    trait
#> 1 site1      t2          1        79 47.76196 46.10733        46.93465 50.64409
#> 2 site1      t2          2        79 47.76196 46.10733        46.93465 44.09239
#> 3 site1      t2          3        79 47.76196 46.10733        46.93465 48.50888
#> 4 site1      t2          4        79 47.76196 46.10733        46.93465 47.59176
#> 5 site1      t2          5        79 47.76196 46.10733        46.93465 43.93408
#> 6 site1      t2          6        79 47.76196 46.10733        46.93465 44.05996
```
