# Individuals from a simulated metacommunity

Draws individuals from the populations of a simulated metacommunity and
gives each one a trait value, so that intraspecific variation and local
adaptation can be studied with the T-statistics of Violle et al. (2012)
([`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)).
The trait is the environmental optimum of
[`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md),
the one that sets growth in
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).
The value of individual `k` of species `j` at site `i` is

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
  env = NULL
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

## Value

Data frame with one row per individual and the columns `site`,
`species`, `individual` (its number within the population), `abundance`
(of the population), `env` (of the site), `optimum` (of the species),
`population_mean` (the expected value of the population,
`theta_j + h_j (E_i - theta_j) + d_ij`, after character displacement)
and `trait`. Populations with no individual measured are absent.

## Details

`z_ijk = theta_j + h_j (E_i - theta_j) + d_ij + e_ijk`,

where `theta_j` is the optimum of the species, `E_i` the environment of
the site, `h_j` (`local_adaptation`) the fraction of the gap between
them that the population closes, `d_ij ~ N(0, population_sd^2)` a
deviation of the population that does not track the environment (drift,
founder effects), and `e_ijk ~ N(0, intraspecific_sd^2)` the variation
among individuals of a population. With `local_adaptation = 0` and
`population_sd = 0` all populations of a species share its mean; with
`local_adaptation = 1` every population sits on the environment of its
site and species no longer differ within a site.

With `displacement > 0`, the populations of each site also diverge from
one another, as in ecological character displacement. Their means are
the equilibrium of

`mu_ij = m_ij + c sum_k p_ik (mu_ij - mu_ik) exp(-(mu_ij - mu_ik)^2 / (2 w^2))`,

where `m_ij` is the mean above (`theta_j + h_j (E_i - theta_j) + d_ij`),
`c` is `displacement`, `p_ik` the share of species `k` among the
individuals of the site and `w` (`displacement_width`) the width of a
Gaussian competition kernel. Each population is pushed away from the
populations close to it in trait, more so from abundant ones, and held
near `m_ij`, so `c` sets the strength of competition relative to the
selection that keeps the population near `m_ij`. The push from one
competitor is largest at a distance `w` and at most `c w exp(-1/2)`,
about `0.6 c w`. The equilibrium is found by damped fixed-point
iteration from `m_ij`, with steps of `1 / (1 + c)`. Populations with the
same `m_ij`, as with `local_adaptation = 1` and no `population_sd`, are
not pushed apart, since the push is symmetric.

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

## See also

[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md) for
the T-statistics.

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
