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
  min_n = 1,
  founder_sd = 0,
  monopolization = 0,
  gene_flow = 1,
  source = c("traits", "evolution"),
  session_sd = 0
)
```

## Arguments

- x:

  An `mc_sim` object from
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
  (or an `mc_evolving` one from
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)).

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
  uses the abundances at the last step. Several steps simulate
  collections on several dates: each date is sampled on its own (new
  individuals, drawn from the populations then), and the result gets a
  column `step`. The deviations of each population from `population_sd`
  (and `founder_sd`) are kept across dates, so populations change
  between dates only through the simulation. For evolving
  metacommunities, the records keep only the lineages extant at the end,
  so an earlier date misses species that died out later.

- env:

  Optional environment of each site at the sampled step, or, with
  several steps, a matrix with one column per step (a vector is used for
  every date). The default is the environment of the landscape for every
  date; give it for simulations with `env_time`. With
  `source = "evolution"` the environment comes from the simulation.

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

- founder_sd:

  Standard deviation of the founder deviation of a population founded by
  a single individual (\\s_f\\ in the section Founder effects and
  monopolization); the deviation of a population founded by \\n\\
  individuals has standard deviation `founder_sd / sqrt(n)`. One value
  or one per species; 0 (the default) for none.

- monopolization:

  Rate per time step (\\\rho\\) at which a population closes the gap
  between the optimum of its species and the environment of its site
  after it arrives, so that early colonists end up better adapted than
  late ones. One value or one per species; 0 (the default) for none.

- gene_flow:

  Rate (\\\gamma\\) at which immigration erodes founder deviations and
  adaptation since arrival, per unit of migrant load: 1 (the default)
  when every newcomer mixes fully with the residents, 0 for no erosion.
  One value or one per species. Has no effect without `founder_sd` or
  `monopolization`.

- source:

  Where population means come from: `"traits"` (the default) builds them
  from the optimum of each species with the arguments above;
  `"evolution"` takes them from a simulation run with `evolution`
  ([`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md)),
  as the genetic mean plus the plastic part of each population at the
  sampled step, and draws individuals around them with the genetic plus
  environmental variance of each species (`genetic_var + env_var`, plus,
  for a species of several lineages, the variance among them in the
  site), so that the T-statistics of
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  reflect the variation that evolved. `intraspecific_sd`, when given,
  then adds variation such as measurement error, and `local_adaptation`,
  `population_sd`, `displacement`, `founder_sd`, `monopolization` and
  `env` must keep their defaults (`env` is the environment at the
  sampled step). At a `step` before the last, the simulation needs
  `record_steps` with that step, which keep the genetic means (and, for
  evolving metacommunities, of the lineages extant at the end). Evolving
  metacommunities from
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)
  are accepted too. With several evolving traits (one per niche axis),
  individuals are drawn from a multivariate normal with covariance \\G +
  E\\ of their species, and `intraspecific_sd` takes one value or one
  per trait; the output then has, for each trait, the columns
  `env_<trait>` and `mean_<trait>` and the trait itself (named after its
  axis) instead of `env`, `optimum`, `population_mean` and `trait`, with
  the trait names in the attribute `"traits"`. `session_sd` shifts each
  trait on its own.

- session_sd:

  Standard deviation of an error shared by all individuals measured on
  the same date (`step`), such as a difference between measurers or
  calibrations; it is added to `trait`, not to `population_mean`. 0 (the
  default) for none.

## Value

Data frame with one row per individual and the columns `step` (with
several steps), `site`, `species`, `individual` (its number within the
population), `abundance` (of the population), `env` (of the site),
`optimum` (of the species), `population_mean` (the expected value of the
population, \\m\_{ij}\\ above, after character displacement) and
`trait`, and with `founder_sd` or `monopolization` also `arrival` (the
step at which the population was founded), before `trait`. Populations
with no individual measured are absent.

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

## Founder effects and monopolization

`population_sd` and `local_adaptation` treat every population of a
species alike, whatever its history. Two further mechanisms make the
mean of a population depend on when and how it was founded, from the
colonisation history recorded by
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
(see
[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md)).
Both need the last step (`step = NULL`). For a population founded at
step \\a\_{ij}\\ by \\n\_{ij}\\ individuals, with migrant load
\\L\_{ij}\\ since then and \\T\\ steps in the simulation:

- **Founder effects** (`founder_sd`, \\s_f\\): the population starts
  from the mean of its founders, a deviation \\f\_{ij} \sim N(0, s_f^2 /
  n\_{ij})\\ from the species mean, so few founders deviate more.
  `founder_sd` is thus the deviation of a population founded by a single
  individual. The deviation persists, and gene flow erodes it toward the
  species mean: it is multiplied by \\g\_{ij} = \exp(-\gamma L\_{ij})\\,
  where \\\gamma\\ is `gene_flow`.

- **Monopolization** (`monopolization`, a rate \\\rho\\ per step): after
  arrival the population adapts to its site, closing a share \\u\_{ij}\\
  of the gap between the optimum of the species and the environment.
  Without gene flow, \\u\_{ij} = 1 - \exp(-\rho t\_{ij})\\, where
  \\t\_{ij} = T - a\_{ij}\\ is the age of the population. Newcomers pull
  it back toward the species mean, so \\du/dt = \rho (1 - u) - \gamma
  \bar m u\\, with \\\bar m = L\_{ij} / t\_{ij}\\ the mean migrant load
  per step since arrival, whose solution is \\u\_{ij} = \rho / (\rho +
  \gamma \bar m) \[1 - \exp(-\rho t\_{ij} - \gamma L\_{ij})\]\\:
  adaptation tends to a balance between selection and migration instead
  of reaching the environment. Early colonists have had longer to adapt
  and end up closer to the environment of their site than late arrivals,
  which is the evolution-mediated priority effect of the monopolization
  hypothesis (Urban & De Meester 2009; De Meester et al. 2016). The gap
  between early and late colonists is largest when \\1 / \rho\\ is close
  to the spread of arrival times; much faster adaptation brings every
  population to the balance and erases it. It combines with
  `local_adaptation` (\\h_j\\) so that the share of the gap closed is
  \\h_j + (1 - h_j) u\_{ij}\\.

The mean of the population is then \\\theta_j + \[h_j + (1 - h_j)
u\_{ij}\] (E_i - \theta_j) + d\_{ij} + g\_{ij} f\_{ij}\\, before
character displacement. Use `population_sd` for deviations unrelated to
history (drift at equilibrium, unmeasured microhabitat),
`local_adaptation` for adaptation that every population reaches whatever
its age, `founder_sd` for random deviations set at colonisation and kept
by isolation, and `monopolization` for adaptation that depends on
arrival order. The two arrival-dependent mechanisms leave different
signatures: founder deviations are random in direction and largest in
populations founded by few individuals and receiving few immigrants,
while monopolization moves populations toward the environment of their
site, more so the earlier they arrived. With `gene_flow = 0` neither is
eroded. Erosion of founder deviations multiplies them by the product of
the shares of residents over the steps, approximated by \\\exp(-L)\\;
monopolization assumes newcomers arrived at a constant rate since
arrival, \\\bar m\\, rather than following the steps in which they
actually came.

Limits of the two mechanisms:

- They need the last step (`step = NULL`): the colonisation history is
  recorded at the end of the simulation only, so a recorded `step` other
  than the last gives an error.

- They need arrivals spread over time. Without extirpation and with
  little dispersal, every population is founded at the start or during
  `seed_steps`, all have about the same age, and monopolization cannot
  separate early from late colonists. Use `extirpation > 0` and some
  `dispersal`, so that populations are lost and founded again throughout
  the simulation, and check the spread of arrival times with
  [`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md).

- The gap between early and late colonists does not grow without bound
  with `monopolization`: it is largest when \\1 / \rho\\ is close to the
  spread of arrival times, and much faster adaptation brings every
  population to the balance with migration and erases it. In a
  simulation of 200 steps, the difference in the share of the gap closed
  between colonists that arrived before and after the median rose from
  0.14 at \\\rho = 0.002\\ to 0.40 at 0.02 and fell to 0.11 at 0.2.
  Choose the range of `monopolization` (for instance, a prior in
  approximate Bayesian computation) from the arrival times of the
  simulation.

- Like the other trait deviations, they do not feed back on the
  dynamics: early colonists are better matched to their site but do not
  grow faster or compete better, so the demographic part of the
  monopolization hypothesis is not simulated.

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

Urban, M. C., & De Meester, L. (2009). Community monopolization: local
adaptation enhances priority effects in an evolving metacommunity.
*Proceedings of the Royal Society B: Biological Sciences*, 276,
4129–4138.
[doi:10.1098/rspb.2009.1382](https://doi.org/10.1098/rspb.2009.1382)

De Meester, L., Vanoverbeke, J., Kilsdonk, L. J., & Urban, M. C. (2016).
Evolving perspectives on monopolization and priority effects. *Trends in
Ecology & Evolution*, 31, 136–146.
[doi:10.1016/j.tree.2015.12.009](https://doi.org/10.1016/j.tree.2015.12.009)

## See also

[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md)
for the arrival history;
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
[`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md),
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
# arrival-dependent deviations need colonisation after the start: extirpation and dispersal
sim2 <- simulate_metacommunity(landscape_implicit(20), species_traits(tree),
                               competition_matrix(tree, "phylo"), n_steps = 200,
                               dispersal = 0.005, extirpation = 0.01)
# founder effects: deviations set at colonisation, larger with few founders, eroded by
# immigration
founders <- sample_individuals(sim2, founder_sd = 2, intraspecific_sd = 0.5)
# monopolization: early colonists close more of the gap to the environment of their site
mono <- sample_individuals(sim2, monopolization = 0.02, intraspecific_sd = 0)
closed <- (mono$population_mean - mono$optimum) / (mono$env - mono$optimum)
early <- mono$arrival <= stats::median(mono$arrival)
c(early = mean(closed[early]), late = mean(closed[!early]))
#>     early      late 
#> 0.7385448 0.3317174 
```
