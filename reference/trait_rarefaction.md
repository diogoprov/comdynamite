# Rarefaction of a statistic over sampling effort

Computes a statistic of the individuals, such as the mean T-statistics
of [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
or the decomposition of
[`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md),
on repeated subsamples of decreasing effort, with the schemes of
[`subsample_individuals()`](https://provetelab.org/comdynamite/reference/subsample_individuals.md).
By default effort is the number of individuals per population: `n`
individuals are drawn from each population (populations smaller than `n`
keep all their individuals), for a series of `n`. With `cover`, effort
is instead the share of the abundance of each site covered by the
species measured, from the most abundant down; all individuals of those
species are kept, or `n` of them when a single `n` is given. The curve
shows how much effort is needed before the statistic stops changing.
This rarefies individuals within populations or species within sites;
the rarefaction of species richness (for instance
[`vegan::rarefy()`](https://vegandevs.github.io/vegan/reference/rarefy.html)
or the package iNEXT) answers another question.

## Usage

``` r
trait_rarefaction(
  individuals,
  statistic = NULL,
  n = NULL,
  cover = NULL,
  reps = 20,
  level = 0.95,
  min_n = 2,
  site = "site",
  species = "species",
  abundance = NULL
)
```

## Arguments

- individuals:

  Data frame with one row per individual.

- statistic:

  Function of a data frame of individuals that returns a named numeric
  vector. The default gives the mean over sites of the three
  T-statistics of
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  for the column `trait`.

- n:

  Numbers of individuals per population; the default runs from 2 to the
  median size of the populations. With `cover`, a single number of
  individuals per population, or `NULL` (the default) to keep all of
  them.

- cover:

  Shares (0 to 1) of abundance covered by the species measured, for a
  rarefaction over species instead of individuals; `NULL` (the default)
  for a rarefaction over `n`.

- reps:

  Number of subsamples for each level of effort. With `cover` and no `n`
  the subsamples are all equal, and a single one is computed.

- level:

  Coverage of the interval of the statistic over subsamples.

- min_n, site, species:

  Passed to
  [`subsample_individuals()`](https://provetelab.org/comdynamite/reference/subsample_individuals.md).
  With `min_n = 2` (the default) populations of one individual are left
  out, as in
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  analyses.

- abundance:

  Name of a column with the abundance of each population, used with
  `cover` (passed to
  [`subsample_individuals()`](https://provetelab.org/comdynamite/reference/subsample_individuals.md)).

## Value

A data frame of class `trait_rarefaction` with the columns `effort` (`n`
or `cover`, as the attribute `by` says), `statistic` (the names returned
by `statistic`), `mean`, `low` and `high` (quantiles over subsamples)
and `reps` (subsamples with a finite value), with a
[plot()](https://provetelab.org/comdynamite/reference/plot.trait_rarefaction.md)
method. The values of all subsamples are in the attribute `values`.

## See also

[`subsample_individuals()`](https://provetelab.org/comdynamite/reference/subsample_individuals.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md).

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
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md),
[`subsample_individuals()`](https://provetelab.org/comdynamite/reference/subsample_individuals.md),
[`trait_densities()`](https://provetelab.org/comdynamite/reference/trait_densities.md),
[`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md),
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md),
[`trait_gradient_brms()`](https://provetelab.org/comdynamite/reference/trait_gradient_brms.md),
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md),
[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md),
[`trait_gradient_weights()`](https://provetelab.org/comdynamite/reference/trait_gradient_weights.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
[`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md)

## Examples

``` r
set.seed(1)
ind <- data.frame(site = rep(paste0("s", 1:5), each = 40),
                  species = sample(c("a", "b", "c"), 200, replace = TRUE))
ind$trait <- c(a = 10, b = 12, c = 15)[ind$species] + rnorm(200)
rare <- trait_rarefaction(ind, n = c(2, 4, 8), reps = 10)
rare
#>   effort statistic      mean       low      high reps
#> 1      2   T_IP.IC 0.1792578 0.1114679 0.2825401   10
#> 2      2   T_IC.IR 1.1253155 1.0673701 1.1558713   10
#> 3      2   T_PC.PR 1.3538095 1.2787934 1.3944207   10
#> 4      4   T_IP.IC 0.1761934 0.1233734 0.2363405   10
#> 5      4   T_IC.IR 1.0593501 1.0317182 1.0713293   10
#> 6      4   T_PC.PR 1.3799166 1.3390941 1.3978650   10
#> 7      8   T_IP.IC 0.1863448 0.1504747 0.2177141   10
#> 8      8   T_IC.IR 1.0330574 1.0312667 1.0339801   10
#> 9      8   T_PC.PR 1.3972044 1.3943099 1.3987061   10
# over the share of abundance covered by the species measured
trait_rarefaction(ind, cover = c(0.5, 0.8, 1))
#>   effort statistic      mean       low      high reps
#> 1    0.5   T_IP.IC 0.3638244 0.3638244 0.3638244    1
#> 2    0.5   T_IC.IR 0.7635113 0.7635113 0.7635113    1
#> 3    0.5   T_PC.PR 1.2750472 1.2750472 1.2750472    1
#> 4    0.8   T_IP.IC 0.1958225 0.1958225 0.1958225    1
#> 5    0.8   T_IC.IR 1.0114336 1.0114336 1.0114336    1
#> 6    0.8   T_PC.PR 1.3996447 1.3996447 1.3996447    1
#> 7    1.0   T_IP.IC 0.1958225 0.1958225 0.1958225    1
#> 8    1.0   T_IC.IR 1.0114336 1.0114336 1.0114336    1
#> 9    1.0   T_PC.PR 1.3996447 1.3996447 1.3996447    1
```
