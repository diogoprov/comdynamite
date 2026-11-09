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

## Examples

``` r
set.seed(1)
ind <- data.frame(site = rep(paste0("s", 1:5), each = 40),
                  species = sample(c("a", "b", "c"), 200, replace = TRUE))
ind$trait <- c(a = 10, b = 12, c = 15)[ind$species] + rnorm(200)
rare <- trait_rarefaction(ind, n = c(2, 4, 8), reps = 10)
rare
#>   effort statistic      mean        low      high reps
#> 1      2   T_IP.IC 0.1725745 0.08750932 0.2771844   10
#> 2      2   T_IC.IR 1.1399778 1.11425692 1.1519776   10
#> 3      2   T_PC.PR 1.3735930 1.34136002 1.3884881   10
#> 4      4   T_IP.IC 0.1783959 0.14080405 0.2292743   10
#> 5      4   T_IC.IR 1.0653356 1.05554298 1.0717343   10
#> 6      4   T_PC.PR 1.3888587 1.37405313 1.3984714   10
#> 7      8   T_IP.IC 0.1867279 0.16365975 0.2042643   10
#> 8      8   T_IC.IR 1.0328705 1.03099574 1.0340827   10
#> 9      8   T_PC.PR 1.3968893 1.39384630 1.3988592   10
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
