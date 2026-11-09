# Null intervals of the tests of each site

Summarises, for each site, trait and statistic, the null distribution
kept by
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md) or
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md):
its median, an outer interval (`level`) and an inner one (`inner`), with
the observed value, the standardised effect size and the one-sided tests
of
[`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md).
This is the table drawn by the
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods of
these objects, and can be used to build other figures.

## Usage

``` r
null_intervals(x, level = 0.95, inner = 0.5)
```

## Arguments

- x:

  Result of
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  or
  [`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md)
  with `n_perm > 0`.

- level:

  Coverage of the outer null interval, and `1 - level` the size of the
  two-sided test that defines `deviation`.

- inner:

  Coverage of the inner null interval.

## Value

A data frame with one row per site, trait and statistic (for
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
`statistic` is the name of the index) and the columns `observed`,
`null_median`, `null_low`, `null_high`, `inner_low`, `inner_high`,
`ses`, `p_low`, `p_high` and `deviation` (`"low"`, `"high"`, `"none"` or
`NA` for sites that cannot be tested).

## See also

[`plot.tstats()`](https://provetelab.org/comdynamite/reference/plot.tstats.md),
[`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md).

Other trait variation among individuals:
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
[`neighbour_metrics`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
[`nested_variance()`](https://provetelab.org/comdynamite/reference/nested_variance.md),
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
[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
[`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md)

## Examples

``` r
set.seed(1)
ind <- data.frame(site = rep(paste0("s", 1:6), each = 20),
                  species = sample(c("a", "b", "c"), 120, replace = TRUE))
ind$trait <- c(a = 10, b = 20, c = 30)[ind$species] + rnorm(120)
head(null_intervals(tstats(ind, n_perm = 99)))
#>   site trait statistic    observed null_median  null_low null_high inner_low
#> 1   s1 trait   T_IP.IC 0.009373208    1.045136 0.6343477  1.132411 0.9543454
#> 2   s2 trait   T_IP.IC 0.022171031    1.034900 0.6993979  1.218829 0.9206989
#> 3   s3 trait   T_IP.IC 0.019857632    1.032582 0.5455449  1.537149 0.7933045
#> 4   s4 trait   T_IP.IC 0.012470624    1.043910 0.6006849  1.199357 0.9239751
#> 5   s5 trait   T_IP.IC 0.016391403    1.026634 0.6838355  1.131407 0.9220763
#> 6   s6 trait   T_IP.IC 0.014484968    1.016252 0.7254197  1.137528 0.9367842
#>   inner_high       ses p_low p_high deviation
#> 1   1.085213 -7.629985  0.01      1       low
#> 2   1.119922 -7.228212  0.01      1       low
#> 3   1.171645 -3.853109  0.01      1       low
#> 4   1.084610 -6.616754  0.01      1       low
#> 5   1.077998 -7.816353  0.01      1       low
#> 6   1.073531 -8.396158  0.01      1       low
```
