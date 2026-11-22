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
  two-sided test that defines `deviation`. One number between 0 and 1,
  both excluded.

- inner:

  Coverage of the inner null interval. One number between 0 and 1,
  smaller than `level`.

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
[`drift_test()`](https://provetelab.org/comdynamite/reference/drift_test.md),
[`neighbour_metrics`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
[`nested_variance()`](https://provetelab.org/comdynamite/reference/nested_variance.md),
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
ind <- data.frame(site = rep(paste0("s", 1:6), each = 20),
                  species = sample(c("a", "b", "c"), 120, replace = TRUE))
ind$trait <- c(a = 10, b = 20, c = 30)[ind$species] + rnorm(120)
head(null_intervals(tstats(ind, n_perm = 99)))
#> tstats() implements the method of:
#>   Violle, C., Enquist, B. J., McGill, B. J., Jiang, L., Albert, C. H., Hulshof, C., Jung, V., & Messier, J. (2012). The return of the variance: intraspecific variability in community ecology. Trends in Ecology & Evolution, 27, 244-252. https://doi.org/10.1016/j.tree.2011.11.014
#>   Taudiere, A., & Violle, C. (2016). cati: an R package using functional traits to detect and quantify multi-level community assembly processes. Ecography, 39, 699-708. https://doi.org/10.1111/ecog.01433
#> Please cite these papers together with comdynamite (citation("comdynamite")).
#> This message is shown once per session; options(comdynamite.cite = FALSE) turns it off.
#>   site trait statistic    observed null_median  null_low null_high inner_low
#> 1   s1 trait   T_IP.IC 0.009373208    1.042140 0.7625785  1.125826 0.9789126
#> 2   s2 trait   T_IP.IC 0.022171031    1.056030 0.7206097  1.280613 0.9467706
#> 3   s3 trait   T_IP.IC 0.019857632    1.056567 0.6433469  1.442143 0.8319205
#> 4   s4 trait   T_IP.IC 0.012470624    1.039040 0.6851734  1.216771 0.9302977
#> 5   s5 trait   T_IP.IC 0.016391403    1.040285 0.8143519  1.148197 0.9521490
#> 6   s6 trait   T_IP.IC 0.014484968    1.020742 0.5792220  1.176768 0.9433539
#>   inner_high       ses p_low p_high deviation
#> 1   1.074819 -8.780631  0.01      1       low
#> 2   1.139538 -7.149138  0.01      1       low
#> 3   1.161632 -4.517790  0.01      1       low
#> 4   1.082762 -7.078147  0.01      1       low
#> 5   1.086457 -9.734629  0.01      1       low
#> 6   1.079700 -6.916428  0.01      1       low
```
