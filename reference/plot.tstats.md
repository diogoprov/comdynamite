# Plot the null-model tests of each site

Draws the tests of
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md) or
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md)
site by site, in one panel per trait and statistic. With
`type = "intervals"` (the default) each site is a row with its null
distribution as two nested intervals (thin: `level`; thick: `inner`),
the null median as a tick and the observed value as a point coloured by
the result of the test: below the null, above it, or within it. Sites
that cannot be tested (such as `T_IP.IC` at a site with one species) are
drawn as hollow points. With `type = "ses"` the standardised effect
sizes are drawn instead, with the band of \\\pm z\_{1 - (1 -
\mathrm{level})/2}\\; sites that cannot be tested have no standardised
effect size and are left empty. The colours still give the result of the
permutation test, which does not assume a normal null distribution. With
`hatch = TRUE` (the default) the area expected by chance (the outer null
interval of each site, or the band of standardised effect sizes) is
hatched.

## Usage

``` r
# S3 method for class 'tstats'
plot(
  x,
  type = c("intervals", "ses"),
  level = 0.95,
  inner = 0.5,
  statistics = NULL,
  traits = NULL,
  sites = NULL,
  order = c("data", "ses", "observed"),
  colours = c("#0072B2", "grey55", "#D55E00"),
  point_size = 2,
  interval_width = 2.5,
  hatch = TRUE,
  hatch_density = 60,
  scales = "free_x",
  ...
)

# S3 method for class 'community_index'
plot(x, ...)
```

## Arguments

- x:

  Result of
  [`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
  or
  [`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md)
  with `nperm > 0`.

- type:

  `"intervals"` or `"ses"`.

- level, inner:

  Coverage of the outer and inner null intervals; `1 - level` is also
  the size of the two-sided test (as in
  [`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md)).

- statistics, traits, sites:

  Subsets to draw; `NULL` (the default) draws all.

- order:

  Order of the sites within each panel's axis: `"data"` (as in `x`),
  `"ses"` (by the mean standardised effect size over the panels drawn)
  or `"observed"` (by the mean observed value).

- colours:

  Three colours for sites below the null, within it and above it. The
  defaults are blue, grey and vermillion of the Okabe-Ito palette.

- point_size, interval_width:

  Size of the observed points and width of the thick inner interval.

- hatch:

  Hatch the area expected by chance; with `FALSE` it is drawn as a thin
  line (`"intervals"`) or a grey band (`"ses"`).

- hatch_density:

  Number of hatching lines across the width of a panel.

- scales:

  Passed to
  [`ggplot2::facet_grid()`](https://ggplot2.tidyverse.org/reference/facet_grid.html);
  `"free_x"` (the default) lets each statistic keep its own range.

- ...:

  Not used.

## Value

A ggplot object.

## Details

The plot is built from
[`null_intervals()`](https://provetelab.org/comdynamite/reference/null_intervals.md),
and the ggplot object can be changed further with ggplot2: for instance,
`+ ggplot2::facet_wrap(~ statistic)`, `+ ggplot2::theme_bw()` or
`+ ggplot2::scale_x_log10()`.

## See also

[`null_intervals()`](https://provetelab.org/comdynamite/reference/null_intervals.md),
[`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md),
[`plot.trait_flex()`](https://provetelab.org/comdynamite/reference/plot.trait_flex.md).

Other trait variation among individuals:
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
[`neighbour_metrics`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
[`nested_variance()`](https://provetelab.org/comdynamite/reference/nested_variance.md),
[`null_intervals()`](https://provetelab.org/comdynamite/reference/null_intervals.md),
[`plot.nested_variance()`](https://provetelab.org/comdynamite/reference/plot.nested_variance.md),
[`plot.trait_flex()`](https://provetelab.org/comdynamite/reference/plot.trait_flex.md),
[`plot.trait_rarefaction()`](https://provetelab.org/comdynamite/reference/plot.trait_rarefaction.md),
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
ind <- data.frame(site = rep(paste0("s", 1:8), each = 20),
                  species = sample(c("a", "b", "c"), 160, replace = TRUE))
ind$size <- c(a = 10, b = 20, c = 30)[ind$species] + rnorm(160, sd = 3)
res <- tstats(ind, traits = "size", nperm = 199)
plot(res)

plot(res, type = "ses", order = "ses")

# an index of communities
plot(community_index(ind, stats::var, traits = "size", null = "regional.ind", nperm = 199))
```
