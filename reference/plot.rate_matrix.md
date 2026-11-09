# Plot evolutionary rate matrices as ellipses

Draws each rate matrix as the ellipse of the changes of two traits per
unit time (one standard deviation by default), one ellipse per regime.
The shape shows the evolutionary correlation, the size the rates:
proportional matrices give nested ellipses of the same shape, shared
eigenvectors give ellipses with the same axes, equal variances give
ellipses inscribed in the same box.

## Usage

``` r
# S3 method for class 'rate_matrix'
plot(x, traits = c(1, 2), level = 1, colours = NULL, legend = "topleft", ...)
```

## Arguments

- x:

  A `rate_matrix` object from
  [`rate_matrix()`](https://provetelab.org/comdynamite/reference/rate_matrix.md).

- traits:

  The two traits to draw (names or numbers); with more traits, the
  ellipse is that of the marginal covariance of the two.

- level:

  Number of standard deviations of the ellipse.

- colours:

  Colours of the regimes, in the order of the regimes.

- legend:

  Where to put the legend of the regimes, or `NULL` for none.

- ...:

  Further arguments to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`NULL`, invisibly; called for the plot.

## See also

[`rate_matrix()`](https://provetelab.org/comdynamite/reference/rate_matrix.md),
[`plot_traits()`](https://provetelab.org/comdynamite/reference/plot_traits.md).

Other species traits:
[`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md),
[`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md),
[`food_web_rates()`](https://provetelab.org/comdynamite/reference/food_web_rates.md),
[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md),
[`paint_regimes()`](https://provetelab.org/comdynamite/reference/paint_regimes.md),
[`partner_turnover()`](https://provetelab.org/comdynamite/reference/partner_turnover.md),
[`plot_rates()`](https://provetelab.org/comdynamite/reference/plot_rates.md),
[`plot_traits()`](https://provetelab.org/comdynamite/reference/plot_traits.md),
[`rate_matrix()`](https://provetelab.org/comdynamite/reference/rate_matrix.md),
[`region_pool()`](https://provetelab.org/comdynamite/reference/region_pool.md),
[`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md),
[`tip_history_metrics()`](https://provetelab.org/comdynamite/reference/tip_history_metrics.md),
[`trait_bm()`](https://provetelab.org/comdynamite/reference/trait_bm.md),
[`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md),
[`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md),
[`traits_ou()`](https://provetelab.org/comdynamite/reference/traits_ou.md),
[`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md)

## Examples

``` r
r <- rate_matrix(c(1, 2), correlation = 0.6, constraint = "shared",
                 eigenvalues = list(slow = c(1, 0.2), fast = c(1, 2)))
plot(r)
```
