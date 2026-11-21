# Correlation between the axes of a niche

Attaches to a list of traits (one set per environmental axis, as for
[`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md)
and
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md))
a correlation between the responses of each species to the axes, so that
the niche is a multivariate Gaussian with breadths `sigma` on each axis
and this correlation. A species with a positive correlation between two
axes grows best where both depart from its optima in the same direction,
as a species that tolerates warm water only if it is also well
oxygenated.

## Usage

``` r
niche_correlation(traits, cor)
```

## Arguments

- traits:

  List of data frames from
  [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md),
  one per axis, with the same species.

- cor:

  A correlation matrix (axes x axes), shared by all species, or a list
  of them named by species. Each must be symmetric, with ones on the
  diagonal and positive definite.

## Value

`traits`, with the attribute `niche_cor`.

## Details

With `evolution` in
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
and in
[`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md),
the correlation enters the matrix of Gaussian selection, \\\Omega = D R
D\\ with \\D\\ the diagonal matrix of the breadths (see the section
Several traits of
[`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md)):
growth, the response to selection, the selection of immigrants and the
overlap of niches (`niche_overlap`, and
`competition_matrix(type = "trait")`) all use it. New species take the
correlation of their mother. An axis with an infinite breadth exerts no
selection, whatever its correlation; a breadth of 0 cannot be combined
with a correlation.

## See also

Other species traits:
[`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md),
[`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md),
[`food_web_rates()`](https://provetelab.org/comdynamite/reference/food_web_rates.md),
[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md),
[`paint_regimes()`](https://provetelab.org/comdynamite/reference/paint_regimes.md),
[`partner_turnover()`](https://provetelab.org/comdynamite/reference/partner_turnover.md),
[`plot.rate_matrix()`](https://provetelab.org/comdynamite/reference/plot.rate_matrix.md),
[`plot_rates()`](https://provetelab.org/comdynamite/reference/plot_rates.md),
[`plot_traits()`](https://provetelab.org/comdynamite/reference/plot_traits.md),
[`rate_matrix()`](https://provetelab.org/comdynamite/reference/rate_matrix.md),
[`region_pool()`](https://provetelab.org/comdynamite/reference/region_pool.md),
[`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md),
[`tip_history_metrics()`](https://provetelab.org/comdynamite/reference/tip_history_metrics.md),
[`trait_bm()`](https://provetelab.org/comdynamite/reference/trait_bm.md),
[`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md),
[`trait_ou()`](https://provetelab.org/comdynamite/reference/trait_ou.md),
[`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md),
[`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md)

## Examples

``` r
set.seed(1)
tr <- list(temperature = species_traits(2, optima = c(40, 60), sigma = c(10, 10),
                                        r_max = c(5, 5), rescale_optima = FALSE),
           oxygen = species_traits(2, optima = c(40, 60), sigma = c(10, 10),
                                   r_max = c(1, 1), rescale_optima = FALSE))
env <- cbind(temperature = c(30, 50), oxygen = c(30, 50))
growth_rates(env, tr)
#>           sp1         sp2
#> [1,] 1.839397 0.000617049
#> [2,] 1.839397 1.839397206
growth_rates(env, niche_correlation(tr, matrix(c(1, 0.8, 0.8, 1), 2)))
#>           sp1        sp2
#> [1,] 2.868767 0.03368973
#> [2,] 2.868767 2.86876710
```
