# Turnover of interaction partners in a simulated metacommunity

Counts, for each species, how often the composition around it changes:
the number of other species that enter or leave the sites where the
species is present, between successive recorded steps of
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).
The count is divided by the number of sites and time steps over which
the species was present at both ends of an interval, so it is a rate of
change of partners per site and per time step. It measures the entries
and exits of partners of scenario 5 of the notes on branch rates, and
feeds
[`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md).

## Usage

``` r
partner_turnover(x, partners = c("all", "competitors"))
```

## Arguments

- x:

  An `mc_sim` object with `history`, that is, run with `record_steps`.

- partners:

  `"all"` counts every other species as a partner; `"competitors"`
  counts only the species with a non-zero competition coefficient with
  the focal species in `x$alpha` (in either direction).

## Value

Named numeric vector, one value per species. Species never present at
two successive recorded steps get zero.

## See also

[`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md),
[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md)

Other species traits:
[`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md),
[`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md),
[`food_web_rates()`](https://provetelab.org/comdynamite/reference/food_web_rates.md),
[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md),
[`paint_regimes()`](https://provetelab.org/comdynamite/reference/paint_regimes.md),
[`plot.rate_matrix()`](https://provetelab.org/comdynamite/reference/plot.rate_matrix.md),
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
set.seed(1)
tr <- treeshapesim::sim_tree_shape(10)
sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tr),
                              competition_matrix(tr, "none"), n_steps = 300,
                              record_steps = seq(100, 300, 20), extirpation = 0.02)
partner_turnover(sim)
#>         t1         t2         t3         t4         t5         t6         t7 
#> 0.05564516 0.08333333 0.05588235 0.09000000 0.04500000 0.07142857 0.00000000 
#>         t8         t9        t10 
#> 0.00000000 0.06093750 0.06106870 
```
