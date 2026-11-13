# Sites where each species can occur, from its region

Gives the `pool` of
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md)
from a history of areas: each species can occur only in the sites of the
region (or regions, in a range such as `"A+B"`) it occupies at the tips
of `history`.

## Usage

``` r
region_pool(history, regions)
```

## Arguments

- history:

  A `simmap` tree from
  [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)
  whose states are the regions.

- regions:

  Region of each site. Names, when given, name the sites; otherwise they
  are `site1`, `site2` and so on, as in
  [`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).

## Value

Logical matrix of sites x species.

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
tree <- treeshapesim::sim_tree_shape(15)
areas <- trait_mk(tree, q = 0.5, states = c("A", "B"))
regions <- rep(c("A", "B"), each = 10)
sim <- simulate_metacommunity(landscape_implicit(20), species_traits(tree),
                              competition_matrix(tree, "none"), n_steps = 200,
                              pool = region_pool(areas, regions))
input <- as_herodotools(sim, regions, history = areas)
#> Warning: sites without species dropped: site1
str(input, max.level = 1)
#> List of 4
#>  $ W             : num [1:19, 1:15] 0 0 1 0 0 0 0 1 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ tree          :List of 4
#>   ..- attr(*, "class")= chr "phylo"
#>  $ biogeo        :'data.frame':  19 obs. of  1 variable:
#>  $ ancestral.area:'data.frame':  14 obs. of  1 variable:
arrival_ages(areas, sim, regions)[c(1, 11), ]
#>        t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12 t13 t14 t15
#> site1  NA NA NA NA NA NA NA NA NA  NA  NA  NA  NA  NA  NA
#> site11 NA NA NA NA NA  1 NA NA  1  NA  NA   1   1   1  NA
# with Herodotools installed: do.call(Herodotools::calc_age_arrival, input)
```
