# True ages at which lineages arrived in each region

For each species in each site, the time before the present at which its
lineage last entered the region of the site, from the history of areas:
the time spent in that region along the path from the tip towards the
root, until the first period in another region, plus the time from the
tip to the present when the tree is not ultrametric. A lineage that
stayed in the region since the root gets the age of the root. With
ranges of several regions (`"A+B"`, from
[`range_history()`](https://provetelab.org/comdynamite/reference/range_history.md)),
the age in each region is the time since the range last came to include
it. For an evolving metacommunity with protracted speciation, the tips
of
[`range_history()`](https://provetelab.org/comdynamite/reference/range_history.md)
have the range of the lineage that represents each species, so cells
occupied only by its incipient lineages, elsewhere, get `NA`. These are
the true values of what `Herodotools::calc_age_arrival()` estimates from
the areas of the nodes, which can only see arrivals at nodes and so
gives younger ages.

## Usage

``` r
arrival_ages(history, x, regions)
```

## Arguments

- history:

  A `simmap` tree from
  [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)
  whose states are the regions.

- x:

  An `mc_sim`, `mc_evolving` or `mc_data` object, or a community matrix
  with sites in rows and species in columns, named by species.

- regions:

  Region of each site, by site name or in the order of the sites.

## Value

Matrix of sites x species with the age of arrival, `NA` where the
species is absent or where its tip is in another region.

## See also

Other species traits:
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
