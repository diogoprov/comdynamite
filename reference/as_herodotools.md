# Simulations as input for Herodotools, with the known history

Herodotools (Nakamura et al. 2024) estimates the biogeographic history
of assemblages from a tree, the occurrence of species in sites, the
region of each site and the area of each ancestor. `as_herodotools()`
gives these inputs from a simulation or an `mc_data` object, with the
names of the arguments of `Herodotools::calc_age_arrival()`,
`Herodotools::calc_dispersal_from()` and
`Herodotools::calc_insitu_diversification()`, so the result can go to
them through [`do.call()`](https://rdrr.io/r/base/do.call.html). When
the areas evolved as a discrete trait from
[`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md),
whose states are the regions, `history` gives the true area of each
ancestor; otherwise reconstruct them elsewhere, for example with
BioGeoBEARS.
[`region_pool()`](https://provetelab.org/comdynamite/reference/region_pool.md)
limits each species to the sites of its region in the simulation, and
[`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md)
gives the true ages at which lineages arrived in each region, to compare
with what Herodotools estimates from the nodes.

## Usage

``` r
as_herodotools(x, regions, history = NULL, tree = NULL)
```

## Arguments

- x:

  An `mc_sim` or `mc_data` object, or a community matrix with sites in
  rows and species in columns, named by species.

- regions:

  Region of each site, by site name or in the order of the sites.

- history:

  Optional `simmap` tree from
  [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)
  whose states are the regions, with the same species as `x`.

- tree:

  The tree, when `x` has none and `history` is not given.

## Value

A list with `W`, the incidence matrix of sites x species; `tree`, a
plain `phylo` object; `biogeo`, a one-column data frame with the region
of each site; and, with `history`, `ancestral.area`, a one-column data
frame with the area of each internal node, in the order of their numbers
in `tree$edge` (the root first), as Herodotools reads it.

## Details

Herodotools finds a region among the areas of the ancestors by matching
its name as a pattern, so no region name may be contained in another,
such as `"A"` and `"AB"`. Sites without species are dropped with a
warning, since Herodotools needs at least one species in each site.

## References

Nakamura, G., Rodrigues, A. V., Luza, A. L., Maestri, R., Debastiani, V.
J., & Duarte, L. D. S. (2024). Herodotools: An R package to integrate
macroevolution, biogeography and community ecology. *Journal of
Biogeography*, 51, 533–543.
[doi:10.1111/jbi.14610](https://doi.org/10.1111/jbi.14610)

## See also

[`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md),
[`region_pool()`](https://provetelab.org/comdynamite/reference/region_pool.md),
[`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md),
[`tip_history_metrics()`](https://provetelab.org/comdynamite/reference/tip_history_metrics.md)

Other species traits:
[`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md),
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
[`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md),
[`traits_ou()`](https://provetelab.org/comdynamite/reference/traits_ou.md),
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
