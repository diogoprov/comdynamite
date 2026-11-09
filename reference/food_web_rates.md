# Rates of trait evolution from the trophic level of each lineage in a food web

Turns the rates of change of the mean traits of the resource, the
consumer and the predator in
[`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md)
into rates of evolution of a continuous trait on the branches of a tree
whose tips belong to these levels, for `trait_rates(edge_rates = )`. It
implements the expectation that lineages at different trophic positions
evolve at different speeds, and that the rates change when a level is
added (scenario 4 of the notes on branch rates).

## Usage

``` r
food_web_rates(x, tree, levels, sigma2 = 1, k = 1, window = NULL)
```

## Arguments

- x:

  Result of
  [`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md).

- tree:

  A tree of class `phylo`.

- levels:

  Trophic level of each tip, `"R"`, `"N"` or `"P"`, named by tip label
  or in the order of `tree$tip.label`.

- sigma2:

  Rate outside the stretch informed by the metacommunity.

- k:

  Increase of the rate, relative to `sigma2`, per unit of speed (trait
  units per unit of time of the food web).

- window:

  Length of the end of each tip branch represented by the metacommunity,
  in the units of the branch lengths; `NULL` uses the whole tip branch.
  A window longer than a tip branch is cut to its length.

## Value

Numeric vector with the rate of each branch, in the order of the rows of
`tree$edge`, with the speed of each level as attribute `"speed"`.

## Details

The speed of a level is the mean of its rate of trait change, \\\|d\bar
z / dt\|\\, over the steps at which it is present. Like
[`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md),
the food web runs in ecological time and only informs the tip branches:
over the last `window` units of time of each tip branch the rate is
\\\mathrm{sigma2}\\ (1 + k \cdot \mathrm{speed})\\ of the level of the
tip, and elsewhere it is `sigma2`. These speeds are directional changes
of the mean trait under selection, not the variance per unit time of
Brownian motion, so they give relative rates among levels, scaled by
`k`, not an equivalence between the two.

## References

McPeek, M. A. (2017). *Evolutionary Community Ecology*. Princeton
University Press.
[doi:10.23943/princeton/9780691088778.001.0001](https://doi.org/10.23943/princeton/9780691088778.001.0001)

## See also

[`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md),
[`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md),
[`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md)

Other species traits:
[`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md),
[`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md),
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
web <- food_web_dynamics(seq(0, 4000, by = 20), init = c(R = 50, N = 0.1, P = 0.1),
                         traits = c(R = 10, N = 4, P = 6), introductions = c(P = 2000),
                         dt = 0.5)
tr <- treeshapesim::sim_tree_shape(9)
levels <- rep(c("R", "N", "P"), each = 3)
# speeds are small (trait units per time step), so k is large
rates <- food_web_rates(web, tr, levels, k = 500)
attr(rates, "speed")
#>           R           N           P 
#> 0.001989566 0.002973562 0.000143305 
x <- trait_rates(tr, edge_rates = rates)
```
