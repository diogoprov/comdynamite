# Assign selective regimes to clades

Paints the branches of a tree with selective regimes for
[`traits_ou()`](https://provetelab.org/comdynamite/reference/traits_ou.md):
every branch starts in regime `"0"`, and the branches descending from
each node in `nodes` take a new regime (`"1"`, `"2"`, ...). Nodes are
painted in order, so a later node nested in an earlier one overrides it.

## Usage

``` r
paint_regimes(tree, nodes, stem = FALSE)
```

## Arguments

- tree:

  A rooted `phylo` tree.

- nodes:

  Internal node numbers whose descendant branches get a new regime.

- stem:

  Should the branch leading to each node also take the new regime (the
  shift happens at the start of the stem rather than at the node)?

## Value

Character vector with the regime of each branch, in the order of
`tree$edge` after reordering it "cladewise" (the order
[`traits_ou()`](https://provetelab.org/comdynamite/reference/traits_ou.md)
uses).

## See also

Other species traits:
[`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md),
[`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md),
[`food_web_rates()`](https://provetelab.org/comdynamite/reference/food_web_rates.md),
[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md),
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
tr <- treeshapesim::sim_tree_shape(8)
table(paint_regimes(tr, nodes = 10))
#> 
#> 0 1 
#> 6 8 
```
