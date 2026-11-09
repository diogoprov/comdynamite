# Plot rates that vary among branches

Draws the tree as a rectangular phylogram with each branch coloured by
its rate. Works with any rate per branch: the truth from
[`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md),
[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md),
[`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md)
or
[`food_web_rates()`](https://provetelab.org/comdynamite/reference/food_web_rates.md),
or estimates such as those of RRphylo, as long as they follow the order
of the rows of `tree$edge`.

## Usage

``` r
plot_rates(
  tree,
  rates = NULL,
  log = TRUE,
  palette = "viridis",
  tip_labels = FALSE,
  linewidth = 1,
  title = NULL
)
```

## Arguments

- tree:

  A `phylo` object with branch lengths, or a result of
  [`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md),
  which carries its tree and rates.

- rates:

  Rate of each branch, in the order of the rows of `tree$edge`, or a
  result of
  [`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md),
  whose `edge_rate` is used. Left out when `tree` is a result of
  [`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md).

- log:

  Colour by the logarithm of the rate (the default), since rates are
  multiplicative. Needs positive rates; use `FALSE` for estimates that
  can be negative.

- palette:

  Name of a palette of
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html).

- tip_labels:

  Draw the names of the tips.

- linewidth:

  Width of the branches.

- title:

  Title of the plot, or `NULL`.

## Value

A ggplot object.

## See also

[`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md),
[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md),
[`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md).

Other species traits:
[`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md),
[`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md),
[`food_web_rates()`](https://provetelab.org/comdynamite/reference/food_web_rates.md),
[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md),
[`paint_regimes()`](https://provetelab.org/comdynamite/reference/paint_regimes.md),
[`partner_turnover()`](https://provetelab.org/comdynamite/reference/partner_turnover.md),
[`plot.rate_matrix()`](https://provetelab.org/comdynamite/reference/plot.rate_matrix.md),
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
tree <- treeshapesim::sim_tree_shape(30)
x <- trait_rates(tree, "shift", shifts = 2, shift_size = c(8, 1 / 8))
plot_rates(x)

habitat <- trait_mk(tree, q = 0.5, states = c("forest", "open"))
plot_rates(tree, history_rates(habitat, "state", state_rates = c(forest = 1, open = 5)))
```
