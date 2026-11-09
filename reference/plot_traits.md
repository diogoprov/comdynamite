# Plot traits on the tree: phylomorphospace or traitgram

Draws simulated traits with the tree that generated them. With two
traits, a phylomorphospace: the tips and internal nodes in the plane of
the two traits, joined by the branches of the tree. With one trait, a
traitgram: the trait against time from the root. Branches are coloured
by selective regime, with the change of colour at the point of each
shift in a mapped history.

## Usage

``` r
plot_traits(
  x,
  tree = NULL,
  traits = c(1, 2),
  regimes = NULL,
  colours = NULL,
  tip_labels = FALSE,
  legend = TRUE,
  point_size = 1.8,
  linewidth = 0.7,
  ...,
  labels = NULL
)
```

## Arguments

- x:

  Traits of the tips from
  [`trait_ou()`](https://provetelab.org/comdynamite/reference/trait_ou.md):
  a matrix with one row per tip, named by tip label.

- tree:

  The tree of `x`, a `phylo` object; or `NULL` to draw only the tips.

- traits:

  The two columns of `x` to draw (names or numbers), or one for a
  traitgram.

- regimes:

  Regimes of the branches, as in
  [`trait_ou()`](https://provetelab.org/comdynamite/reference/trait_ou.md):
  a vector in the cladewise order of the edges (from
  [`paint_regimes()`](https://provetelab.org/comdynamite/reference/paint_regimes.md))
  or a tree with a mapped history (class `simmap`, from
  [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)).
  `NULL` draws every branch in one colour.

- colours:

  Colours of the regimes, named by regime or in the order of the sorted
  regimes; one colour when there are no regimes.

- tip_labels:

  Should the tips be labelled?

- legend:

  Draw the legend of the regimes?

- point_size:

  Size of the tips.

- linewidth:

  Width of the branches.

- labels, ...:

  Deprecated: `labels` is now `tip_labels`; graphical arguments of base
  graphics are no longer used (modify the ggplot object instead).

## Value

A ggplot object.

## Details

The values of the internal nodes are those kept by
[`trait_ou()`](https://provetelab.org/comdynamite/reference/trait_ou.md)
in the attribute `nodes`, so the branches show the simulated path
between nodes, not a reconstruction. Without that attribute, only the
tips are drawn.

## See also

[`trait_ou()`](https://provetelab.org/comdynamite/reference/trait_ou.md),
[`rate_matrix()`](https://provetelab.org/comdynamite/reference/rate_matrix.md)
and
[`plot.rate_matrix()`](https://provetelab.org/comdynamite/reference/plot.rate_matrix.md).

Other species traits:
[`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md),
[`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md),
[`food_web_rates()`](https://provetelab.org/comdynamite/reference/food_web_rates.md),
[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md),
[`paint_regimes()`](https://provetelab.org/comdynamite/reference/paint_regimes.md),
[`partner_turnover()`](https://provetelab.org/comdynamite/reference/partner_turnover.md),
[`plot.rate_matrix()`](https://provetelab.org/comdynamite/reference/plot.rate_matrix.md),
[`plot_rates()`](https://provetelab.org/comdynamite/reference/plot_rates.md),
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
set.seed(2)
tr <- treeshapesim::sim_tree_shape(30)
habitat <- trait_mk(tr, q = 0.5, states = c("forest", "open"))
r <- rate_matrix(c(1, 1), correlation = list(forest = 0.8, open = -0.8),
                 constraint = "variance")
x <- trait_ou(habitat, sigma2 = r, regimes = habitat)
plot_traits(x, habitat, regimes = habitat)

# traitgram of the first trait
plot_traits(x, habitat, traits = 1, regimes = habitat)
```
