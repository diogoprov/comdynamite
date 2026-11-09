# Trait evolution by Brownian motion with Pagel's lambda

Draws from a multivariate normal distribution with covariance \\\lambda
C + (1 - \lambda)\\ \mathrm{diag}(C)\\, where \\C\\ is the phylogenetic
covariance matrix of `tree`.

## Usage

``` r
trait_bm(tree, lambda = 1, nodes = FALSE)
```

## Arguments

- tree:

  A `phylo` tree with branch lengths.

- lambda:

  Pagel's lambda, in \[0, 1\].

- nodes:

  Also draw the values at the internal nodes? They are drawn after the
  tips, given them, so the tips are the same as with `nodes = FALSE` and
  the same seed. Under lambda, the nodes are those of the tree whose
  internal branches are multiplied by lambda and whose terminal branches
  are lengthened to keep the tips where they were.

## Value

Named numeric vector with the simulated trait value of each tip. With
`nodes = TRUE`, the attribute `nodes` holds the values at the internal
nodes, named by their numbers in `tree$edge`, the root first.

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
[`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md),
[`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md),
[`traits_ou()`](https://provetelab.org/comdynamite/reference/traits_ou.md),
[`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md)

## Examples

``` r
set.seed(1)
trait_bm(treeshapesim::sim_tree_shape(5), lambda = 0.5)
#>           t1           t2           t3           t4           t5 
#> -0.005767173  2.275670676  0.763593461 -0.554890164 -1.062009857 
```
