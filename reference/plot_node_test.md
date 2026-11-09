# Plot the results of a node test on the phylogeny

Draws the tree with the associations found by
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
at each node, as in Figure 3 of Leibold et al. (2010). Blue circles are
predictors and red triangles are biogeographic events, each with the
names of the variables. Solid symbols remain significant after the false
discovery rate, and open symbols are significant only as isolated tests.

## Usage

``` r
plot_node_test(
  x,
  tree,
  test = c("max", "site", "joint"),
  type = c("phylogram", "fan"),
  level = 0.05,
  colours = c("#0072B2", "#D55E00"),
  text_size = 2.8,
  legend = TRUE,
  ...,
  alpha = NULL,
  col = NULL,
  cex = NULL
)
```

## Arguments

- x:

  Result of
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md).

- tree:

  The `phylo` object that gave the contrasts of
  [`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md).

- test:

  Permutation test to show: `"max"` for the larger p-value of sites and
  species, `"site"` for sites alone, or `"joint"` for the permutation of
  the article. See
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  for what each one tests.

- type:

  `"phylogram"`, as in the article, or `"fan"`, for very large trees,
  with the names of the variables along the radii so that they stay
  clear of the species.

- level:

  Significance level.

- colours:

  Colours of predictors and biogeographic events.

- text_size:

  Size of the text, in mm (as in
  [`ggplot2::geom_text()`](https://ggplot2.tidyverse.org/reference/geom_text.html)).

- legend:

  Draw the legend?

- alpha:

  Deprecated, use `level`.

- col, cex, ...:

  Deprecated: `col` is now `colours` and `cex` is `text_size` (in mm,
  about 3.5 times the old value); arguments of
  [`ape::plot.phylo()`](https://rdrr.io/pkg/ape/man/plot.phylo.html) are
  no longer used.

## Value

A ggplot object, with the table of symbols in the attribute `marks`: one
row per symbol, with `node`, `set`, `variable`, `p`, `p_adjusted` and
`corrected`, which is `TRUE` for solid symbols.

## Details

The article does not say which tests form a family for the correction.
Here, the correction of Benjamini & Hochberg (1995) uses all nodes and
variables of each set, so predictors and events are corrected
separately.

As the smallest p-value is 1 / (`n_perm` + 1), solid symbols need many
permutations. With few species, `test = "max"` may leave no symbols at
all (see
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)).

## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
rate: a practical and powerful approach to multiple testing. *Journal of
the Royal Statistical Society: Series B*, 57, 289–300.
[doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
phylogenetics: separating the roles of environmental filters and
historical biogeography. *Ecology Letters*, 13, 1290–1299.
[doi:10.1111/j.1461-0248.2010.01523.x](https://doi.org/10.1111/j.1461-0248.2010.01523.x)

## See also

Other clades and biogeography:
[`biogeography_matrix()`](https://provetelab.org/comdynamite/reference/biogeography_matrix.md),
[`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md),
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md),
[`plot_node_heatmap()`](https://provetelab.org/comdynamite/reference/plot_node_heatmap.md)

## Examples

``` r
nodes <- node_contrasts(daphniids$tree)
set.seed(1)
out <- node_test(daphniids$comm, nodes, daphniids$env, daphniids$biogeography,
                 n_perm = 49)
plot_node_test(out, daphniids$tree, test = "site")

# \donttest{
# solid symbols need more permutations
out <- node_test(daphniids$comm, nodes, daphniids$env, daphniids$biogeography,
                 n_perm = 999)
plot_node_test(out, daphniids$tree, test = "site")

# }
```
