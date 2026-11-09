# Plot the semipartial correlations of a node test

Draws a table of nodes by variables with the semipartial correlations of
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md),
as in Table S1 of Leibold et al. (2010). Bold values with an asterisk
remain significant after the false discovery rate, as in
[`plot_node_test()`](https://provetelab.org/comdynamite/reference/plot_node_test.md);
other bold values are significant only as isolated tests. Blank cells
could not be calculated.

## Usage

``` r
plot_node_heatmap(
  x,
  test = c("max", "site", "joint"),
  level = 0.05,
  colours = c("#0072B2", "#D55E00"),
  text_size = 2.5,
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

- test:

  Permutation test to show: `"max"` for the larger p-value of sites and
  species, `"site"` for sites alone, or `"joint"` for the permutation of
  the article. See
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
  for what each one tests.

- level:

  Significance level.

- colours:

  Colours of negative and positive values.

- text_size:

  Size of the text, in mm (as in
  [`ggplot2::geom_text()`](https://ggplot2.tidyverse.org/reference/geom_text.html)).

- alpha:

  Deprecated, use `level`.

- col, cex, ...:

  Deprecated: `col` is now `colours` and `cex` is `text_size`; arguments
  of [`graphics::image()`](https://rdrr.io/r/graphics/image.html) are no
  longer used.

## Value

A ggplot object, with the matrix of semipartial correlations (predictors
followed by biogeographic events) in the attribute `values`.

## Details

Rows are the node numbers of the tree, as in `tree$edge`;
[`ape::nodelabels()`](https://rdrr.io/pkg/ape/man/nodelabels.html) on a
plot of
[`ape::plot.phylo()`](https://rdrr.io/pkg/ape/man/plot.phylo.html) shows
where they are.

## References

Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
phylogenetics: separating the roles of environmental filters and
historical biogeography. *Ecology Letters*, 13, 1290–1299.
[doi:10.1111/j.1461-0248.2010.01523.x](https://doi.org/10.1111/j.1461-0248.2010.01523.x)

## See also

Other clades and biogeography:
[`biogeography_matrix()`](https://provetelab.org/comdynamite/reference/biogeography_matrix.md),
[`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md),
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md),
[`plot_node_test()`](https://provetelab.org/comdynamite/reference/plot_node_test.md)

## Examples

``` r
nodes <- node_contrasts(daphniids$tree)
set.seed(1)
out <- node_test(daphniids$comm, nodes, daphniids$env, daphniids$biogeography,
                 n_perm = 49)
plot_node_heatmap(out, test = "site")
```
