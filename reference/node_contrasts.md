# Node contrasts of a phylogeny

Codes each node of a tree as a contrast between its two descendant
clades, as in Leibold et al. (2010). One clade starts at -1 and the
other at 1, and the value is halved at each node on the way to the tips.
Species outside the node get zero.

## Usage

``` r
node_contrasts(tree, cherries = TRUE, branch_lengths = FALSE)
```

## Arguments

- tree:

  A `phylo` object with dichotomous nodes.

- cherries:

  Keep nodes that join two tips? The article tests them (its Figure 3),
  but the function `CreateNodeContrasts` of the authors leaves them out.

- branch_lengths:

  Use the branch lengths of `tree`?

## Value

Matrix with one row per tip and one column per node, named by the node
numbers of `tree`.

## Details

With `branch_lengths = TRUE`, values follow section III of the
supplement of Leibold et al. (2010): each branch between a species and
the node counts its length divided by the number of species below it,
and the sum is divided by the total length of that side of the node.
Species that share history get less weight, and each side still sums
to 1. The supplement proposes this coding for trees in which all species
are equally distant from the root.

## Note

Adapted from the MATLAB function `CreateNodeContrasts` in the supplement
of Leibold et al. (2010).

## References

Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
phylogenetics: separating the roles of environmental filters and
historical biogeography. *Ecology Letters*, 13, 1290–1299.
[doi:10.1111/j.1461-0248.2010.01523.x](https://doi.org/10.1111/j.1461-0248.2010.01523.x)

## See also

Other clades and biogeography:
[`biogeography_matrix()`](https://provetelab.org/comdynamite/reference/biogeography_matrix.md),
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md),
[`plot_node_heatmap()`](https://provetelab.org/comdynamite/reference/plot_node_heatmap.md),
[`plot_node_test()`](https://provetelab.org/comdynamite/reference/plot_node_test.md)

## Examples

``` r
node_contrasts(ape::read.tree(text = "((a,b),(c,(d,e)));"))
#>       6  7    8  9
#> a -0.50 -1  0.0  0
#> b -0.50  1  0.0  0
#> c  0.50  0 -1.0  0
#> d  0.25  0  0.5 -1
#> e  0.25  0  0.5  1

# Figure S2 of the supplement
tree <- ape::read.tree(text = "((((A:0.2,B:0.2):0.65,C:0.85):0.15,D:1):0.4,(E:0.5,F:0.5):0.9);")
round(node_contrasts(tree, branch_lengths = TRUE), 3)
#>        7      8    9 10 11
#> A -0.196 -0.280 -0.5 -1  0
#> B -0.196 -0.280 -0.5  1  0
#> C -0.290 -0.439  1.0  0  0
#> D -0.319  1.000  0.0  0  0
#> E  0.500  0.000  0.0  0 -1
#> F  0.500  0.000  0.0  0  1
```
