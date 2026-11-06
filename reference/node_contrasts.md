# Node contrasts of a phylogeny

Codes each node of a tree as a contrast between its two descendant
clades, as in Leibold et al. (2010). One clade starts at -1 and the
other at 1, and the value is halved at each node on the way to the tips.
Species outside the node get zero. Nodes that join two tips are left
out, as in the original code, because they contrast two single species.

## Usage

``` r
node_contrasts(tree)
```

## Arguments

- tree:

  A `phylo` object with dichotomous nodes.

## Value

Matrix with one row per tip and one column per node, named by the node
numbers of `tree`.

## References

Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
phylogenetics: separating the roles of environmental filters and
historical biogeography. *Ecology Letters*, 13, 1290–1299.
[doi:10.1111/j.1461-0248.2010.01523.x](https://doi.org/10.1111/j.1461-0248.2010.01523.x)

## Examples

``` r
node_contrasts(ape::read.tree(text = "((a,b),(c,(d,e)));"))
#>       6    8
#> a -0.50  0.0
#> b -0.50  0.0
#> c  0.50 -1.0
#> d  0.25  0.5
#> e  0.25  0.5
```
