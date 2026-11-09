# Random partitions of a tree like a given one

`random_regimes()` draws a partition of the edges of `tree` with the
same regimes as `regimes` and clades of similar size, placed at random:
the null partitions of
[`damocles_partition_test()`](https://provetelab.org/comdynamite/reference/damocles_partition_test.md).
Each regime other than that of the root is taken as one or more clades,
each starting at an edge where the regime changes (as a shift on that
branch in BAMM, or a clade with its stem in
[`clade_regimes()`](https://provetelab.org/comdynamite/reference/clade_regimes.md));
each clade is moved to an edge, drawn at random, whose subtree holds
between `1 / tolerance` and `tolerance` times as many edges, from the
largest clade to the smallest, so that nested clades stay nested when
sizes allow it. Draws in which a later clade covers an earlier one
completely are discarded.

## Usage

``` r
random_regimes(tree, regimes, tolerance = 2, max_overlap = 0.25)
```

## Arguments

- tree:

  A rooted `phylo` tree.

- regimes:

  Regime of each edge, in the order of `tree$edge`.

- tolerance:

  Largest ratio between the sizes of a clade and its random replacement.

- max_overlap:

  Largest Jaccard index between the edges of a clade and those of its
  random replacement, between 0 and 1.

## Value

A factor with one regime per edge, with the levels of `regimes`.

## Details

Edges whose subtree largely coincides with the original clade, such as
the edge just above it or those of its large daughter clades, give
partitions that are nearly the original one. `max_overlap` excludes
them: the Jaccard index between the edges of the original clade and
those of its replacement (shared edges over edges in either) must not
exceed it. With `max_overlap = 1` any edge of the right size can be
drawn.

## See also

[`damocles_partition_test()`](https://provetelab.org/comdynamite/reference/damocles_partition_test.md),
[`clade_regimes()`](https://provetelab.org/comdynamite/reference/clade_regimes.md).

## Examples

``` r
set.seed(8)
tr <- treeshapesim::sim_tree_shape(30)
reg <- clade_regimes(tr, list(A = ape::Ntip(tr) + 2))
table(reg)
#> reg
#> background          A 
#>         30         28 
table(random_regimes(tr, reg))
#> 
#> background          A 
#>         32         26 
```
