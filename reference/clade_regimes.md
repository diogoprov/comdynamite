# Regimes of the edges of a tree from a list of clades

`clade_regimes()` assigns each edge of `tree` to a regime named after
the clade it lies in, for the `regimes` of
[`fit_damocles()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md)
and
[`damocles_test()`](https://provetelab.org/comdynamite/reference/damocles_test.md):
a partition of the tree defined a priori, such as families or
biogeographic groups, against which partitions from methods that
estimate shifts in rates of diversification can be compared.

## Usage

``` r
clade_regimes(tree, clades, background = "background", stem = FALSE)
```

## Arguments

- tree:

  A rooted `phylo` tree.

- clades:

  A named list of character vectors, the tips of each clade.

- background:

  Name of the regime of the edges outside every clade.

- stem:

  If `TRUE`, the stem edge of each clade belongs to it.

## Value

A factor with one regime per edge of `tree`, in the order of
`tree$edge`, with the regimes that occur as levels: `background` first,
then the clades.

## Details

A clade is the crown group of the most recent common ancestor of its
tips: the edges below that node belong to the clade, and its stem edge
too when `stem = TRUE`; a clade of one tip is its terminal edge. If the
tips given are not monophyletic the clade covers every descendant of
their common ancestor, with a warning. Clades later in the list override
earlier ones, so nested clades go after the clades that hold them. Edges
in no clade take the regime `background`.

## See also

[`fit_damocles()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md),
[`damocles_lrt()`](https://provetelab.org/comdynamite/reference/damocles_lrt.md).

Other dynamic null model:
[`confint.mc_damocles_fit()`](https://provetelab.org/comdynamite/reference/confint.mc_damocles_fit.md),
[`damocles_loglik()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md),
[`damocles_lrt()`](https://provetelab.org/comdynamite/reference/damocles_lrt.md),
[`damocles_test()`](https://provetelab.org/comdynamite/reference/damocles_test.md),
[`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md)

## Examples

``` r
set.seed(4)
tr <- treeshapesim::sim_tree_shape(20)
# the clade of the first node after the root
reg <- clade_regimes(tr, list(A = ape::extract.clade(tr, ape::Ntip(tr) + 2)$tip.label))
table(reg)
#> reg
#> background          A 
#>          2         36 
```
