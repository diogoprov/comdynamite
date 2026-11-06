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

## Examples

``` r
tr <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
paint_regimes(tr, nodes = 6)
#> [1] "0" "1" "1" "0" "0" "0"
#> attr(,"order")
#> [1] "cladewise"
```
