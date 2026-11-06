# Site permutations that keep the structure of a torus

The translations of a `side` x `side` torus combined with the eight
symmetries of the square keep all distances among sites, and therefore
the spatial autocorrelation of the environment. They give a restricted
null model for spatially explicit metacommunities (see the vignette).

## Usage

``` r
torus_permutations(side)
```

## Arguments

- side:

  Number of cells per side, as in
  [`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md).

## Value

Integer matrix with one permutation of the site indices per row,
excluding the identity (`8 * side^2 - 1` rows when `side > 2`).

## Examples

``` r
p <- torus_permutations(4)
dim(p)
#> [1] 127  16
```
