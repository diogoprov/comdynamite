# Site permutations that keep the structure of a torus

The translations of a torus combined with the symmetries of the grid
keep all distances among sites, and therefore the spatial
autocorrelation of the environment. They give a restricted null model
for spatially explicit metacommunities (see the vignette), as the
torus-translation test of Harms et al. (2001). A square grid has eight
symmetries and a rectangular one four. Sites are numbered with the first
coordinate running fastest, as in `expand.grid(x, y)` and
[`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md).
A real landscape has edges, and this null treats opposite edges as
adjacent: see
[`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md)
and
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md).

## Usage

``` r
torus_permutations(side, symmetries = TRUE, n_perm = NULL)
```

## Arguments

- side:

  Number of cells per side, as in
  [`landscape_grid()`](https://provetelab.org/comdynamite/reference/landscape_grid.md),
  or two numbers for a rectangular grid, the first for the coordinate
  that runs fastest.

- symmetries:

  Combine the translations with the reflections and rotations of the
  grid?

- n_perm:

  Number of permutations to keep, drawn without replacement. The default
  keeps all of them.

## Value

Integer matrix with one permutation of the site indices per row,
excluding the identity (\\8\\ \mathrm{side}^2 - 1\\ rows on a square
grid with `side > 2`).

## References

Harms, K. E., Condit, R., Hubbell, S. P., & Foster, R. B. (2001).
Habitat associations of trees and shrubs in a 50-ha neotropical forest
plot. *Journal of Ecology*, 89, 947–959.

## Examples

``` r
p <- torus_permutations(4)
dim(p)
#> [1] 127  16
dim(torus_permutations(c(6, 4), symmetries = FALSE))
#> [1] 23 24
```
