# Permutations that move the patches of biogeographic events

A restricted null model for the events of
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
in landscapes with edges, such as refugia in a biome or islands in an
archipelago. Each area is split into patches of neighbouring sites. In
each permutation, every patch is placed around a site drawn at random,
on the sites nearest to it that are still free, from the smallest patch
to the largest. The number and size of the patches are kept, and so is
the fact that sites of a patch are neighbours, while their place
changes. Permutations of a torus cannot do this, because a landscape
with edges does not wrap around.

## Usage

``` r
patch_permutations(coords, biogeography, n_perm = 999, gap = NULL)
```

## Arguments

- coords:

  Matrix or data frame with the two coordinates of each site, in the
  order of the rows of `biogeography`.

- biogeography:

  Matrix of biogeographic events with one row per site and values of 1,
  -1 and 0, as given to
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md).

- n_perm:

  Number of permutations.

- gap:

  Largest distance between neighbouring sites of a patch. The default,
  1.5 times the median distance to the nearest site, joins diagonal
  neighbours in a regular grid.

## Value

List with one integer matrix per event, each with one permutation per
row and the names of the sites as column names, for `area_permutations`
in
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md).
Rows only move the areas of that event, so they should not be used as
`site_permutations`.

## Details

Patches take a compact form near their random site, so their shape is
not kept, and a patch next to an edge or to another patch may take an
elongated form or be split. Distances are Euclidean, so coordinates
should be projected. A landscape with holes or with parts far apart,
such as islands, may give patches that cross the gaps. With sparse
sites, as cells with records in a large region, the default `gap` splits
areas into many patches of one site, which move almost freely, so `gap`
should follow the scale of the area. Results should be compared across
values of `gap`.

This method was built for this package and has no published evaluation.
In `data-raw/validation_edges.R`, with communities that followed a
hidden autocorrelated environment and areas without effect, it rejected
3.0% to 6.3% of the tests at the 5% level in a square, a strip and an
irregular band, against 44% to 55% with free permutations (see the
article on analysing clades). With autocorrelation at scales up to that
of the landscape it rejected 2.2% to 6.0%, and free permutations 49% to
75%. Areas much larger than those simulated were not evaluated.

## See also

[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md)
for a torus.
[`biogeography_matrix()`](https://provetelab.org/comdynamite/reference/biogeography_matrix.md)
codes the events.

## Examples

``` r
coords <- expand.grid(x = 1:10, y = 1:6)
event <- ifelse(coords$x <= 3 & coords$y <= 3, 1, ifelse(coords$x >= 8, -1, 0))
permutations <- patch_permutations(coords, cbind(barrier = event), n_perm = 3)

# where the first permutation places the two sides
matrix(event[permutations$barrier[1, ]], 10)
#>       [,1] [,2] [,3] [,4] [,5] [,6]
#>  [1,]    0    0    0    0    0    0
#>  [2,]    0    0    0   -1   -1    0
#>  [3,]    0   -1   -1   -1   -1   -1
#>  [4,]    0   -1   -1   -1   -1   -1
#>  [5,]    0   -1   -1   -1   -1    0
#>  [6,]    0    0   -1   -1    0    0
#>  [7,]    0    0    1    1    1    0
#>  [8,]    0    0    1    1    1    0
#>  [9,]    0    0    1    1    1    0
#> [10,]    0    0    0    0    0    0
```
