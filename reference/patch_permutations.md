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
  1.5 times the median distance from a site to its fourth nearest site,
  joins diagonal neighbours in a full grid where most sites are not on
  an edge, and the four nearest sites are one step away, and grows with
  the spacing of sparse sites, where it follows the distance between
  neighbouring sites rather than the nearest pair.

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
sites, as cells with records in a large region, the default `gap`
follows the spacing of the sites, so that neighbouring sites of an area
form a patch instead of moving one by one. Results should still be
compared across values of `gap`.

This method was built for this package and has no published evaluation.
In `data-raw/validation_edges.R`, with communities that followed a
hidden autocorrelated environment and areas without effect placed at
random, it rejected 3.0% to 6.3% of the tests at the 5% level in a
square, a strip and an irregular band, against 44% to 55% with free
permutations (see the article on analysing clades). With autocorrelation
at scales up to that of the landscape it rejected 2.2% to 6.0%, and free
permutations 49% to 75%.

These rates average over places of areas, and a real analysis has one.
Patches are placed independently of each other, so the distance between
the two sides of an event is not kept, and sides tend to end close
together. With areas fixed in a square of 20 by 20 sites, the rate was
0.4% to 2% for adjacent areas, and 6% to 10% for areas in opposite
corners or halves of the landscape. The test is therefore liberal for
sides far apart, as across a barrier, and conservative for neighbouring
sides. With 25% to 40% of the sites of that square and areas in opposite
corners (`data-raw/validation_gap.R`, 150 replicates and 99
permutations), it rejected 15% to 17% with 1.5 times the distance to the
nearest site as `gap`, the former default, and 5.6% to 5.9% with the
current default, which follows the fourth nearest site; those two cases
were used to choose it. In sparse landscapes held out from that choice
(irregular areas with 40% of the sites, adjacent areas with 60%, random
blocks with 80% and irregular areas on seven clusters of sites) it
rejected 2.1% to 5.9%, against 2.6% to 8.4% with the former default. It
is still far better than free permutations, which rejected 48% to 74%,
but p-values near the threshold need caution.
[`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md)
keeps those distances, and is the better choice for sides close enough
to have room to move.

## See also

[`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md)
keeps the shape of the sides and the distances between them, and
[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md)
gives a torus.
[`biogeography_matrix()`](https://provetelab.org/comdynamite/reference/biogeography_matrix.md)
codes the events.

Other restricted null models:
[`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md),
[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md)

## Examples

``` r
coords <- expand.grid(x = 1:10, y = 1:6)
event <- ifelse(coords$x <= 3 & coords$y <= 3, 1, ifelse(coords$x >= 8, -1, 0))
permutations <- patch_permutations(coords, cbind(barrier = event), n_perm = 3)

# where the first permutation places the two sides
matrix(event[permutations$barrier[1, ]], 10)
#>       [,1] [,2] [,3] [,4] [,5] [,6]
#>  [1,]    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0
#>  [3,]    0    0    0    0    0    0
#>  [4,]    0    0    0    0    0   -1
#>  [5,]    0    0    0   -1   -1   -1
#>  [6,]    0    0    0   -1   -1   -1
#>  [7,]    1    0   -1   -1   -1   -1
#>  [8,]    1    1    0   -1   -1   -1
#>  [9,]    1    1    1   -1   -1   -1
#> [10,]    1    1    1    0    0   -1
```
