# Permutations that move biogeographic events as rigid bodies

A restricted null model for the events of
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
that keeps the shape of each side and the distances between the sides.
Each permutation moves all the sites of an event together, by a rotation
or a reflection and a translation. A move is kept when its sites land
inside the landscape: within `gap` of a site, or in a hole between
sites, but not beyond an edge. Sites that land where there is no site,
as in a hole of a sparse landscape, go to the nearest free site within
twice `gap`, and a fraction `outside` of the sites of an event may land
beyond the landscape and be brought in the same way. Moves are never
carried across the landscape, as in a torus, which would split a side.
This is the torus-translation test of Harms et al. (2001) without the
wrap-around.

## Usage

``` r
rigid_permutations(
  coords,
  biogeography,
  n_perm = 999,
  gap = NULL,
  outside = 0.1,
  independent = Inf,
  max_tries = 20 * n_perm
)
```

## Arguments

- coords:

  Matrix or data frame with the two coordinates of each site, in the
  order of the rows of `biogeography`. Distances are Euclidean, so they
  should be projected.

- biogeography:

  Matrix of biogeographic events with one row per site and values of 1,
  -1 and 0, as given to
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md).

- n_perm:

  Largest number of permutations. Fewer are returned when fewer distinct
  placements of the event exist.

- gap:

  A moved site counts as inside the landscape when it is closer than
  `gap` to a site, or in a hole narrower than twice `gap` between sites
  (four times `gap` off a lattice, where sites leave larger voids
  between neighbours), but not beyond an edge. A site that lands inside
  but not on a site goes to the nearest free site within twice `gap`.
  The default is the median distance from a site to its fourth nearest
  site: the step of a full grid, where it makes a site just beyond the
  edge count as outside, and the spacing of sparse sites, so that holes
  count as inside. It is not the `gap` of
  [`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md),
  which is 1.5 times this one by default.

- outside:

  Largest fraction of the sites of an event that may land outside the
  landscape in a move. Those sites go to the nearest free site within
  twice `gap`, and the move is dropped when one has none.

- independent:

  Distance beyond which parts of an event move independently of each
  other, while staying at least that far apart. The default, `Inf`,
  moves the whole event together and assumes nothing about the range of
  autocorrelation.

- max_tries:

  Largest number of moves tried: of the enumerated moves of one part, in
  random order, when placements fail or repeat a layout, or of the
  combinations of several independent parts, which cannot all be
  enumerated.

## Value

List with one integer matrix per event, each with one permutation per
row and the names of the sites as column names, for `area_permutations`
in
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md).
Rows only move the sites of that event, so they should not be used as
`site_permutations`.

## Details

The moves are enumerated, so the permutations are distinct placements of
the event, up to the sign of the contrast, which does not change the
test, and drawn at random when there are more than `n_perm`. Sites on a
regular lattice get all the moves of the lattice: its symmetries with
every translation in cells. Other coordinates get rotations in steps
that move the site farthest from the centre of the event by about `gap`,
at least 4 and at most 360 of them, each with and without a reflection,
and the translations that put the site nearest the centre on each site
of the landscape; the other sites go to the nearest free sites. Which
moves fit is then fixed by the landscape, and only the draw among them
changes between runs. In an event without zeros where one side has up to
a third of the sites, only that side moves and the other is the rest of
the landscape, as when one area is compared with all the others.

An event that spans the landscape, as two halves, two strips along
opposite edges or areas in opposite corners, has few placements,
sometimes none, because a move would push a side out. The function then
returns fewer rows than `n_perm`, with a warning, or stops. The smallest
p-value in
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
is 1 / (rows + 1), so with fewer than 19 rows no node can reach the 5%
level. This is a limit of the event, not of the function: a permutation
test that keeps the shape of two halves has almost nothing to permute.

With `independent`, parts of the event farther apart than that distance
move independently of each other, but keep at least that distance
between them. This gives placements to sides far apart, as across a
barrier, at the price of assuming that the community is not spatially
autocorrelated beyond `independent`. When that assumption fails the test
is liberal, and no value of `independent` can be checked from the data
alone, so results should be compared across values. A large part in an
irregular region may have few moves that fit. It then stays nearly in
place while the small parts move, which gives many placements that say
little, and the function warns.

This method was built for this package and has no published evaluation.
In `data-raw/validation_rigid.R`, with communities that followed a
hidden autocorrelated environment and areas without effect, 150
replicates and 199 permutations, the default rejected 3.7% to 4.7% of
the tests at the 5% level with areas placed at random in a square of 20
by 20 sites and in a strip of 50 by 8. With areas fixed in the square,
it rejected 2.2% to 3.5% for adjacent areas and for several small areas,
and 4.2% to 4.8% for adjacent areas, an L shape and two areas at a short
distance when the hidden environment had an exponential covariance.
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
rejected 1.4% to 6.4% in those cases, and is conservative for adjacent
areas. Where the event had room to move, power was equal to or higher
than that of
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md),
up to twice for adjacent areas (24% against 13%); with few placements it
was lower (4% against 13% for blocks in the band).

With sparse sites (`data-raw/validation_gap.R`, 150 replicates and 99
permutations), the default rejected 5.5% to 6.1% where the event had
room to move: irregular areas with 40% of the sites of the square,
adjacent areas with 60% and random blocks with 80%. Areas in opposite
corners with 25% to 40% of the sites, and irregular areas on seven
clusters of sites, had no placement or fewer than five, so the default
does not test them.

Off a lattice (`data-raw/validation_offgrid.R`, 150 replicates and 99
permutations), the default rejected 4.7% for adjacent areas in the
square with jitter and 5.9% in 400 sites placed at random, and 6.8% to
8.0% for two blocks placed at random in those landscapes, in the band
with jitter and in a crescent of random sites, where it is liberal.
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
rejected 1.2% to 6.3% in those cases. With a second generator in an
L-shaped region of random sites (`data-raw/validation_nulls_check.R`,
150 replicates and 199 permutations), it rejected 7.0% for two fixed
areas and 8.7% for areas placed at random, against 5.5% and 6.6% with
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md).
Off a lattice this function is thus liberal, and
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
is the safer choice unless its conservatism for adjacent areas matters.
Power for adjacent areas with jitter was 23% against 12%. With 400
sites, the enumeration and 99 placements took 0.3 to 2.8 s; with 750
sites placed at random and an event of 57 sites, 999 placements took 9
s, and with 750 sites on a lattice with jitter and an event of 225
sites, 701 placements took 12 min, since every site of a large event has
to be moved to a free site in each try.

Sides far apart remain a problem. In an irregular band most random areas
had no placement or fewer than 20, so the default does not test them.
With `independent = 10`, random areas in the three landscapes gave 3.9%
to 5.5%, but areas in opposite corners of the square gave 6% to 7% with
autocorrelation over 3 or 4 sites and 12% over 12 sites, strips along
opposite edges 8% to 9% (2% to 3% with
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)),
and areas in opposite corners with sparse sites 10% to 12%.

## References

Harms, K. E., Condit, R., Hubbell, S. P., & Foster, R. B. (2001).
Habitat associations of trees and shrubs in a 50-ha neotropical forest
plot. *Journal of Ecology*, 89, 947–959.

## See also

[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
moves each patch on its own, and
[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md)
gives a torus.
[`biogeography_matrix()`](https://provetelab.org/comdynamite/reference/biogeography_matrix.md)
codes the events.

Other restricted null models:
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md),
[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md)

## Examples

``` r
coords <- expand.grid(x = 1:10, y = 1:6)
event <- ifelse(coords$x <= 3 & coords$y <= 3, 1, ifelse(coords$x >= 5 & coords$x <= 6, -1, 0))
permutations <- rigid_permutations(coords, cbind(barrier = event), n_perm = 3)

# where the first permutation places the two sides
matrix(event[permutations$barrier[1, ]], 10)
#>       [,1] [,2] [,3] [,4] [,5] [,6]
#>  [1,]    0    0    0    0    0    0
#>  [2,]    0    0    0    0    0    0
#>  [3,]    0    0    0    0    0    0
#>  [4,]    0    0    0    0    0    0
#>  [5,]   -1   -1    0    0    0    0
#>  [6,]   -1   -1    0    0    0    0
#>  [7,]   -1   -1    0    0    0    0
#>  [8,]   -1   -1    0    1    1    1
#>  [9,]   -1   -1    0    1    1    1
#> [10,]   -1   -1    0    1    1    1
```
