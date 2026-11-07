# Matrix of biogeographic events

Codes the events of
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
from the area of each site, for example the refugium, basin or island
where it lies. Each event gets 1 for sites on one side, -1 for sites on
the other and 0 for sites on neither.

## Usage

``` r
biogeography_matrix(
  area,
  contrast = "others",
  outside = c("zero", "other"),
  min_sites = 5
)
```

## Arguments

- area:

  Area of each site, with `NA` for sites outside all areas, or a matrix
  or data frame with one row per site and one column per area, with 1 or
  `TRUE` for the sites of each area. Polygons can give it with
  [`sf::st_join()`](https://r-spatial.github.io/sf/reference/st_join.html)
  or
  [`sf::st_intersects()`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html).

- contrast:

  Named list of events, each a list with the areas of the two sides, or
  `"others"` for one event per area.

- outside:

  For `contrast = "others"`, sites outside all areas are left out of
  each event with `"zero"` or count as the other side with `"other"`.

- min_sites:

  A warning is given for events with fewer sites on one side. This is a
  rule of thumb, not a value from the literature.

## Value

Matrix with one row per site and one column per event, with values of 1,
-1 and 0.

## Details

Events should come from a hypothesis stated in advance. Leibold et al.
(2010) contrasted areas on the two sides of historical barriers, and
Leibold & Chase (2018, chapter 9) a published biogeographic hypothesis.
A named list in `contrast` gives this design, and nested events are
written as nested groups of areas. `contrast = "others"` compares each
area with the rest, which does not need that hypothesis, but asks a
different question and is better seen as exploratory.

A site on both sides of an event gets 0.

The sites of an area are often neighbours, so their number overstates
the information in an event. See `site_permutations` in
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
and
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md).

## References

Leibold, M. A., & Chase, J. M. (2018). *Metacommunity Ecology*.
Monographs in Population Biology 59. Princeton University Press.
[doi:10.2307/j.ctt1wf4d24](https://doi.org/10.2307/j.ctt1wf4d24)

Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
phylogenetics: separating the roles of environmental filters and
historical biogeography. *Ecology Letters*, 13, 1290–1299.
[doi:10.1111/j.1461-0248.2010.01523.x](https://doi.org/10.1111/j.1461-0248.2010.01523.x)

## Examples

``` r
area <- rep(c("north", "centre", "south", NA), c(8, 6, 6, 10))

# a barrier between the north and the rest, and a later one between centre and south
events <- list(barrier = list("north", c("centre", "south")),
               river = list("centre", "south"))
colSums(biogeography_matrix(area, events) == 1)
#> barrier   river 
#>       8       6 

# each area against the others
head(biogeography_matrix(area, "others"))
#>      centre north south
#> [1,]     -1     1    -1
#> [2,]     -1     1    -1
#> [3,]     -1     1    -1
#> [4,]     -1     1    -1
#> [5,]     -1     1    -1
#> [6,]     -1     1    -1
```
