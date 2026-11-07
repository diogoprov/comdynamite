# Preparing your data

The analyses of the package take the same objects: a community matrix, a
tree, a table of predictors and, for some of them, coordinates and
areas. This article goes from tables like those of a field survey to a
first test. The data below are simulated, to keep the article
self-contained.

``` r

library(comdynamite)

set.seed(1)
tree <- treeshapesim::sim_tree_shape(10)
tree$tip.label <- paste0("Genus_species", 1:10)
sites <- expand.grid(x = 1:8, y = 1:5)
sites$site <- paste0("s", seq_len(nrow(sites)))
sites$depth <- sites$x * 10 + rnorm(nrow(sites), sd = 5)
sites$habitat <- rep(c("forest", "open"), length.out = nrow(sites))
traits <- species_traits(tree, lambda = 1, optima_range = c(0, 90), rmax_range = c(0, 30))
wide <- simulate_static(traits, sites$depth, incidence = TRUE)$abundance
records <- data.frame(site = sites$site[row(wide)[wide > 0]],
                      species = sub("_", " ", tree$tip.label[col(wide)[wide > 0]]))
head(records)
#>   site        species
#> 1   s3 Genus species1
#> 2   s4 Genus species1
#> 3   s5 Genus species1
#> 4   s6 Genus species1
#> 5   s7 Genus species1
#> 6   s8 Genus species1
```

## The community matrix

`comm` has one row per site and one column per species, with numbers
only, so names of sites go to the row names. A table of records becomes
a matrix with [`table()`](https://rdrr.io/r/base/table.html).

``` r

comm <- records |>
  with(table(site, species)) |>
  unclass()
comm[1:3, 1:3]
#>      species
#> site  Genus species1 Genus species10 Genus species2
#>   s1               0               0              1
#>   s10              0               0              1
#>   s11              1               0              1
```

A site without records does not appear in this matrix.
[`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md),
[`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md),
[`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
and
[`phylo_beta()`](https://provetelab.org/comdynamite/reference/phylo_beta.md)
need at least one individual in every site, so sites without any species
stay out of those analyses.

## The tree

Species are matched to the tips by name, in any order. The most frequent
problem is the spelling of names, as here, where the table has spaces
and the tree has underscores.

``` r

setdiff(colnames(comm), tree$tip.label)[1:3]
#> [1] "Genus species1"  "Genus species10" "Genus species2"
colnames(comm) <- gsub(" ", "_", colnames(comm))
setdiff(colnames(comm), tree$tip.label)
#> character(0)
```

Species in only one of the two are dropped with a warning that names
them, so read these warnings before reading results. Analyses that use
distances among species need branch lengths, whereas
[`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md)
also works with a topology.

## Predictors

`env` has one row per site and numeric columns. Factors become dummy
variables with
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html), without
the intercept column. Missing values are not accepted, so sites without
data need to leave `comm` and `env` together.

``` r

rownames(sites) <- sites$site
sites <- sites[rownames(comm), ]
env <- cbind(depth = sites$depth, open = model.matrix(~ habitat, sites)[, "habitatopen"])
rownames(env) <- rownames(comm)
head(env, 3)
#>        depth open
#> s1  14.10611    0
#> s10 12.64624    1
#> s11 27.60925    0
```

When `comm` and the other tables all have row names, sites are matched
by name. Otherwise they are matched by position, and the rows need to be
in the same order.

## Areas and coordinates

For historical hypotheses,
[`biogeography_matrix()`](https://provetelab.org/comdynamite/reference/biogeography_matrix.md)
builds the events of
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
from the area of each site, with `NA` for sites outside all areas.
Polygons give the area with
[`sf::st_join()`](https://r-spatial.github.io/sf/reference/st_join.html).
Events should come from a hypothesis stated in advance.

``` r

area <- ifelse(sites$x <= 3, "west", ifelse(sites$x >= 6, "east", NA))
names(area) <- rownames(comm)
events <- biogeography_matrix(area, list(barrier = list("west", "east")))
table(events)
#> events
#> -1  0  1 
#> 15 10 15
```

Sites of an area are usually neighbours, and free permutations then
reject too often.
[`patch_permutations()`](https://provetelab.org/comdynamite/reference/patch_permutations.md)
gives a null model for landscapes with edges, and
[`rigid_permutations()`](https://provetelab.org/comdynamite/reference/rigid_permutations.md)
another that keeps the distance between the sides of an event, for sides
with room to move (see the article on analysing clades). Coordinates
should be projected, because distances are Euclidean;
[`sf::st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html)
does this.

``` r

coords <- as.matrix(sites[, c("x", "y")])
area_moves <- patch_permutations(coords, events, n_perm = 99)
```

## A first test

``` r

nodes <- node_contrasts(tree)
out <- node_test(comm, nodes, env, events, n_perm = 99, area_permutations = area_moves)
out
#> <mc_node_test> 9 nodes, 9 tested
#> permutations: 99 free of sites, 99 of species, of areas barrier 99; smallest possible p: 0.01
#> nodes with p (max) <= 0.05: 0 before correction, 0 after
```

The number of permutations here keeps the article fast. Use more than
this in an analysis, because the smallest possible p-value depends on it
and it is printed with the results.

## What each analysis needs

| Function | `comm` | Tree | Predictors | Other |
|----|----|----|----|----|
| [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md) | sites by species | contrasts from [`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md) | `env` | `biogeography`, optional |
| [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md) | no empty sites | with branch lengths | `env` |  |
| [`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md) | no empty sites | with branch lengths | `env` |  |
| [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md) | no empty sites | with branch lengths | `env` |  |
| [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md) | sites by species | with branch lengths | one predictor |  |
| [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md) | sites by species | traits or [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md) scores | `env` | `space` |
