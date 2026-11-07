# Metacommunity data in one object

Keeps the community and what describes its species and sites together,
with names checked once. The analyses accept an `mc_data` object as
their first argument and take from it what the call leaves out, so
`pcps_test(d)` uses `d$tree` and `d$env`. Simulations convert with
`as_mc_data()`.

## Usage

``` r
mc_data(
  comm,
  tree = NULL,
  traits = NULL,
  env = NULL,
  coords = NULL,
  areas = NULL,
  series = NULL
)

as_mc_data(x, ...)

# Default S3 method
as_mc_data(x, ...)

# S3 method for class 'mc_data'
as_mc_data(x, ...)

# S3 method for class 'mc_sim'
as_mc_data(x, ...)

# S3 method for class 'mc_data'
print(x, ...)

# S3 method for class 'mc_data'
x[i, j, drop = FALSE]

# S3 method for class 'mc_data'
summary(object, ...)

# S3 method for class 'mc_data'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- comm:

  Matrix or data frame with one row per site and one column per species,
  named by species. Sites without names are called `site1`, `site2` and
  so on.

- tree:

  Optional `phylo` object.

- traits:

  Optional data frame, matrix or named vector with one row per species,
  named by species.

- env:

  Optional data frame, matrix or vector with one row per site.

- coords:

  Optional matrix or data frame of coordinates with one row per site.

- areas:

  Optional vector with the area or biogeographic event of each site, for
  [`biogeography_matrix()`](https://provetelab.org/comdynamite/reference/biogeography_matrix.md).

- series:

  Optional array of sites x species x steps, such as the `history` of a
  simulation.

- x, object:

  An `mc_data` object, or for `as_mc_data()` an object to convert.

- ...:

  Passed to `mc_data()` by `as_mc_data()`; not used elsewhere.

- i, j:

  Sites and species to keep, by position, name or logical vector.

- drop, row.names, optional:

  Not used.

## Value

`mc_data()` and `as_mc_data()` give a list of class `mc_data` with
`comm`, `tree`, `traits` (data frame), `env` (data frame), `coords`
(matrix), `areas` and `series`, each `NULL` when absent. `[` gives the
same for a subset, [`summary()`](https://rdrr.io/r/base/summary.html)
gives
[`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md)
and [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) a
long table with `site`, `species`, `step` (`NA` without `series`) and
`abundance`, without absences.

## Details

Species are matched by name between `comm`, `tree` and `traits`, and
those missing from any of them are dropped with a warning, as in the
analyses. Species without individuals stay, since they belong to the
regional pool. Sites are matched by name when `env`, `coords`, `areas`
or `series` have names, and by position otherwise.

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(12)
traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100), rmax_range = c(0, 30))
env <- seq(0, 100, length.out = 20)
comm <- simulate_static(traits, env, incidence = TRUE)$abundance
d <- mc_data(comm, tree, traits["optimum"], env = data.frame(env = env))
d
#> <mc_data> 20 sites x 12 species
#> tree: 12 tips
#> traits: optimum
#> env: env
d[1:10, ]
#> <mc_data> 10 sites x 12 species
#> tree: 12 tips
#> traits: optimum
#> env: env
pcps_test(d[rowSums(d$comm) > 0, ], n_perm = 9)
#> <mc_pcps> principal coordinates of phylogenetic structure
#> permutations: 9 free of sites, 9 of species; smallest possible p: 0.1
#>   test share    r2 statistic p_site p_species p_max
#>  whole 1.000 0.235      5.53    0.1       0.5   0.5
#>  pcps1 0.277 0.735     49.95    0.1       0.2   0.2
#>  pcps2 0.237 0.093      1.84    0.1       0.6   0.6
```
