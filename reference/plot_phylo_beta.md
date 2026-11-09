# Community matrix ordered by phylogenetic beta diversity

The heatmap that Peres-Neto & Kembel (2015) suggest for
[`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md):
sites follow an axis of the matrix predicted by the environment and
species follow the tree, so that a link between environment and
phylogeny appears as clades that replace each other along the sites.

## Usage

``` r
plot_phylo_beta(
  x,
  comm,
  tree,
  axis = 1,
  colours = "#0072B2",
  text_size = 2.5,
  ...,
  col = NULL,
  cex = NULL
)
```

## Arguments

- x:

  Result of
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md).

- comm:

  The community matrix given to
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md).

- tree:

  The `phylo` object given to
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md).

- axis:

  Axis of `x$scores` that orders the sites.

- colours:

  Colour of the largest value, with white for absences.

- text_size:

  Size of the names of species, in mm.

- col, cex, ...:

  Deprecated: `col` is now `colours` and `cex` is `text_size`; arguments
  of [`graphics::image()`](https://rdrr.io/r/graphics/image.html) are no
  longer used.

## Value

A ggplot object, with the ordered matrix in the attribute `ordered`.

## References

Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient
analysis: environmental drivers of phylogenetic variation across
ecological communities. *Plant Ecology*, 216, 709–724.
[doi:10.1007/s11258-014-0405-0](https://doi.org/10.1007/s11258-014-0405-0)

## See also

Other the whole tree:
[`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md),
[`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md),
[`phylo_beta()`](https://provetelab.org/comdynamite/reference/phylo_beta.md),
[`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md),
[`phylo_eigenvectors()`](https://provetelab.org/comdynamite/reference/phylo_eigenvectors.md),
[`plot.mc_pcps()`](https://provetelab.org/comdynamite/reference/plot.mc_pcps.md)

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(20)
traits <- species_traits(tree, lambda = 1, optima_range = c(-1, 101),
                         rmax_range = c(0, 30))
env <- runif(50, 0, 100)
comm <- simulate_static(traits, env, incidence = TRUE)$abundance
keep <- rowSums(comm) > 0
out <- phylo_beta_test(comm[keep, ], tree, env[keep], n_perm = 9)
plot_phylo_beta(out, comm[keep, ], tree)
```
