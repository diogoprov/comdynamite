# Community matrix ordered by phylogenetic beta diversity

The heatmap that Peres-Neto & Kembel (2015) suggest for
[`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md):
sites follow an axis of the matrix predicted by the environment and
species follow the tree, so that a link between environment and
phylogeny appears as clades that replace each other along the sites.

## Usage

``` r
plot_phylo_beta(x, comm, tree, axis = 1, col = "#0072B2", cex = 0.7, ...)
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

- col:

  Colour of the largest value, with white for absences.

- cex:

  Size of the text.

- ...:

  Passed to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html).

## Value

The ordered matrix, invisibly.

## References

Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient
analysis: environmental drivers of phylogenetic variation across
ecological communities. *Plant Ecology*, 216, 709–724.
[doi:10.1007/s11258-014-0405-0](https://doi.org/10.1007/s11258-014-0405-0)

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
