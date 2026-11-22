# Plot the phylogenetic eigenvectors related to a predictor

Draws the fourth-corner correlation of each phylogenetic eigenvector
with the predictor, in the order of the eigenvalues, as computed by
[`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md).
The eigenvectors summed into `pev` are drawn in `colour` and the others
in grey; the subtitle gives the largest p-value of the global test and
the correlations of `pev` and `pev_brown`.

## Usage

``` r
# S3 method for class 'mc_phylo_vector'
plot(x, colour = "#0072B2", ...)
```

## Arguments

- x:

  Result of
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md).

- colour:

  Colour of the retained eigenvectors.

- ...:

  Not used.

## Value

A ggplot object.

## See also

[`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md)

Other environment and space:
[`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md),
[`plot_trait_partition()`](https://provetelab.org/comdynamite/reference/plot_trait_partition.md),
[`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md)

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(20)
traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100),
                         rmax_range = c(0, 30))
env <- runif(50, 0, 100)
comm <- simulate_static(traits, env, incidence = TRUE)$abundance
keep <- rowSums(comm) > 0
plot(phylo_vector(comm[keep, ], tree, env[keep], n_perm = 99, n_sim = 100))
```
