# Plot a test of phylogenetic community composition

Draws the adjusted R2 of each phylogenetic eigenvector regressed on the
predictors, as computed by
[`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md),
with the mean of the positive values (`r2_mean`) as a dashed line and
the largest p-value in the subtitle. The test statistic is the variance
of these bars.

## Usage

``` r
# S3 method for class 'mc_pcc'
plot(x, colour = "#0072B2", ...)
```

## Arguments

- x:

  Result of
  [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md).

- colour:

  Colour of the bars.

- ...:

  Not used.

## Value

A ggplot object.

## See also

[`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md)

Other the whole tree:
[`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md),
[`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md),
[`phylo_beta()`](https://provetelab.org/comdynamite/reference/phylo_beta.md),
[`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md),
[`phylo_eigenvectors()`](https://provetelab.org/comdynamite/reference/phylo_eigenvectors.md),
[`plot.mc_pcps()`](https://provetelab.org/comdynamite/reference/plot.mc_pcps.md),
[`plot_phylo_beta()`](https://provetelab.org/comdynamite/reference/plot_phylo_beta.md)

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(20)
traits <- species_traits(tree, lambda = 1, optima_range = c(-1, 101),
                         rmax_range = c(0, 30))
env <- runif(50, 0, 100)
comm <- simulate_static(traits, env, incidence = TRUE)$abundance
keep <- rowSums(comm) > 0
plot(pcc_test(comm[keep, ], tree, env[keep], n_perm = 99))
```
