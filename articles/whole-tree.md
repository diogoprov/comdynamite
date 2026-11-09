# Testing the whole tree

[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
asks which clades respond to the environment. Peres-Neto & Kembel (2015)
proposed two frameworks that ask whether phylogenetic composition
follows the environment across the whole tree, and both are in the
package. The results below are simulated.

``` r

library(comdynamite)

set.seed(1)
tree <- treeshapesim::sim_tree_shape(20)
traits <- species_traits(tree, lambda = 1, optima_range = c(-1, 101), rmax_range = c(0, 30))
env <- runif(50, 0, 100)
comm <- simulate_static(traits, env, incidence = TRUE)$abundance
keep <- rowSums(comm) > 0
comm <- comm[keep, ]
env <- env[keep]
```

## Phylogenetic community composition

[`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md)
takes the principal coordinates of the distances among species, averages
them within sites, and regresses each average on the predictors. The
test statistic is the variance of the adjusted R2 across eigenvectors,
and the value to report is the mean of the positive adjusted R2.

``` r

pcc <- pcc_test(comm, tree, env, n_perm = 199)
pcc
#> <mc_pcc> phylogenetic community composition
#> mean positive adjusted R2: 0.0965 
#> permutations: 199 free of sites, 199 of species; smallest possible p: 0.005
#> p (sites): 0.005  p (species): 0.005  p (max): 0.005
```

## Phylogenetic beta diversity

[`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
starts from the distance between communities given by
[`phylo_beta()`](https://provetelab.org/comdynamite/reference/phylo_beta.md),
the mean phylogenetic distance between their species, and fits a
distance-based redundancy analysis.

``` r

beta <- phylo_beta_test(comm, tree, env, n_perm = 199)
beta
#> <mc_phylo_beta> phylogenetic beta diversity
#> adjusted R2: 0.221  pseudo-F: 14.9 
#> permutations: 199 free of sites, 199 of species; smallest possible p: 0.005
#> p (sites): 0.005  p (species): 0.005  p (max): 0.005
```

`scores` holds the sites on the principal coordinates of the predicted
matrix, which the authors suggest using to order sites in a heatmap of
the community matrix.
[`plot_phylo_beta()`](https://provetelab.org/comdynamite/reference/plot_phylo_beta.md)
draws it, with species in the order of the tree, so that clades replace
each other along the sites when composition follows the environment.

``` r

plot_phylo_beta(beta, comm, tree)
```

![](whole-tree_files/figure-html/unnamed-chunk-5-1.png)

## Principal coordinates of phylogenetic structure

[`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md)
is the method of the package PCPS. Each species shares its relative
abundance with its relatives, in proportion to their similarity on the
tree, and the principal coordinates of the resulting matrix order sites
by phylogenetic composition. The test covers all the dissimilarity and,
separately, each axis in `axes`, whose `share` is its part of the
variation.

``` r

pcps <- pcps_test(comm, tree, env, n_perm = 199)
pcps
#> <mc_pcps> principal coordinates of phylogenetic structure
#> permutations: 199 free of sites, 199 of species; smallest possible p: 0.005
#>   test share     r2 statistic p_site p_species p_max
#>  whole 1.000 0.4989    47.793  0.005     0.025 0.025
#>  pcps1 0.541 0.9126   500.958  0.005     0.075 0.075
#>  pcps2 0.137 0.0133     0.648  0.445     0.735 0.735
```

Unlike the method of
[`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md),
which works on eigenvectors of the tree, this one takes eigenvectors of
the sites, so `scores` can be plotted against the environment.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the sites
on the first two axes, coloured by the environment, with each species as
an arrow from the correlation of its column of the matrix P with the
axes.

``` r

plot(pcps)
```

![](whole-tree_files/figure-html/unnamed-chunk-7-1.png)

When sites are neighbours and the environment and the communities are
both autocorrelated, free permutations of sites reject far too often.
`site_permutations` takes a restricted null model, such as
[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md)
on a complete grid, in
[`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md),
[`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
and
[`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md)
alike. The help of
[`pcps_test()`](https://provetelab.org/comdynamite/reference/pcps_test.md)
gives the rates, including those in a landscape with edges, where the
torus still rejects too often when fields are broad.

## Reading the two p-values

All these functions permute sites, which breaks the link between
composition and environment, and species across the tips, which breaks
the link between composition and phylogeny. The association requires
both links, so the p-value to report is `p_max`, the larger of the two.
A small `p_site` with a large `p_species` says that composition follows
the environment, but that the phylogeny adds nothing to it.

The scripts `data-raw/validation_pcc.R` and
`data-raw/validation_phylo_beta.R` compare the rejection rates of these
tests with those published.

## Sites that are neighbours

Permutations of sites treat them as exchangeable. When they are
neighbours, `site_permutations` takes a restricted null model, such as
those of
[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md)
for a grid on a torus.

## Reference

Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient
analysis: environmental drivers of phylogenetic variation across
ecological communities. *Plant Ecology*, 216, 709–724.
<https://doi.org/10.1007/s11258-014-0405-0>
