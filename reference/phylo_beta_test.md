# Test phylogenetic beta diversity against the environment

Implements the phylogenetic beta-diversity framework of Peres-Neto &
Kembel (2015), a distance-based redundancy analysis (Legendre & Anderson
1999; McArdle & Anderson 2001) of
[`phylo_beta()`](https://provetelab.org/comdynamite/reference/phylo_beta.md)
on the predictors. The test statistic is the pseudo-F. It is tested by
permuting sites, which breaks the link between composition and
environment, and by permuting species across the tips, which breaks the
link between composition and phylogeny. The association requires both,
so the p-value to report is the larger of the two.

## Usage

``` r
phylo_beta_test(comm, tree, env, n_perm = 999, site_permutations = NULL)
```

## Arguments

- comm:

  Matrix with one row per site and one column per species, named as the
  tips of `tree`. Every site needs at least one individual. Or an
  [`mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md)
  object, whose parts fill the arguments left out.

- tree:

  A `phylo` object with the same species as `comm`.

- env:

  Predictor, or matrix or data frame of predictors, with one row per
  site.

- n_perm:

  Number of permutations of sites and of species.

- site_permutations:

  Optional matrix of site indices with one permutation per row, for
  example from
  [`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md),
  to use instead of free permutations of sites. Empty sites are then
  allowed, because they are part of the landscape: they leave the
  analysis but their predictors can reach occupied sites.

## Value

List of class `mc_phylo_beta` with `r2_adjusted`, the adjusted R2;
`statistic`, the pseudo-F; `p_site`, `p_species` and `p_max`; `scores`,
the sites on the principal coordinates of the predicted matrix, which
the authors suggest using to order sites in a heatmap of the community
matrix; and `permutations`.

## Details

The equations printed in the article omit the 1 in the adjusted R2, the
transpose in the hat matrix and the centring of the matrix that is
regressed. Here the predictors and the matrix are centred, and the
residual degrees of freedom are \\n - p - 1\\, with \\n - 1\\ for the
total, instead of the printed \\n - p\\ and \\n\\. The pseudo-F and R2
therefore equal those of a distance-based redundancy analysis in vegan.
The printed degrees of freedom would multiply the pseudo-F by a
constant, which does not change the permutation tests.

The authors used presences and absences. Abundances are accepted, as
weights in
[`phylo_beta()`](https://provetelab.org/comdynamite/reference/phylo_beta.md).

## References

Legendre, P., & Anderson, M. J. (1999). Distance-based redundancy
analysis: testing multispecies responses in multifactorial ecological
experiments. *Ecological Monographs*, 69, 1–24.
[doi:10.1890/0012-9615(1999)069\[0001:DBRATM\]2.0.CO;2](https://doi.org/10.1890/0012-9615%281999%29069%5B0001%3ADBRATM%5D2.0.CO%3B2)

McArdle, B. H., & Anderson, M. J. (2001). Fitting multivariate models to
community data: a comment on distance-based redundancy analysis.
*Ecology*, 82, 290–297.
[doi:10.1890/0012-9658(2001)082\[0290:FMMTCD\]2.0.CO;2](https://doi.org/10.1890/0012-9658%282001%29082%5B0290%3AFMMTCD%5D2.0.CO%3B2)

Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient
analysis: environmental drivers of phylogenetic variation across
ecological communities. *Plant Ecology*, 216, 709–724.
[doi:10.1007/s11258-014-0405-0](https://doi.org/10.1007/s11258-014-0405-0)

## See also

[`plot_phylo_beta()`](https://provetelab.org/comdynamite/reference/plot_phylo_beta.md)
shows the result.
[`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md)
is the other framework of the same article.

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(20)
traits <- species_traits(tree, lambda = 1, optima_range = c(-1, 101),
                         rmax_range = c(0, 30))
env <- runif(50, 0, 100)
comm <- simulate_static(traits, env, incidence = TRUE)$abundance
keep <- rowSums(comm) > 0
phylo_beta_test(comm[keep, ], tree, env[keep], n_perm = 99)
#> <mc_phylo_beta> phylogenetic beta diversity
#> adjusted R2: 0.221  pseudo-F: 14.9 
#> permutations: 99 free of sites, 99 of species; smallest possible p: 0.01
#> p (sites): 0.01  p (species): 0.01  p (max): 0.01 
```
