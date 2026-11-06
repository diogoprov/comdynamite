# Test clades against the environment, node by node

Implements the regression of Leibold et al. (2010). For each node, sites
without species of that node are removed, the contrast is standardised
with species weights and summed within sites, and the result is
regressed on the predictors standardised with site weights, using site
richness as weights. Each node shows whether the two clades that split
there differ in their response to the predictors.

## Usage

``` r
node_test(comm, nodes, env, biogeography = NULL, n_perm = 9999)
```

## Arguments

- comm:

  Matrix with one row per site and one column per species. Species are
  matched to the rows of `nodes` by name, in any order, and those in
  only one of the two are dropped with a warning.

- nodes:

  Matrix of contrasts with one row per species, from
  [`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md).

- env:

  Predictor, or matrix or data frame of predictors, with one row per
  site.

- biogeography:

  Optional matrix of biogeographic events with one row per site and
  values of 1, -1 and 0.

- n_perm:

  Number of permutations of sites and of species.

## Value

List of class `mc_node_test` with `nodes`, a data frame with one row per
node with `node`, `n_species`, `n_sites`, `r2`, `r2_adjusted`,
`statistic`, `p_site`, `p_taxa`, `p_max` and `p_adjusted`;
`semipartial`, a matrix of nodes by predictors; `p_partial`, with the
larger of the two p-values for each of them; and `biogeography`, a list
with the same two matrices for the events, or `NULL`.

## Details

The test statistic is the explained over the residual sum of squares. It
is tested by permuting sites, and by replacing the species of the node
with a random set of the same size from all species. The article
permutes both at once; following Peres-Neto et al. (2012) and ter Braak
et al. (2012), the p-value here is the larger of the two tests. With
many nodes these p-values need a correction, and `p_adjusted` gives the
false discovery rate that the article uses. This correction needs many
permutations, as the smallest p-value is 1 / (`n_perm` + 1), so the
default is the 9999 of the article.

Each predictor gets a semipartial correlation, the square root of the
loss of R2 when it is removed, with the sign of its slope when it is
alone in the model.

Biogeographic events follow the article: a column of `biogeography` has
1 for sites in one area, -1 for sites in the other and 0 for sites in
neither. Each event is tested only in its sites, controlling for `env`.
Other contrasts are allowed, for example a column without zeros, which
compares one area with all the others. Events with few sites in one area
give unstable values.

R2, adjusted R2 and the statistic equal those of the MATLAB function
`MetacommunityPhylogenetics` of the authors, checked in Octave. Three
parts of that function differ from this one:

- In the original, the species permutations shuffle only the first
  columns of the community matrix, as many as there are species in the
  node, so they do not draw from all species. Here they do.

- The original keeps the total of the observed data after permuting
  species. Here it is recalculated.

- The semipartial values of the original use the raw predictors, without
  an intercept, and multiply by the absolute slope. The article uses
  standardised predictors and divides by it, to keep only the sign,
  which is what we do.

As in the original, R2 sums squares without site weights, although the
fit uses them, so it is not the R2 of a weighted regression, and the
adjusted R2 uses the degrees of freedom of an unweighted model.

Nodes for which a predictor has no variance among the remaining sites,
or with too few sites, get `NA`. The original stops with an error.

## References

Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
phylogenetics: separating the roles of environmental filters and
historical biogeography. *Ecology Letters*, 13, 1290–1299.
[doi:10.1111/j.1461-0248.2010.01523.x](https://doi.org/10.1111/j.1461-0248.2010.01523.x)

Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the
effects of spatial contingency and environmental filtering on
metacommunity phylogenetics. *Ecology*, 93, S14–S30.
[doi:10.1890/11-0494.1](https://doi.org/10.1890/11-0494.1)

ter Braak, C. J. F., Cormont, A., & Dray, S. (2012). Improved testing of
species traits-environment relationships in the fourth-corner problem.
*Ecology*, 93, 1525–1526.

## See also

[`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md)
and
[`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md)
test the whole tree at once.

## Examples

``` r
set.seed(1)
tree <- ape::rcoal(12)
traits <- species_traits(tree, lambda = 1, optima_range = c(-1, 101),
                         rmax_range = c(0, 30))
env <- runif(40, 0, 100)
comm <- simulate_static(traits, env, incidence = TRUE)$abundance
keep <- rowSums(comm) > 0
node_test(comm[keep, ], node_contrasts(tree), env[keep], n_perm = 19)
#> <mc_node_test> 7 nodes, 7 tested
#> nodes with p (max) <= 0.05: 0 before correction, 0 after
```
