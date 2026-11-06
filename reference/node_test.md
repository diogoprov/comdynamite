# Test clades against the environment, node by node

Implements the regression of Leibold et al. (2010). For each node, sites
without species of that node are removed, the contrast is standardised
with species weights and summed within sites, and the result is
regressed on the predictors standardised with site weights, using site
richness as weights. Each node shows whether the two clades that split
there differ in their response to the predictors.

## Usage

``` r
node_test(comm, nodes, env, n_perm = 999)
```

## Arguments

- comm:

  Matrix with one row per site and one column per species, named as the
  rows of `nodes`.

- nodes:

  Matrix of contrasts with one row per species, from
  [`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md).

- env:

  Predictor, or matrix or data frame of predictors, with one row per
  site.

- n_perm:

  Number of permutations of sites and of species.

## Value

List of class `mc_node_test` with `nodes`, a data frame with one row per
node with `node`, `n_species`, `n_sites`, `r2`, `r2_adjusted`,
`statistic`, `p_site`, `p_taxa` and `p_max`; `partial_r2`, a matrix of
nodes by predictors; and `p_partial`, with the larger of the two
p-values for each of them.

## Details

The test statistic is the explained over the residual sum of squares. It
is tested by permuting sites, and by replacing the species of the node
with a random set of the same size from all species. Following
Peres-Neto et al. (2012) and ter Braak et al. (2012), the p-value to
report is the larger of the two. With many nodes, these p-values need a
correction for multiple tests.

R2, adjusted R2 and the statistic equal those of the MATLAB function
`MetacommunityPhylogenetics` of the authors, checked in Octave. Three
parts differ:

- In the original, the species permutations shuffle only the first
  columns of the community matrix, as many as there are species in the
  node, so they do not draw from all species. Here they do.

- The original keeps the total of the observed data after permuting
  species. Here it is recalculated.

- The semipartial values of the original use the raw predictors, without
  an intercept, whereas the full model uses standardised predictors.
  Here `partial_r2` is the loss of R2 when a predictor is removed, with
  standardised predictors in both models.

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
#> nodes with p (max) <= 0.05, before correction: 0 
```
