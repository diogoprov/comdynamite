# Represent a phylogeny by one vector related to a predictor

The phylogenetic representation of Peres-Neto et al. (2012).
Eigenvectors of the phylogenetic distances are extracted with species
weights. A global test relates all of them to the predictor. Only if it
is significant, eigenvectors enter in the order of their correlations
until the first that is not significant, and they are summed with their
correlations as weights (`pev`). Brownian vectors are then added to
`pev`, and the sum with the largest correlation is kept (`pev_brown`).

## Usage

``` r
phylo_vector(
  comm,
  tree,
  predictor,
  n_perm = 999,
  alpha = 0.05,
  n_sim = 1000,
  site_permutations = NULL
)
```

## Arguments

- comm:

  Matrix with one row per site and one column per species. Every site
  needs at least one individual.

- tree:

  A `phylo` object.

- predictor:

  Numeric vector with the value of each site. The authors estimate one
  vector for each predictor.

- n_perm:

  Number of permutations of sites and of species.

- alpha:

  Significance level of the global test and of the selection.

- n_sim:

  Number of Brownian vectors.

- site_permutations:

  Optional matrix of site indices with one permutation per row, for
  example from
  [`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md),
  to use instead of free permutations.

## Value

List of class `mc_phylo_vector` with `correlation`, the correlations of
`pev` and `pev_brown`; `p_site`, `p_taxa` and `p_max` of the global
test; `selected`, the columns of `eigenvectors` retained; `pev` and
`pev_brown`, named vectors or `NULL` when the global test is not
significant or no eigenvector is retained; `eigenvectors`;
`eigenvalues`; `eigenvector_correlation`; and `permutations`. Check
`pev` before combining gradients with
[`cbind()`](https://rdrr.io/r/base/cbind.html), which drops a `NULL`
without warning.

## Details

The correlations are fourth-corner correlations, and the p-value to
report is the larger of the tests that permute sites and species (ter
Braak et al. 2012). For a test between a trait and the environment, see
[`ade4::fourthcorner()`](https://adeverse.github.io/ade4/reference/fourthcorner.html).
The two vectors are chosen to fit the predictor, so a further test of
them is not valid on its own. They are meant for
[`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md)
after a significant global test.

Species are matched to the tips of `tree` by name, in any order, and
those in only one of the two are dropped with a warning. Species without
individuals have no weight and are left out.

## Note

Adapted from the MATLAB functions `TraitEnvCor`,
`ModelSelectionStochastic` and `PCoAWeigthed` in the supplement of
Peres-Neto et al. (2012).

## Article and original code

The functions `TraitEnvCor`, `ModelSelectionStochastic` and
`PCoAWeigthed` and the script `ScriptDemo` of the supplement differ from
the article and its appendices in the points below. Here we follow the
text, except in the addition of Brownian vectors and where noted.

- Weights of the eigenvectors: the text uses the correlations and the
  code their absolute values, so eigenvectors with opposite signs
  cancel.

- Brownian vectors: the text multiplies each one by `pev` and the code
  adds them. We follow the code, as the product of a random vector with
  mean zero would not keep the structure of `pev`.

- Number of Brownian vectors: 1000 in the text and 100 in the code.

- Distances: the text uses the square root of the path length, and the
  script one minus the covariance.

- Centring: `Q` in the text only centres if the weights sum to one, and
  the script passes the sums of the species. `PCoAWeigthed` also does
  not divide the coordinates by the square root of the weights.
  `DoubleCenterWeigth`, which the script does not call, does both and
  agrees with this function.

- Global test: the text regresses each eigenvector on all predictors,
  and the code accepts one. With one predictor both are the same
  statistic. We also take one predictor.

- Forward selection: the text tests each eigenvector given the remaining
  predictors, and the code tests its simple correlation. We follow the
  code.

- Significance level: Appendix B uses 0.05, and the code uses alpha / 2
  for every permutation strategy, although Appendix A gives alpha / 2
  only for the stacked test. Here the larger p-value is compared with
  `alpha`.

- p-value of the global test with row permutations: the code divides by
  twice the number of permutations, which halves it.

- No significant result: the code returns random normal vectors, which
  can then enter later analyses. Here `pev` and `pev_brown` are `NULL`.

- Choice of vector: the article recommends the one with the highest R2,
  and the script always takes the Brownian vector. Here both are
  returned.

- Combination of permutations: Appendix A stacks the two null
  distributions, adds one to the count and uses alpha / 2, and
  `TraitEnvCor` adds two. The stacked p-value is close to the mean of
  the two p-values, so any test that it rejects at alpha / 2 is also
  rejected by the larger p-value at alpha, which we use.

- Permuted variables: Appendix A recalculates the statistic, and
  `TraitEnvCor` keeps the observed mean and standard deviation, so a
  permuted variable no longer has weighted mean zero. Here it is
  standardised again, as in the code of the partition.

- The code runs `n_perm` - 1 permutations in each null distribution.

- The supplement calls `eigen`, which it does not include, and the
  example 2 of `ScriptDemo` removes rows of `PhyloDistance` before
  creating it.

## References

Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the
effects of spatial contingency and environmental filtering on
metacommunity phylogenetics. *Ecology*, 93, S14–S30.
[doi:10.1890/11-0494.1](https://doi.org/10.1890/11-0494.1)

ter Braak, C. J. F., Cormont, A., & Dray, S. (2012). Improved testing of
species traits–environment relationships in the fourth-corner problem.
*Ecology*, 93, 1525–1526.
[doi:10.1890/12-0126.1](https://doi.org/10.1890/12-0126.1)

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(20)
traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100),
                         rmax_range = c(0, 30))
env <- runif(50, 0, 100)
comm <- simulate_static(traits, env, incidence = TRUE)$abundance
keep <- rowSums(comm) > 0
phylo_vector(comm[keep, ], tree, env[keep], n_perm = 99, n_sim = 100)
#> <mc_phylo_vector> 18 eigenvectors, 2 retained
#> permutations: 99 free of sites, 99 of species; smallest possible p: 0.01 
#> global test: p (sites): 0.01  p (species): 0.01  p (max): 0.01 
#> correlation: PEV 0.61  PEV-B 0.588 
```
