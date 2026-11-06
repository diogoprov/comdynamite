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
`eigenvalues`; and `eigenvector_correlation`.

## Details

The p-values of the global test are the ones to report. The two vectors
are chosen to fit the predictor, so a test of them with
[`fourth_corner_test()`](https://provetelab.org/comdynamite/reference/fourth_corner_test.md)
is not valid on its own. They are meant for
[`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md)
after a significant global test.

Species are matched to the tips of `tree` by name, in any order, and
those in only one of the two are dropped with a warning. Species without
individuals have no weight and are left out.

## Article and original code

The code and the text differ. Here we follow the text in the first three
points and the code in the last:

- The text sums eigenvectors weighted by their correlations. The code
  uses the absolute correlations, and eigenvectors with opposite signs
  then cancel.

- The text uses the square root of the path length as distance. The code
  example uses one minus the covariance.

- The text draws 1000 Brownian vectors, and the code 100.

- The text multiplies each Brownian vector by `pev`, and the code adds
  them.

In addition, `PCoAWeigthed` of the authors is called with weights that
do not sum to one and returns coordinates without dividing by the square
root of the weights, so its centring is not the one in the text.
`DoubleCenterWeigth`, which is not used in their example, does both, and
agrees with this function.

## References

Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the
effects of spatial contingency and environmental filtering on
metacommunity phylogenetics. *Ecology*, 93, S14–S30.
[doi:10.1890/11-0494.1](https://doi.org/10.1890/11-0494.1)

## Examples

``` r
set.seed(1)
tree <- ape::rcoal(20)
traits <- species_traits(tree, lambda = 1, optima_range = c(0, 100),
                         rmax_range = c(0, 30))
env <- runif(50, 0, 100)
comm <- simulate_static(traits, env, incidence = TRUE)$abundance
keep <- rowSums(comm) > 0
phylo_vector(comm[keep, ], tree, env[keep], n_perm = 99, n_sim = 100)
#> <mc_phylo_vector> 18 eigenvectors, 1 retained
#> global test: p (sites): 0.01  p (species): 0.01  p (max): 0.01 
#> correlation: PEV 0.513  PEV-B 0.602 
```
