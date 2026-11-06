# Test phylogenetic community composition against the environment

Implements the phylogenetic community composition (PCC) framework of
Peres-Neto & Kembel (2015). Phylogenetic eigenvectors are standardised
with species weights and averaged within sites, and each resulting
variable is regressed on the predictors with site weights. The test
statistic is the variance of the adjusted R2 across eigenvectors. It is
tested by permuting sites, which breaks the link between composition and
environment, and by permuting species across the tips, which breaks the
link between composition and phylogeny. The association requires both,
so the p-value to report is the larger of the two.

## Usage

``` r
pcc_test(
  comm,
  tree,
  env,
  n_perm = 999,
  squared = FALSE,
  site_permutations = NULL
)
```

## Arguments

- comm:

  Matrix with one row per site and one column per species. Species are
  matched to the tips of `tree` by name, in any order, and those in only
  one of the two are dropped with a warning. Every site needs at least
  one individual.

- tree:

  A `phylo` object.

- env:

  Predictor, or matrix or data frame of predictors, with one row per
  site.

- n_perm:

  Number of permutations of sites and of species.

- squared:

  Passed to
  [`phylo_eigenvectors()`](https://provetelab.org/comdynamite/reference/phylo_eigenvectors.md).

- site_permutations:

  Optional matrix of site indices with one permutation per row, for
  example from
  [`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md),
  to use instead of free permutations of sites.

## Value

List of class `mc_pcc` with `r2`, the mean of the positive adjusted R2
values, which is the value the authors recommend reporting; `statistic`,
the variance of the adjusted R2 across eigenvectors; `p_site`, `p_taxa`
and `p_max`; and `r2_adjusted`, with one value per eigenvector.

## Details

The authors used presences and absences. Abundances are accepted, with
species and site totals as weights.

## References

Peres-Neto, P. R., & Kembel, S. W. (2015). Phylogenetic gradient
analysis: environmental drivers of phylogenetic variation across
ecological communities. *Plant Ecology*, 216, 709–724.
[doi:10.1007/s11258-014-0405-0](https://doi.org/10.1007/s11258-014-0405-0)

## Examples

``` r
set.seed(1)
tree <- ape::rcoal(20)
traits <- species_traits(tree, lambda = 1, optima_range = c(-1, 101),
                         rmax_range = c(0, 30))
env <- runif(50, 0, 100)
comm <- simulate_static(traits, env, incidence = TRUE)$abundance
keep <- rowSums(comm) > 0
pcc_test(comm[keep, ], tree, env[keep], n_perm = 99)
#> <mc_pcc> phylogenetic community composition
#> mean positive adjusted R2: 0.0893 
#> p (sites): 0.01  p (species): 0.01  p (max): 0.01 
```
