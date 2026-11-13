# Phylogenetic models of beta and alpha, with their sampling error

Fits models of trait evolution to the \\\beta\\, \\\alpha\\ and mean
trait of the species from
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md),
with the sampling error of each value added to the phylogenetic
covariance, and compares a model without phylogeny, Brownian motion and
an Ornstein-Uhlenbeck process by AIC. It tells whether the position of
species along the gradient (\\\beta\\) and their difference from
co-occurring species (\\\alpha\\) are conserved or labile.

## Usage

``` r
trait_gradient_phylo(
  x,
  tree,
  responses = c("beta", "alpha", "mean"),
  models = c("white", "BM", "OU"),
  n_boot = 199,
  min_sites = 3
)
```

## Arguments

- x:

  Result of
  [`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md).

- tree:

  Phylogeny (`phylo`) with branch lengths; its tips are matched to
  species by name.

- responses:

  Which of `"beta"`, `"alpha"` and `"mean"` to model.

- models:

  Which of `"white"`, `"BM"` and `"OU"` to fit.

- n_boot:

  Number of bootstrap samples for the standard errors.

- min_sites:

  Smallest number of sites of a species with a value (sites where it is
  alone do not count); species with fewer sites, or not in `tree`, are
  left out.

## Value

An object of class `mc_trait_gradient_phylo`, a list with `models` (one
row per response and model: `response`, `model`, `mu`, `sigma2`, `a`,
`half_life` (\\\ln 2 / a\\, in units of branch length), `at_bound` (an
OU fit with \\a\\ at a bound), `loglik`, `n_par`, `aic` and `delta_aic`
within each response), `data` (one row per species used: values and
standard errors) and `tree` (the tree pruned to those species).

## Details

The standard error of each value comes from a bootstrap of the sites of
each species: the sites (with their trait values and abundances) are
resampled with replacement, and \\\bar t_i\\, \\\beta_i\\ and
\\\alpha_i\\ recomputed with the site means of the observed data
(without the focal species, if so chosen). The site means are held fixed
because the question is how well the sites of a species measure its
position, not how the gradient itself would change.

For each response \\y\\ the model is \\y \sim N(\mu 1, \sigma^2 C +
S)\\, with \\S\\ the diagonal matrix of squared standard errors and
\\C\\ the identity (`"white"`), the shared branch lengths (`"BM"`), or,
for `"OU"`, the covariance of an Ornstein-Uhlenbeck process with a root
at \\\mu\\ and strength \\a\\, \$\$C\_{ij} = e^{-a(d_i + d_j - 2
s\_{ij})} \frac{1 - e^{-2 a s\_{ij}}}{2a},\$\$ where \\d_i\\ is the
distance from the root to tip \\i\\ and \\s\_{ij}\\ the length shared by
\\i\\ and \\j\\. The parameters are estimated by restricted maximum
likelihood, so AIC compares the covariance models with the same mean.
Without the error term, the noise of species seen in few sites passes
for evolution near the tips, which an OU process fits better: in
simulations of \\\beta\\ under Brownian motion (100 replicates, 30
species, 25 sites), AIC chose OU in 22 without the error term and in 9
with it.

\\a\\ is bounded between \\10^{-3}\\ and 50 divided by the height of the
tree. At the lower bound OU is Brownian motion and at the upper bound it
is the model without phylogeny; such fits are marked in `at_bound` and
their half-life is `NA`. AIC differences among three models fitted to
few species are unstable, and the bootstrap standard error of a species
seen in \\n\\ sites is biased low by about \\\sqrt{(n - 1)/n}\\.

## References

Ackerly, D. D., & Cornwell, W. K. (2007). A trait-based approach to
community assembly: partitioning of species trait values into within-
and among-community components. *Ecology Letters*, 10, 135–145.
[doi:10.1111/j.1461-0248.2006.01006.x](https://doi.org/10.1111/j.1461-0248.2006.01006.x)

## See also

[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md),
[`trait_ou()`](https://provetelab.org/comdynamite/reference/trait_ou.md)
to simulate traits under these models.

Other trait variation among individuals:
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
[`drift_test()`](https://provetelab.org/comdynamite/reference/drift_test.md),
[`neighbour_metrics`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
[`nested_variance()`](https://provetelab.org/comdynamite/reference/nested_variance.md),
[`null_intervals()`](https://provetelab.org/comdynamite/reference/null_intervals.md),
[`plot.mc_trait_gradient()`](https://provetelab.org/comdynamite/reference/plot.mc_trait_gradient.md),
[`plot.nested_variance()`](https://provetelab.org/comdynamite/reference/plot.nested_variance.md),
[`plot.trait_flex()`](https://provetelab.org/comdynamite/reference/plot.trait_flex.md),
[`plot.trait_rarefaction()`](https://provetelab.org/comdynamite/reference/plot.trait_rarefaction.md),
[`plot.tstats()`](https://provetelab.org/comdynamite/reference/plot.tstats.md),
[`plot_population_traits()`](https://provetelab.org/comdynamite/reference/plot_population_traits.md),
[`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md),
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md),
[`subsample_individuals()`](https://provetelab.org/comdynamite/reference/subsample_individuals.md),
[`trait_densities()`](https://provetelab.org/comdynamite/reference/trait_densities.md),
[`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md),
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md),
[`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md),
[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md),
[`trait_gradient_weights()`](https://provetelab.org/comdynamite/reference/trait_gradient_weights.md),
[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
[`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md)

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(20)
sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree, lambda = 1),
                              competition_matrix(tree, "none"), n_steps = 100,
                              seed_steps = 1:5)
ind <- sample_individuals(sim, n = 5, local_adaptation = 0.5, intraspecific_sd = 1)
tg <- trait_gradient(ind, abundance = "abundance")
trait_gradient_phylo(tg, tree, n_boot = 99)
#> Phylogenetic models of the trait-gradient analysis: 14 species
#>  response model     mu sigma2     a half_life at_bound loglik n_par   aic
#>      beta white 64.600   90.9    NA        NA       NA  -49.9     2 104.0
#>      beta    BM 61.700  201.0    NA        NA       NA  -49.7     2 103.0
#>      beta    OU 62.900  450.0 1.820      0.38    FALSE  -49.2     3 104.0
#>     alpha white  0.605   38.7    NA        NA       NA  -44.0     2  92.0
#>     alpha    BM -1.000   75.6    NA        NA       NA  -41.4     2  86.9
#>     alpha    OU -1.000   75.7 0.001        NA     TRUE  -41.4     3  88.9
#>      mean white 64.600  236.0    NA        NA       NA  -55.4     2 115.0
#>      mean    BM 60.000  538.0    NA        NA       NA  -52.9     2 110.0
#>      mean    OU 60.600  722.0 0.692      1.00    FALSE  -52.8     3 112.0
#>  delta_aic
#>       0.40
#>       0.00
#>       1.09
#>       5.08
#>       0.00
#>       2.00
#>       4.92
#>       0.00
#>       1.72
```
