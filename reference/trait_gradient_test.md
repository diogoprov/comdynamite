# Null models of the trait-gradient analysis

Compares the values of
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md)
with those of null communities. The default null (`"species"`) keeps
every community as observed and exchanges the trait means of species
among species, so it asks whether the traits of species follow the sites
they share. The curveball null (`"curveball"`, Strona et al. 2014) keeps
the traits and moves species among sites, with the richness of each site
and the number of sites of each species fixed, so it asks whether the
observed co-occurrences matter given the traits.

## Usage

``` r
trait_gradient_test(
  x,
  n_null = 999,
  null = c("species", "curveball", "permute"),
  burn_in = NULL,
  thin = NULL
)
```

## Arguments

- x:

  Result of
  [`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md).

- n_null:

  Number of null communities.

- null:

  `"species"` (the default), `"curveball"` or `"permute"`.

- burn_in, thin:

  Trades of the curveball chain before the first null community and
  between null communities; by default 50 and 5 times the number of
  sites, so that each site takes part in about 10 trades between null
  communities.

## Value

An object of class `mc_trait_gradient_test`, a list with `statistics`
(one row per statistic: `observed`, `null_mean`, `null_sd`, `ses` and
`p`), `species` (one row per species: observed \\\beta\\ with `ses_beta`
and `p_beta`, and with `null = "species"` also \\\alpha\\ with
`ses_alpha` and `p_alpha`), `null` (a matrix with one row per null
community and one column per statistic) and the settings.

## Details

Statistics: the correlation of \\\beta\\ and \\\alpha\\, the variance of
\\\beta\\ (how far species spread along the gradient) and the variance
of \\\alpha\\, and the \\\beta\\ of each species. Tests are two-sided.

**`null = "species"`.** The trait values of species \\i\\ become
\\t\_{ij} - \bar t_i + \bar t\_{\pi(i)}\\ for a random permutation
\\\pi\\ of species: each species carries the mean of another and keeps
its own variation among sites. Sites, abundances and co-occurrences stay
as observed. This is the null for the variances: a large variance of
\\\beta\\ means that species with similar traits live in similar
communities.

**`null = "curveball"`.** Populations move among sites with their trait
values and abundances by curveball trades; the chain starts at the
observed communities and runs `burn_in` trades and then `thin` trades
between null communities. Because it also breaks the turnover of species
along the gradient, the variances of \\\beta\\ and \\\alpha\\ differ
from the null whenever composition changes among sites, even when traits
are unrelated to it: species that share sites share their site means.
Use the curveball for the correlation only.

## Validation

In simulated metacommunities (25 species, 20 sites, 150 replicates per
scenario, 5% level, leave-one-out site means) with trait means drawn
independently of the niches, the species null rejected in 4.7% (variance
of \\\beta\\), 5.3% (variance of \\\alpha\\) and 7.3% (correlation) of
replicates; the curveball rejected in 5.3% for the correlation but in
26–66% for the variances. With traits that sort species along the
gradient, all tests rejected in 98–100%. When the trait and the niches
both evolve along the tree, the species null rejected in 10–17% although
the trait plays no part in assembly: related species share both sites
and trait values, and the test detects that association, not its cause.

**`null = "permute"`**, the null of Ackerly & Cornwell (2007), permutes
the site of each population. It keeps the number of populations of each
site but can place a species twice in one site, and it breaks the
turnover of species as the curveball does.

The \\\alpha\\ of a species is not tested under `"curveball"` and
`"permute"`: the species mean does not change, so its test would mirror
that of \\\beta\\.

## References

Ackerly, D. D., & Cornwell, W. K. (2007). A trait-based approach to
community assembly: partitioning of species trait values into within-
and among-community components. *Ecology Letters*, 10, 135–145.
[doi:10.1111/j.1461-0248.2006.01006.x](https://doi.org/10.1111/j.1461-0248.2006.01006.x)

Strona, G., Nappo, D., Boccacci, F., Fattorini, S., & San-Miguel-Ayanz,
J. (2014). A fast and unbiased procedure to randomize ecological binary
matrices with fixed row and column totals. *Nature Communications*, 5,
4114. [doi:10.1038/ncomms5114](https://doi.org/10.1038/ncomms5114)

## See also

[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md)

Other trait variation among individuals:
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
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
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md),
[`trait_gradient_weights()`](https://provetelab.org/comdynamite/reference/trait_gradient_weights.md),
[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
[`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md)

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(20)
sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
                              competition_matrix(tree, "none"), n_steps = 100,
                              seed_steps = 1:5)
ind <- sample_individuals(sim, n = 5, intraspecific_sd = 1)
tg <- trait_gradient(ind, abundance = "abundance")
trait_gradient_test(tg, n_null = 99)
#> Null model of the trait-gradient analysis: 99 communities, species 
#>       statistic observed null_mean null_sd    ses    p
#>  cor_beta_alpha    0.697    -0.417    0.25  4.460 0.02
#>        var_beta   68.200    35.300   43.90  0.749 0.32
#>       var_alpha  200.000   504.000   90.60 -3.350 0.02
# the correlation of beta and alpha against the curveball
trait_gradient_test(tg, n_null = 99, null = "curveball")$statistics[1, ]
#>        statistic  observed null_mean   null_sd      ses    p
#> 1 cor_beta_alpha 0.6969889 -0.500578 0.2130449 5.621195 0.02
```
