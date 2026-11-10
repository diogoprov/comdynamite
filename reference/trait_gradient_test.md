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
  null = c("species", "phylo", "bm", "curveball", "permute"),
  burn_in = NULL,
  thin = NULL,
  tree = NULL,
  phylo_k = 1.1
)
```

## Arguments

- x:

  Result of
  [`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md).

- n_null:

  Number of null communities.

- null:

  `"species"` (the default), `"phylo"`, `"bm"`, `"curveball"` or
  `"permute"`.

- burn_in, thin:

  Trades of the curveball chain before the first null community and
  between null communities; by default 50 and 5 times the number of
  sites, so that each site takes part in about 10 trades between null
  communities.

- tree:

  Phylogeny (`phylo`) with branch lengths, for `null = "phylo"` and
  `"bm"`; every species of `x` must be a tip.

- phylo_k:

  The constant \\k \ge 1\\ of the phylogenetic permutations of
  `null = "phylo"`: closer to 1, more exchanges among relatives only.

## Value

An object of class `mc_trait_gradient_test`, a list with `statistics`
(one row per statistic: `observed`, `null_mean`, `null_sd`, `ses` and
`p`), `species` (one row per species: observed \\\beta\\ with `ses_beta`
and `p_beta`, and with `null = "species"`, `"phylo"` or `"bm"` also
\\\alpha\\ with `ses_alpha` and `p_alpha`), `null` (a matrix with one
row per null community and one column per statistic) and the settings.

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

**`null = "permute"`**, the null of Ackerly & Cornwell (2007), permutes
the site of each population. It keeps the number of populations of each
site but can place a species twice in one site, and it breaks the
turnover of species as the curveball does.

**`null = "phylo"`.** The phylogenetic permutations of Lapointe &
Garland (2001): as `"species"`, but the mean of a species moves
preferably to a close relative. Patristic distances are scaled by the
largest, \\d^\*\_{ij} = d\_{ij} / d\_{max}\\, turned into similarities
\\s\_{ij} = k - d^\*\_{ij}\\ with \\k\\ = `phylo_k` \\\ge 1\\, and each
row is divided by its sum to give the probability \\p\_{ij}\\ that the
mean of species \\i\\ moves to species \\j\\ (staying is the most likely
move). Species are taken in random order, and each mean moves to a
species not yet occupied, with the probabilities of the remaining
species rescaled to sum to one; the last mean takes the last free
species. Large \\k\\ approaches `"species"` (equally likely
permutations); on an ultrametric tree \\k = 1\\ gives no exchange across
the root. Lapointe & Garland (2001) illustrate \\k\\ from 1.0 to 2. A
species with no close relative keeps its own mean in many nulls (about
half of them in our trials with \\k = 1.1\\), so the test has little
power for such species.

**`null = "bm"`.** The species means are replaced by values simulated by
Brownian motion on `tree`, with the mean and rate estimated from the
observed means by restricted maximum likelihood; each species keeps its
own variation among sites. The null is a trait that evolves along the
tree and has nothing to do with assembly.

Use `"phylo"` when the niches of species are conserved: under
`"species"`, a trait that merely shares the phylogenetic structure of
the niches looks sorted along the gradient. `"bm"` is exact only when
the trait does evolve by Brownian motion: a trait with less phylogenetic
structure than Brownian motion makes its variance tests reject (see
Validation).

The \\\alpha\\ of a species is not tested under `"curveball"` and
`"permute"`: the species mean does not change, so its test would mirror
that of \\\beta\\.

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

In the same design (150 replicates per scenario), for the three
statistics:

|  |  |  |  |
|----|----|----|----|
| null | trait independent of niches | trait and niches on the tree | sorting |
| `"species"` | 4.0–4.7% | 14–23% | 99–100% |
| `"phylo"`, \\k = 1.1\\ | 2.7–4.0% | 4.0–11% | 95–100% |
| `"phylo"`, \\k = 1.01\\ | 1.3–4.7% | 2.0–10% | 91–100% |
| `"bm"` | 10–85% | 2.7–4.7% | 43–100% |

The phylogenetic permutations bring the error back near the nominal
level in both null scenarios, with some excess left for the correlation
and the variance of \\\alpha\\ when trait and niches share the tree;
`"bm"` is exact in that scenario but fails when the trait has no
phylogenetic structure, and loses power for the variance of \\\beta\\.

## References

Ackerly, D. D., & Cornwell, W. K. (2007). A trait-based approach to
community assembly: partitioning of species trait values into within-
and among-community components. *Ecology Letters*, 10, 135–145.
[doi:10.1111/j.1461-0248.2006.01006.x](https://doi.org/10.1111/j.1461-0248.2006.01006.x)

Strona, G., Nappo, D., Boccacci, F., Fattorini, S., & San-Miguel-Ayanz,
J. (2014). A fast and unbiased procedure to randomize ecological binary
matrices with fixed row and column totals. *Nature Communications*, 5,
4114. [doi:10.1038/ncomms5114](https://doi.org/10.1038/ncomms5114)

Lapointe, F.-J., & Garland, T., Jr. (2001). A generalized permutation
model for the analysis of cross-species data. *Journal of
Classification*, 18, 109–127.
[doi:10.1007/s00357-001-0007-0](https://doi.org/10.1007/s00357-001-0007-0)

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
[`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md),
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
# means exchanged among relatives, or simulated along the tree
trait_gradient_test(tg, n_null = 99, null = "phylo", tree = tree)
#> Null model of the trait-gradient analysis: 99 communities, phylo 
#>       statistic observed null_mean null_sd    ses    p
#>  cor_beta_alpha    0.697    -0.401   0.319  3.440 0.02
#>        var_beta   68.200    41.200  41.500  0.651 0.30
#>       var_alpha  200.000   505.000  99.200 -3.080 0.02
trait_gradient_test(tg, n_null = 99, null = "bm", tree = tree)
#> Null model of the trait-gradient analysis: 99 communities, bm 
#>       statistic observed null_mean null_sd      ses    p
#>  cor_beta_alpha    0.697    -0.472   0.196  5.97000 0.02
#>        var_beta   68.200    67.800  56.600  0.00652 0.68
#>       var_alpha  200.000   878.000 452.000 -1.50000 0.02
```
