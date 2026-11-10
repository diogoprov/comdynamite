# Trait-gradient analysis

Splits the trait value of each species into a position along the
gradient of community trait means (\\\beta\\) and a difference from the
species it co-occurs with (\\\alpha\\), as in the trait-gradient
analysis of Ackerly & Cornwell (2007), with three changes: the mean of
each site leaves out the focal species, the weights of the site means
and of the species means are chosen separately, and
[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md)
and
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md)
give a null model that keeps the richness of sites and the frequency of
species and a phylogenetic model that carries the error of \\\beta\\ and
\\\alpha\\.

## Usage

``` r
trait_gradient(
  data,
  trait = "trait",
  site = "site",
  species = "species",
  abundance = NULL,
  site_weights = if (is.null(abundance)) "occurrence" else "abundance",
  species_weights = site_weights,
  leave_one_out = TRUE,
  min_slope = 3
)
```

## Arguments

- data:

  Data frame with one row per individual or per population (species in a
  site), such as from
  [`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md).
  Rows of one population are averaged.

- trait:

  Name of the column with trait values.

- site, species:

  Names of the columns with the site and the species.

- abundance:

  Name of the column with the abundance of the species in the site (one
  value per population), or `NULL` to weight each occurrence equally.

- site_weights, species_weights:

  `"abundance"` or `"occurrence"`: weights of the site means and of the
  species means and \\\beta\\. Both default to `"abundance"` when
  `abundance` is given.

- leave_one_out:

  Leave the focal species out of the site means (the default).

- min_slope:

  Smallest number of sites with measured values to fit the slope
  \\b_i\\.

## Value

An object of class `mc_trait_gradient`, a list with

- species:

  one row per species: `species`, `n_sites`, `n_alone` (sites where the
  species is the only species; left out of its values when
  `leave_one_out = TRUE`), `mean` (\\\bar t_i\\), `beta`, `alpha`,
  `breadth` (\\R_i\\), `slope` (\\b_i\\), `slope_se`, `intercept` and
  `min_site`, `max_site` (range of the site means of its sites).

- sites:

  one row per site: `site`, `n_species` and `mean` (\\p_j\\, with all
  species).

- populations:

  one row per population: `site`, `species`, `abundance` (1 without
  `abundance`), `weight` (the weight of the population in the species
  means), `trait` (with species means filled in), `measured` and
  `site_mean` (\\p_j\\, or \\p_j^{(-i)}\\).

- settings:

  the arguments.

## Details

With \\t\_{ij}\\ the trait value of species \\i\\ in site \\j\\ and
\\a\_{ij}\\ its weight, the mean of site \\j\\ is \$\$p_j = \frac{\sum_i
a\_{ij} t\_{ij}}{\sum_i a\_{ij}},\$\$ the mean of species \\i\\ is
\\\bar t_i = \sum_j a\_{ij} t\_{ij} / \sum_j a\_{ij}\\, its position
along the gradient is \$\$\beta_i = \frac{\sum_j a\_{ij} p_j}{\sum_j
a\_{ij}}\$\$ and its difference from the co-occurring species is
\\\alpha_i = \bar t_i - \beta_i\\. The niche breadth \\R_i\\ is the
range of \\p_j\\ over the sites of the species and the slope \\b_i\\ is
the weighted regression of \\t\_{ij}\\ on \\p_j\\, the response of the
species to the gradient within its range.

## The three changes

**Leave-one-out site means.** In the original analysis \\p_j\\ includes
the focal species, so \\\beta_i\\ contains \\\bar t_i\\: in sites with
few species, \\\beta\\ follows the trait of the species itself, and
\\\beta\\ is correlated with \\\bar t\\ and with \\\alpha\\ even when
traits are unrelated to where species live. With `leave_one_out = TRUE`
(the default) \\\beta_i\\, \\R_i\\ and \\b_i\\ use \\p_j^{(-i)} =
(\sum\_{k \ne i} a\_{kj} t\_{kj}) / \sum\_{k \ne i} a\_{kj}\\, and sites
where species \\i\\ is the only species are left out of its values
(counted in `n_alone`). The jackknife of the original script does the
same, as an option; `leave_one_out = FALSE` gives the original values.

Leaving the species out does not make the raw correlation of \\\beta\\
and \\\alpha\\ readable on its own. Because \\\alpha = \bar t - \beta\\,
without any assembly the correlation is negative, near \\-sd(\beta) /
sd(\alpha)\\, whereas with the focal species in the site means it is
positive. Read it against a null model, as its standardised effect size
in
[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md).

**Weights.** The weights of the site means (`site_weights`) and of the
species means and \\\beta\\ (`species_weights`) are chosen apart, as
abundance or occurrence; the original uses one choice for both. Compare
the two with
[`trait_gradient_weights()`](https://provetelab.org/comdynamite/reference/trait_gradient_weights.md).

**Null model and error.** See
[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md)
and
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md).

Trait values missing in a site where the species occurs (`NA` in
`trait`, with the row kept) are replaced by the mean of the species, as
in the original; the slope \\b_i\\ uses only measured values.

## References

Ackerly, D. D., & Cornwell, W. K. (2007). A trait-based approach to
community assembly: partitioning of species trait values into within-
and among-community components. *Ecology Letters*, 10, 135–145.
[doi:10.1111/j.1461-0248.2006.01006.x](https://doi.org/10.1111/j.1461-0248.2006.01006.x)

## See also

[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md),
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md),
[`trait_gradient_weights()`](https://provetelab.org/comdynamite/reference/trait_gradient_weights.md),
[`plot.mc_trait_gradient()`](https://provetelab.org/comdynamite/reference/plot.mc_trait_gradient.md).

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
[`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md),
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md),
[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md),
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
ind <- sample_individuals(sim, n = 5, local_adaptation = 0.5, intraspecific_sd = 1)
tg <- trait_gradient(ind, abundance = "abundance")
head(tg$species)
#>   species n_sites n_alone      mean     beta       alpha  breadth     slope
#> 1      t1       3       0 67.956541 67.41152   0.5450167 28.45005 0.6848220
#> 2     t10      15       0 51.064276 52.59265  -1.5283787 70.25170 0.6234984
#> 3     t11       4       0 60.927325 61.58517  -0.6578440 21.75283 0.6590566
#> 4     t12       2       0 64.591689 65.46605  -0.8743607 14.06901        NA
#> 5     t13      12       0 29.407255 45.42789 -16.0206324 70.45406 0.3865279
#> 6     t15       1       0  4.898275 22.26285 -17.3645737  0.00000        NA
#>      slope_se intercept min_site max_site
#> 1 0.003962136  21.79165 44.54802 72.99807
#> 2 0.020972694  18.27284 12.66171 82.91341
#> 3 0.050671341  20.33921 51.35616 73.10899
#> 4          NA        NA 51.61021 65.67922
#> 5 0.053463943  11.84811 10.26715 80.72121
#> 6          NA        NA 22.26285 22.26285
# the original, with the focal species in the site means
original <- trait_gradient(ind, abundance = "abundance", leave_one_out = FALSE)
head(cbind(leave_one_out = tg$species$beta, original = original$species$beta))
#>      leave_one_out original
#> [1,]      67.41152 67.45006
#> [2,]      52.59265 52.65582
#> [3,]      61.58517 61.55826
#> [4,]      65.46605 65.38875
#> [5,]      45.42789 37.71938
#> [6,]      22.26285 12.75160
```
