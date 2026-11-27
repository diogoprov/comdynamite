# Statistics of trait variation within and among populations

Summarises traits measured on individuals of several populations of each
species, such as the output of
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md)
or field data, into statistics that carry information on the genetic
variance, on drift and on gene flow, for inference by approximate
Bayesian computation
([`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md))
next to
[`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md).
Each statistic that lacks its inputs is `NA`, so the result always has
the same names.

## Usage

``` r
population_stats(
  individuals,
  trait = "trait",
  site = "site",
  species = "species",
  env = "env",
  abundance = NULL,
  step = NULL,
  coords = NULL,
  near = 1.5,
  far = 3,
  by_species = FALSE
)
```

## Arguments

- individuals:

  Data frame with one row per individual.

- trait, site, species, env:

  Names of the columns with the trait (numeric), the site, the species
  and the environment of the site (numeric), as in the output of
  [`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md).

- abundance:

  Name of a column with the abundance of the population of each
  individual, one value per population, or `NULL` (the default), which
  leaves `within_slope`, `size_cor` and `temporal_slope` as `NA`.

- step:

  Name of a column with the survey of each individual, or `NULL` (the
  default) for a single survey, which leaves `temporal`,
  `temporal_slope`, `temporal_step` and `temporal_step_slope` as `NA`.
  Without `step`, all the individuals of a population are pooled.

- coords:

  Matrix or data frame with the coordinates of the sites in two columns
  and the names of the sites as row names, or `NULL` (the default),
  which leaves `near`, `far` and `dist_cor` as `NA`.

- near, far:

  Distances (non-negative numbers, `near` no larger than `far`) that
  define neighbouring and distant pairs of populations.

- by_species:

  If `TRUE`, a data frame with the statistics of each species and its
  weight instead of the averages.

## Value

Named numeric vector, or with `by_species = TRUE` a data frame with one
row per species.

## Details

A population is a species in a site. The statistics of a species use its
populations with at least two individuals measured at the last survey,
and need three such populations at sites that differ in environment.
`slope` and `r2` are the slope and \\R^2\\ of a linear regression of
population means on the environment, `sd` the standard deviation of the
means and `resid_sd` that of the residuals. `within` is the logarithm of
the mean variance within populations (the phenotypic variance, genetic
plus environmental) and `within_slope` the slope of the variance within
populations on the logarithm of abundance, positive when drift erodes
the genetic variance of small populations. `qst` is the residual
variance among population means, \\V_B\\, relative to \\V_B + 2 V_W\\,
where \\V_W\\ is the mean variance within populations, a divergence
among populations beyond what the environment explains, in the form of
\\Q\_{ST}\\ (Spitze 1993) but with phenotypes. `near` and `far` are the
semivariances of residuals between populations closer than `near` and
farther than `far`, relative to the residual variance, and `dist_cor`
the correlation of squared differences of residuals with distance: drift
makes neighbours differ, gene flow makes them alike. `size_cor` is the
correlation of absolute residuals with the logarithm of abundance,
negative when small populations drift more.

With two or more surveys (`step`), the first and the last are compared,
in the populations with at least two individuals measured in both.
`temporal` is the mean squared change of population means minus its
sampling variance, relative to the variance within populations, and
`temporal_slope` the slope of that corrected change on the inverse of
abundance, the basis of the temporal method of estimating the effective
size (Waples 1989). Over a long interval, selection and gene flow also
move the means, so `temporal_step` uses every pair of consecutive
surveys instead: the mean, over these pairs and their populations
measured in both, of the squared change minus its sampling variance,
divided by the interval between the surveys (in the units of `step`, so
per generation when a step is a generation) and by the variance within
populations at the last survey. It estimates the variance of drift per
generation relative to the phenotypic variance, which is about \\G /
(N_e (G + V_E))\\ (Lande 1976). `temporal_step_slope` is the slope of
the corrected change per unit of time on the inverse of the abundance at
the later survey. Both need a numeric `step` and three changes in all;
with two surveys, `temporal_step` is `temporal` divided by the interval.

The statistics of each species are averaged over species, weighted by
the summed abundance of their populations at the last survey (or by
their number of populations without `abundance`). `slope_min` and
`slope_max` are the smallest and largest slope among species. Means and
variances are those of the measured individuals, so they include
sampling error, as field data would;
[`subsample_individuals()`](https://provetelab.org/comdynamite/reference/subsample_individuals.md)
changes the effort.

## References

Lande, R. (1976). Natural selection and random genetic drift in
phenotypic evolution. *Evolution*, 30, 314–334.
[doi:10.1111/j.1558-5646.1976.tb00911.x](https://doi.org/10.1111/j.1558-5646.1976.tb00911.x)

Spitze, K. (1993). Population structure in *Daphnia obtusa*:
quantitative genetic and allozymic variation. *Genetics*, 135, 367–374.
[doi:10.1093/genetics/135.2.367](https://doi.org/10.1093/genetics/135.2.367)

Waples, R. S. (1989). A generalized approach for estimating effective
population size from temporal changes in allele frequency. *Genetics*,
121, 379–391.
[doi:10.1093/genetics/121.2.379](https://doi.org/10.1093/genetics/121.2.379)

## See also

Other trait variation among individuals:
[`coef.drift_test()`](https://provetelab.org/comdynamite/reference/coef.drift_test.md),
[`coef.mc_trait_gradient()`](https://provetelab.org/comdynamite/reference/coef.mc_trait_gradient.md),
[`coef.population_shifts()`](https://provetelab.org/comdynamite/reference/coef.population_shifts.md),
[`community_index()`](https://provetelab.org/comdynamite/reference/community_index.md),
[`drift_test()`](https://provetelab.org/comdynamite/reference/drift_test.md),
[`neighbour_metrics`](https://provetelab.org/comdynamite/reference/neighbour_metrics.md),
[`nested_variance()`](https://provetelab.org/comdynamite/reference/nested_variance.md),
[`null_intervals()`](https://provetelab.org/comdynamite/reference/null_intervals.md),
[`plot.mc_trait_gradient()`](https://provetelab.org/comdynamite/reference/plot.mc_trait_gradient.md),
[`plot.mc_trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/plot.mc_trait_gradient_bayes.md),
[`plot.mc_trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/plot.mc_trait_gradient_phylo.md),
[`plot.nested_variance()`](https://provetelab.org/comdynamite/reference/plot.nested_variance.md),
[`plot.trait_flex()`](https://provetelab.org/comdynamite/reference/plot.trait_flex.md),
[`plot.trait_rarefaction()`](https://provetelab.org/comdynamite/reference/plot.trait_rarefaction.md),
[`plot.tstats()`](https://provetelab.org/comdynamite/reference/plot.tstats.md),
[`plot_population_traits()`](https://provetelab.org/comdynamite/reference/plot_population_traits.md),
[`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md),
[`print.trait_flex()`](https://provetelab.org/comdynamite/reference/print.trait_flex.md),
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md),
[`subsample_individuals()`](https://provetelab.org/comdynamite/reference/subsample_individuals.md),
[`trait_change()`](https://provetelab.org/comdynamite/reference/trait_change.md),
[`trait_densities()`](https://provetelab.org/comdynamite/reference/trait_densities.md),
[`trait_flex()`](https://provetelab.org/comdynamite/reference/trait_flex.md),
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md),
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
land <- landscape_grid(side = 6, env_bandwidth = 2)
traits <- data.frame(optimum = c(30, 50, 70), sigma = 15, r_max = 5,
                     row.names = c("a", "b", "c"))
alpha <- matrix(0.01, 3, 3, dimnames = list(rownames(traits), rownames(traits)))
diag(alpha) <- 0.02
sim <- simulate_metacommunity(land, traits, alpha, dispersal = 0.01, n_steps = 120,
                              record_steps = c(100, 120),
                              evolution = evolution_control(4, ne_ratio = 0.2, env_var = 1))
ind <- sample_individuals(sim, n = 10, source = "evolution", step = c(100, 120))
coords <- land$coords
rownames(coords) <- rownames(sim$abundance)
round(population_stats(ind, abundance = "abundance", step = "step", coords = coords), 3)
#>               slope                  r2                  sd            resid_sd 
#>               0.204               0.218               8.577               7.463 
#>              within        within_slope                 qst                near 
#>               1.596               0.166               0.846               0.976 
#>                 far            dist_cor            size_cor            temporal 
#>               0.920              -0.042              -0.223               0.414 
#>      temporal_slope       temporal_step temporal_step_slope           slope_min 
#>              29.080               0.021               1.454               0.107 
#>           slope_max 
#>               0.350 
```
