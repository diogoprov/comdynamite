# Plot the Bayesian model of the trait gradient

`type = "intervals"` draws the posterior median and central 95% interval
of each quantity in `signal` of
[`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md):
the phylogenetic share \\\lambda\\ of \\\beta\\, \\\alpha\\ and the
mean, and the phylogenetic correlations of \\\beta\\ with the mean and
with \\\alpha\\. With `species = TRUE` it draws instead the medians and
intervals of \\\beta\\, \\\alpha\\ and the mean of each species.
`type = "trace"` draws the chains of `variables` with
[`bayesplot::mcmc_trace()`](https://mc-stan.org/bayesplot/reference/MCMC-traces.html).
Read neither before checking `diagnostics`.

## Usage

``` r
# S3 method for class 'mc_trait_gradient_bayes'
plot(
  x,
  type = c("intervals", "trace"),
  species = FALSE,
  variables = c("lambda_beta", "lambda_alpha", "lambda_mean"),
  inc_warmup = FALSE,
  colour = "#0072B2",
  ...
)
```

## Arguments

- x:

  Result of
  [`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md).

- type:

  `"intervals"` or `"trace"`.

- species:

  `TRUE` or `FALSE`. With `type = "intervals"`, draw the values of the
  species instead of `signal`.

- variables:

  Names of the variables of the Stan model for `type = "trace"`, such as
  `"lambda_beta"`, `"rho_u"` or `"beta[1]"`.

- inc_warmup:

  `TRUE` or `FALSE`. Include the warmup draws in `type = "trace"`; see
  Details.

- colour:

  Colour of the points and intervals.

- ...:

  Further arguments to
  [`bayesplot::mcmc_trace()`](https://mc-stan.org/bayesplot/reference/MCMC-traces.html)
  for `type = "trace"`, and not used otherwise.

## Value

A ggplot object.

## Details

Warmup draws are shown with `inc_warmup = TRUE` only when the model was
fitted with `save_warmup = TRUE`, passed through the `...` of
[`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md)
to cmdstanr; otherwise the fit has no warmup draws and cmdstanr stops
with an error.

## See also

[`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md),
[`pp_check.mc_trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/pp_check.mc_trait_gradient_bayes.md)

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
[`plot.mc_trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/plot.mc_trait_gradient_phylo.md),
[`plot.nested_variance()`](https://provetelab.org/comdynamite/reference/plot.nested_variance.md),
[`plot.trait_flex()`](https://provetelab.org/comdynamite/reference/plot.trait_flex.md),
[`plot.trait_rarefaction()`](https://provetelab.org/comdynamite/reference/plot.trait_rarefaction.md),
[`plot.tstats()`](https://provetelab.org/comdynamite/reference/plot.tstats.md),
[`plot_population_traits()`](https://provetelab.org/comdynamite/reference/plot_population_traits.md),
[`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md),
[`population_stats()`](https://provetelab.org/comdynamite/reference/population_stats.md),
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
if (FALSE) { # interactive() && requireNamespace("cmdstanr", quietly = TRUE) && requireNamespace("bayesplot", quietly = TRUE) && requireNamespace("treeshapesim", quietly = TRUE)
set.seed(1)
tree <- treeshapesim::sim_tree_shape(30)
sim <- simulate_metacommunity(landscape_implicit(20), species_traits(tree, lambda = 1),
                              competition_matrix(tree, "none"), n_steps = 100,
                              seed_steps = 1:5)
ind <- sample_individuals(sim, n = 5, local_adaptation = 0.5, intraspecific_sd = 1)
tg <- trait_gradient(ind, abundance = "abundance", species_weights = "occurrence")
bayes <- trait_gradient_bayes(tg, tree, parallel_chains = 4, save_warmup = TRUE)
plot(bayes)
plot(bayes, species = TRUE)
plot(bayes, type = "trace", inc_warmup = TRUE)
}
```
