# Posterior predictive check of the trait-gradient model

Overlays the density of the observed site means or trait values on the
densities of replicated data drawn from the posterior of
[`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md).

## Usage

``` r
# S3 method for class 'mc_trait_gradient_bayes'
pp_check(object, resp = c("sitemean", "trait"), ndraws = 50, ...)
```

## Arguments

- object:

  Result of
  [`trait_gradient_bayes()`](https://provetelab.org/comdynamite/reference/trait_gradient_bayes.md).

- resp:

  `"sitemean"` or `"trait"`.

- ndraws:

  Number of posterior draws to plot. One whole number, 1 or more.

- ...:

  Further arguments to
  [`bayesplot::ppc_dens_overlay()`](https://mc-stan.org/bayesplot/reference/PPC-distributions.html).

## Value

A ggplot object from bayesplot.

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
bayes <- trait_gradient_bayes(tg, tree, parallel_chains = 4)
bayesplot::pp_check(bayes, resp = "trait", ndraws = 20)
}
```
