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

  Number of posterior draws to plot.

- ...:

  Further arguments to
  [`bayesplot::ppc_dens_overlay()`](https://mc-stan.org/bayesplot/reference/PPC-distributions.html).

## Value

A ggplot object from bayesplot.
