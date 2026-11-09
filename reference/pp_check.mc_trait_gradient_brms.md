# Posterior predictive check of the trait-gradient model

Posterior predictive check of the trait-gradient model

## Usage

``` r
# S3 method for class 'mc_trait_gradient_brms'
pp_check(object, resp = c("sitemean", "trait"), ...)
```

## Arguments

- object:

  Result of
  [`trait_gradient_brms()`](https://provetelab.org/comdynamite/reference/trait_gradient_brms.md).

- resp:

  `"sitemean"` or `"trait"`.

- ...:

  Further arguments to
  [`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html).

## Value

A ggplot object from bayesplot.
