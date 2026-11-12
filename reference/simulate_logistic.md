# Simulate presences and absences from logistic responses

The two-gradient simulation of Appendix C of Peres-Neto et al. (2012),
after Ives & Helmus (2011). The probability that species \\j\\ occurs at
site \\t\\ along gradient \\x\\ is \\\mathrm{logit}^{-1}(a_j + b_j
x_t)\\, where \\a\\ is a normal deviate and the slope \\b\\ evolves on
the tree. With several gradients, each with its own \\a\\ and \\b\\, the
probabilities multiply. A value drawn for each site from `reduction`
then reduces them, which adds noise.

## Usage

``` r
simulate_logistic(
  tree,
  gradients,
  intercept_mean = 0,
  intercept_sd = 0.5,
  slope_sd = 1,
  lambda = 1,
  reduction = c(0.5, 1),
  draw = c("uniform", "normal")
)
```

## Arguments

- tree:

  A `phylo` tree with branch lengths. Appendix C uses a fully balanced
  tree with 32 species.

- gradients:

  Value of the gradient at each site, or a list or data frame with one
  gradient per element, such as
  [`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md).
  Gradients are centred.

- intercept_mean, intercept_sd:

  Mean and standard deviation of `a`. Appendix C does not give these
  values, and the defaults are those of the code of Ives & Helmus
  (2011).

- slope_sd:

  Standard deviation of the slopes `b` among the species of the tree, as
  in the code of Ives & Helmus (2011). Appendix C does not give this
  value. The slopes are rescaled to this sample standard deviation in
  each simulation, so with `lambda > 0` the variance of the Brownian
  process behind them is larger than `slope_sd^2`; a model that
  estimates that variance, such as a phylogenetic GLMM, should be
  compared with it and not with `slope_sd^2` (see
  `data-raw/validation_pglmm.R`).

- lambda:

  Pagel's lambda of the slopes.

- reduction:

  Range of the uniform distribution of the site reduction.

- draw:

  Distribution of the value compared with each probability.

## Value

List of class `mc_sim` with `abundance`, a 0/1 matrix with one row per
site and one column per species; `probability`, the matrix of
probabilities; `env`, the gradients; and `traits`, a list with the
matrices `intercept` and `slope`, with one row per species and one
column per gradient; and `tree`.

## Details

Appendix C says that a species is present when a normally generated
value is smaller than its probability. Taken literally, a species with a
probability near zero is present half of the time, so the default is a
uniform value (a Bernoulli draw), and `draw = "normal"` follows the
text.

## References

Ives, A. R., & Helmus, M. R. (2011). Generalized linear mixed models for
phylogenetic analyses of community structure. *Ecological Monographs*,
81, 511–525. [doi:10.1890/10-1264.1](https://doi.org/10.1890/10-1264.1)

Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the
effects of spatial contingency and environmental filtering on
metacommunity phylogenetics. *Ecology*, 93, S14–S30.
[doi:10.1890/11-0494.1](https://doi.org/10.1890/11-0494.1)

## See also

Other simulate metacommunities:
[`colonisation_history()`](https://provetelab.org/comdynamite/reference/colonisation_history.md),
[`evolution_control()`](https://provetelab.org/comdynamite/reference/evolution_control.md),
[`growth_rates()`](https://provetelab.org/comdynamite/reference/growth_rates.md),
[`logistic_gradients()`](https://provetelab.org/comdynamite/reference/logistic_gradients.md),
[`mc_dynamics()`](https://provetelab.org/comdynamite/reference/mc_dynamics.md),
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md),
[`simulate_static()`](https://provetelab.org/comdynamite/reference/simulate_static.md)

## Examples

``` r
set.seed(1)
tr <- ape::compute.brlen(ape::stree(32, "balanced"), method = "Grafen")
simulate_logistic(tr, logistic_gradients(1)[c("x", "y")])
#> <mc_sim> 31 sites x 32 species
#> species persisting: 32 of 32 
#> mean local richness: 5.65 
```
