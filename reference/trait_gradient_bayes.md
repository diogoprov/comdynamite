# Bayesian phylogenetic model of the trait gradient, with latent beta

Fits in Stan a bivariate phylogenetic mixed model in which the position
of each species along the gradient (\\\beta\\) and its mean trait
(\\\bar t\\) are latent species effects, estimated from every site of
the species, and returns the phylogenetic signal of \\\beta\\,
\\\alpha\\ and \\\bar t\\ with their credible intervals. Unlike
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md),
which first computes \\\beta\\ and then models it with a standard error,
here the uncertainty of \\\beta\\ enters the model directly.

## Usage

``` r
trait_gradient_bayes(
  x,
  tree,
  ou_alpha = NULL,
  lambda_prior = c(1, 1),
  chains = 4,
  parallel_chains = 1,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  seed = NULL,
  ...
)
```

## Arguments

- x:

  Result of
  [`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md).

- tree:

  Phylogeny (`phylo`) with branch lengths; species not in it are left
  out.

- ou_alpha:

  `NULL` for Brownian motion, or the strength \\a\\ of an
  Ornstein-Uhlenbeck process, in units of 1 / branch length.

- lambda_prior:

  The two shape parameters of the beta prior of both \\\lambda\\.

- chains, parallel_chains, iter_warmup, iter_sampling, adapt_delta,
  seed:

  Passed to the `$sample()` method of cmdstanr.

- ...:

  Further arguments to `$sample()`, such as `max_treedepth` or
  `refresh`.

## Value

An object of class `mc_trait_gradient_bayes`, a list with

- signal:

  the phylogenetic share \\\lambda\\ of \\\beta\\, \\\alpha\\ and the
  mean, and the phylogenetic correlations of \\\beta\\ with the mean and
  with \\\alpha\\: `estimate` (posterior median), `lower`, `upper` (95%
  interval).

- variances:

  the phylogenetic and species variances of the three quantities, in
  squared trait units.

- species:

  posterior median and 95% interval of \\\beta\\, \\\alpha\\ and the
  mean of each species.

- diagnostics:

  largest R-hat, smallest bulk and tail effective sample sizes, number
  of divergent transitions, and `ok`.

- fit:

  the `CmdStanMCMC` object; its draws live in temporary files of the
  session, so save them with `fit$save_object()` to keep them.

- data, A, ou_alpha, lambda_prior:

  the data of the model, the phylogenetic correlation matrix and the
  settings.

## Details

Each population (species \\i\\ in site \\j\\) contributes two
observations: the mean of its site, \\p\_{ij}\\ (without the focal
species, as chosen in
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md)),
and its trait value \\t\_{ij}\\. The model is \$\$p\_{ij} = \mu_p +
u^p_i + v^p_i + e^p\_{ij}, \qquad t\_{ij} = \mu_t + u^t_i + v^t_i +
e^t\_{ij},\$\$ where \\(u^p, u^t)\\ are phylogenetic effects with
covariance \\\Sigma_u \otimes A\\ (\\A\\ the phylogenetic correlation
matrix), \\(v^p, v^t)\\ are species effects without phylogeny with
covariance \\\Sigma_v \otimes I\\, and the residuals are independent. So
\\\beta_i = \mu_p + u^p_i + v^p_i\\, \\\bar t_i = \mu_t + u^t_i +
v^t_i\\ and \\\alpha_i = \bar t_i - \beta_i\\, computed for each draw of
the posterior.

The variance among species of each response is written as a total,
\\\tau^2\\, and a phylogenetic share \\\lambda\\: \\\sigma_u = \tau
\sqrt{\lambda}\\ and \\\sigma_v = \tau \sqrt{1 - \lambda}\\. \\\lambda\\
is between 0 (no phylogenetic signal) and 1 (the structure of `A` only:
Brownian motion, or the OU process given), and its prior is set
directly, `beta(lambda_prior[1], lambda_prior[2])`, uniform by default.
(A model with independent priors on \\\sigma_u\\ and \\\sigma_v\\, as in
brms, puts an implicit U-shaped prior on \\\lambda\\.) For \\\alpha\\,
\\\sigma^2_u(\alpha) = \sigma^2_u(t) + \sigma^2_u(p) - 2 \rho_u
\sigma_u(t) \sigma_u(p)\\, and likewise for \\\sigma^2_v(\alpha)\\; its
\\\lambda\\ follows. The phylogenetic correlation of \\\beta\\ with
\\\bar t\\ is \\\rho_u\\, and with \\\alpha\\ it follows from the same
terms. Both effects are non-centred, so the sampler works when either
share is near 0.

With `ou_alpha`, \\A\\ is the correlation matrix of an
Ornstein-Uhlenbeck process with strength \\a\\ = `ou_alpha` and a fixed
root, as in
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md):
\\C\_{ij} = e^{-a(d_i + d_j - 2 s\_{ij})} (1 - e^{-2 a s\_{ij}}) /
(2a)\\, scaled to correlations. \\a\\ is taken as known, not estimated:
take it from the OU fit of
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md)
for \\\beta\\ (`response == "beta"`, `model == "OU"`, column `a`), which
is in the same units. The same \\a\\ applies to \\\beta\\ and \\\bar
t\\, because their phylogenetic correlation needs one matrix. As \\a\\
goes to 0 the matrix is that of Brownian motion; at large \\a\\ it
approaches the identity, and the phylogenetic and species effects can no
longer be told apart (the function warns). On a tree that is not
ultrametric, scaling the Brownian covariance to correlations is not
Brownian motion, because the tips have different variances.

Both responses are centred on the mean trait value and divided by its
standard deviation, so that \\\alpha\\ keeps its meaning, and the
results are returned in trait units. Weights (from `species_weights` in
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md))
multiply the log-likelihood of each population and are divided by their
mean, so that they change the relative influence of populations without
changing their number. A weighted observation then counts as that many
full observations, so unequal weights, such as abundances, make the
posterior too narrow: the effective sample is \\(\sum w)^2 / \sum w^2\\
rather than the number of populations. The function warns when this is
below 70% of the populations; `species_weights = "occurrence"` avoids
the problem.

**Priors**, on the standardised scale: `normal(0, 1)` for the
intercepts, `exponential(1)` for the total standard deviations among
species and for the residuals, `beta(1, 1)` for both \\\lambda\\ (change
with `lambda_prior`), and `lkj(2)` for the two correlations. Check
convergence (`diagnostics`: R-hat below 1.01, bulk and tail effective
sample sizes above 400 and no divergent transitions, over the
intercepts, standard deviations, \\\lambda\\, correlations, species
values and `lp__`) and the fit, with
[`bayesplot::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
on the result, before reading the estimates.

The model is in
`system.file("stan", "trait_gradient.stan", package = "comdynamite")`.
It is compiled by cmdstanr the first time it is used and kept in
`tools::R_user_dir("comdynamite", "cache")`. The pointwise
log-likelihoods are in `log_lik_p` (site means) and `log_lik_t` (trait
values): cross-validation, for instance to compare values of `ou_alpha`,
needs one response at a time, `fit$loo(variables = "log_lik_p")`,
because the two responses are fitted to different sets of populations.

## Limits

The phylogenetic part is Brownian motion or an Ornstein-Uhlenbeck
process with a given strength, shared by both responses;
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md)
estimates the strength. Site means of species that share a site are
almost the same number, but the model treats them as independent
observations, so the posterior of \\\beta\\ is somewhat too narrow when
species share many sites. The site mean of a population also contains
the trait values of the species it lives with, so when relatives
co-occur, the phylogenetic correlation \\\rho_u\\ of \\\beta\\ and
\\\bar t\\ is inflated by their shared values, not only by niche
conservatism. The residuals of the two responses are taken as
independent, which ignores the slope of the trait on the gradient within
species;
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md)
estimates those slopes. Sites where a species is the only species have
no site mean without it and enter only its trait; a species that is
alone in all its sites has no data on \\\beta\\, and its \\\beta\\ and
\\\alpha\\ are returned as `NA`.

Splitting the variance among species into a phylogenetic and a
non-phylogenetic part needs many species, and the intervals of
\\\lambda\\ stay wide. In 60 trials with 40–45 species, 30 sites,
occurrence weights and the default settings, where the niches evolved
with Pagel's \\\lambda\\ of 0, 0.5 or 1 (20 trials each):

|  |  |  |  |
|----|----|----|----|
|  | \\\lambda = 0\\ | \\\lambda = 0.5\\ | \\\lambda = 1\\ |
| median of \\\lambda\_\beta\\ (range of the 20 medians) | 0.07 (0.05–0.31) | 0.24 (0.04–0.73) | 0.94 (0.38–0.99) |
| 95\\ lower bound above 0.05 / upper bound below 0.95 | 0 / 20 | 5 / 20 | 16 / 1 |
| passed every diagnostic | 12 | 11 | 8 |

A fit took 9–30 s with four chains in parallel. One of the 60 fits had a
divergent transition; the fits that failed the diagnostics did so by an
effective sample size below 400 or an R-hat between 1.01 and 1.03, for
which longer chains (`iter_sampling`) are the first remedy (not tested
here). Intermediate signal is estimated low (interpretation: the
\\\lambda\\ of \\\beta\\ is below that of the niches, because \\\beta\\
also carries the error of the site means), and its interval rarely
excludes 0: a wide interval for \\\lambda\_\beta\\ means that the data
do not tell. An earlier version of this model fitted with brms, with
independent priors on the two standard deviations, had divergent
transitions in most fits and pulled intermediate signal to 0 (median
0.04, interval containing 0.5 in 12 of the same 20 trials). In trials
with 18 species the chains of that version did not mix. Use
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md)
to test the signal as well, and do not read a fit that fails
`diagnostics`.

## See also

[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md),
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md)

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
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md),
[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md),
[`trait_gradient_weights()`](https://provetelab.org/comdynamite/reference/trait_gradient_weights.md),
[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
[`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md)

## Examples

``` r
if (FALSE) { # interactive() && requireNamespace("cmdstanr", quietly = TRUE) && requireNamespace("treeshapesim", quietly = TRUE)
set.seed(1)
tree <- treeshapesim::sim_tree_shape(40)
sim <- simulate_metacommunity(landscape_implicit(30), species_traits(tree, lambda = 1),
                              competition_matrix(tree, "none"), n_steps = 100,
                              seed_steps = 1:5)
ind <- sample_individuals(sim, n = 5, local_adaptation = 0.5, intraspecific_sd = 1)
tg <- trait_gradient(ind, abundance = "abundance", species_weights = "occurrence")
bayes <- trait_gradient_bayes(tg, tree, parallel_chains = 4)
bayes$diagnostics
bayes$signal
bayesplot::pp_check(bayes, resp = "sitemean")
}
```
