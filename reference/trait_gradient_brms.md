# Bayesian phylogenetic model of the trait gradient, with latent beta

Fits with brms a bivariate phylogenetic mixed model in which the
position of each species along the gradient (\\\beta\\) and its mean
trait (\\\bar t\\) are latent species effects, estimated from every site
of the species, and returns the phylogenetic signal of \\\beta\\,
\\\alpha\\ and \\\bar t\\ with their credible intervals. Unlike
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md),
which first computes \\\beta\\ and then models it with a standard error,
here the uncertainty of \\\beta\\ enters the model directly.

## Usage

``` r
trait_gradient_brms(
  x,
  tree,
  ou_alpha = NULL,
  prior = NULL,
  chains = 4,
  iter = 2000,
  warmup = floor(iter/2),
  cores = 1,
  seed = NA,
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

- prior:

  Priors as a `brmsprior` object (from
  [`brms::prior()`](https://paulbuerkner.com/brms/reference/set_prior.html)),
  on the standardised scale, or `NULL` for the defaults in Details.

- chains, iter, warmup, cores, seed:

  Passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

- ...:

  Further arguments to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html),
  such as `backend = "cmdstanr"`; `control` defaults to
  `list(adapt_delta = 0.95)`.

## Value

An object of class `mc_trait_gradient_brms`, a list with

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

  the `brmsfit`.

- data, center, scale:

  the data of the model and the standardisation.

## Details

Each population (species \\i\\ in site \\j\\) contributes two
observations: the mean of its site, \\p\_{ij}\\ (without the focal
species, as chosen in
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md)),
and its trait value \\t\_{ij}\\. The model is \$\$p\_{ij} = \mu_p +
u^p_i + v^p_i + e^p\_{ij}, \qquad t\_{ij} = \mu_t + u^t_i + v^t_i +
e^t\_{ij},\$\$ where \\(u^p, u^t)\\ are phylogenetic effects with
covariance \\\Sigma_u \otimes A\\ (\\A\\ the phylogenetic correlation
matrix under Brownian motion), \\(v^p, v^t)\\ are species effects
without phylogeny with covariance \\\Sigma_v \otimes I\\, and the
residuals are independent. So \\\beta_i = \mu_p + u^p_i + v^p_i\\,
\\\bar t_i = \mu_t + u^t_i + v^t_i\\ and \\\alpha_i = \bar t_i -
\beta_i\\, computed for each draw of the posterior.

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
longer be told apart (the function warns). Comparing values of \\a\\ by
cross-validation needs one response at a time,
`brms::loo(fit, resp = "sitemean")`, because the two responses are
fitted to different subsets of rows; the comparison is then conditional
on the other response.

The phylogenetic share of the variance among species is \\\lambda =
\sigma^2_u / (\sigma^2_u + \sigma^2_v)\\, between 0 (no phylogenetic
signal) and 1 (the structure of `A` only: Brownian motion, or the OU
process given). On a tree that is not ultrametric, scaling the Brownian
covariance to correlations is not Brownian motion, because the tips have
different variances. For \\\alpha\\, \\\sigma^2_u(\alpha) =
\sigma^2_u(t) + \sigma^2_u(p) - 2 \rho_u \sigma_u(t) \sigma_u(p)\\, and
likewise for \\\sigma^2_v(\alpha)\\. The phylogenetic correlation of
\\\beta\\ with \\\bar t\\ is \\\rho_u\\, and with \\\alpha\\ it follows
from the same terms.

Both responses are centred on the mean trait value and divided by its
standard deviation, so that \\\alpha\\ keeps its meaning, and the
results are returned in trait units. Weights (from `species_weights` in
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md))
are divided by their mean, so that they change the relative influence of
populations without changing their number. brms treats a weighted
observation as that many full observations, so unequal weights, such as
abundances, make the posterior too narrow: the effective sample is
\\(\sum w)^2 / \sum w^2\\ rather than the number of populations. The
function warns when this is below 70% of the populations;
`species_weights = "occurrence"` avoids the problem.

**Priors**, on the standardised scale: `normal(0, 1)` for the
intercepts, `exponential(1)` for the standard deviations of species
effects and residuals, and `lkj(2)` for the correlations; give `prior`
to change them. Two independent `exponential(1)` priors on the
phylogenetic and species standard deviations put a U-shaped prior on
\\\lambda\\: a quarter of its mass below 0.1 and a quarter above 0.9.
With few species the posterior of \\\lambda\\ follows that prior. Check
convergence (`diagnostics`: R-hat below 1.01, bulk and tail effective
sample sizes above 400, no divergent transitions, over the
population-level, group-level and species effects, `lp__` and the two
\\\lambda\\) and the fit, with
[`bayesplot::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
on the result, before reading the estimates.

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
independent (no residual correlation), which ignores the slope of the
trait on the gradient within species;
[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md)
estimates those slopes. Sites where a species is the only species have
no site mean without it and enter only its trait; a species that is
alone in all its sites has no data on \\\beta\\, and its \\\beta\\ and
\\\alpha\\ are returned as `NA`.

Splitting the variance among species into a phylogenetic and a
non-phylogenetic part needs many species. In trials with 18 species the
posterior of \\\lambda\\ had modes near 0 and 1 and the chains did not
mix (R-hat 1.04–1.27), even with `adapt_delta = 0.99` and tighter
priors. In 20 trials with 39–45 species, 30 sites and the default
settings (4 chains of 2000 iterations), niches that evolved by Brownian
motion gave \\\lambda\_\beta\\ near 1 (median 0.998–0.9995) in 8 of 10
trials, 5 of them passing every diagnostic; niches without phylogenetic
signal gave \\\lambda\_\beta\\ near 0 (median 0.001–0.005) in 8 of 10
trials, but none passed every diagnostic (bulk ESS 180–400, and up to 4
divergent transitions). In the remaining 4 trials the chains split
between the two ends (R-hat 1.39–1.75). With 4000 iterations and
occurrence weights for species, the chains mixed (R-hat at most 1.03)
but most fits had divergent transitions: all 10 without phylogenetic
signal (2–298 of 8000 draws), where \\\lambda\_\beta\\ was again near 0
(median 0.008–0.04). **Intermediate signal is not recovered.** In 20
trials where the niches followed Pagel's \\\lambda = 0.5\\, the
posterior median of \\\lambda\_\beta\\ ranged from 0.007 to 0.81 (median
across trials 0.04): the 95% interval contained 0.5 in 12 and excluded 0
in only 3, and only 4 fits passed every diagnostic. With about 40
species the model tells strong signal from none, but it pulls
intermediate signal toward 0, and a \\\lambda\_\beta\\ near 0 is no
evidence that \\\beta\\ lacks phylogenetic structure. Use
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md)
to test the signal, and do not read a fit that fails `diagnostics`.

## See also

[`trait_gradient()`](https://provetelab.org/comdynamite/reference/trait_gradient.md),
[`trait_gradient_phylo()`](https://provetelab.org/comdynamite/reference/trait_gradient_phylo.md)

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
[`trait_gradient_test()`](https://provetelab.org/comdynamite/reference/trait_gradient_test.md),
[`trait_gradient_weights()`](https://provetelab.org/comdynamite/reference/trait_gradient_weights.md),
[`trait_rarefaction()`](https://provetelab.org/comdynamite/reference/trait_rarefaction.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md),
[`tstats_summary()`](https://provetelab.org/comdynamite/reference/tstats_summary.md)

## Examples

``` r
if (FALSE) { # interactive() && requireNamespace("brms", quietly = TRUE) && requireNamespace("treeshapesim", quietly = TRUE)
set.seed(1)
tree <- treeshapesim::sim_tree_shape(20)
sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree, lambda = 1),
                              competition_matrix(tree, "none"), n_steps = 100,
                              seed_steps = 1:5)
ind <- sample_individuals(sim, n = 5, local_adaptation = 0.5, intraspecific_sd = 1)
tg <- trait_gradient(ind, abundance = "abundance")
bayes <- trait_gradient_brms(tg, tree, chains = 4, cores = 4)
bayes$diagnostics
bayes$signal
bayesplot::pp_check(bayes, resp = "sitemean")
}
```
