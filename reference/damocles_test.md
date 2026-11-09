# Phylogenetic community structure against the dynamic null model of DAMOCLES

`damocles_test()` compares the mean pairwise distance (MPD) and the mean
nearest taxon distance (MNTD) among the species present in each
community with two null distributions: one from communities simulated
along the tree under the DAMOCLES model (Pigot & Etienne 2015;
[`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md)),
with the rates of colonisation and local extinction estimated from that
community (or shared by all communities) by
[`fit_damocles()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md)
or given, and one from random draws of the same number of species from
the tips, the static null model the dynamic one is compared with in that
study.

## Usage

``` r
damocles_test(
  comm,
  tree,
  n_sim = 999,
  gamma = NULL,
  mu = NULL,
  root = c("both", "absent", "present"),
  shared = FALSE,
  sympatry = 0,
  regimes = NULL,
  edge_covariate = NULL,
  vary = c("both", "gamma", "mu")
)
```

## Arguments

- comm:

  Presences of the species (columns, named as the tips of `tree`) in one
  or more communities (rows): a matrix or data frame, or a named vector
  for one community. Values above zero count as presences. Tips missing
  from `comm` count as absent, with a warning, and species of `comm`
  missing from the tree are dropped.

- tree:

  A rooted, fully dichotomous `phylo` tree with branch lengths.

- n_sim:

  Number of null communities of each model.

- gamma, mu:

  Rates of local colonisation and extinction used to simulate the null
  communities, one value for all communities or one per community, or
  lists with one rate, or one rate per edge of `tree`, per community;
  `NULL` (the default) estimates both with
  [`fit_damocles()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md),
  with `shared`, `sympatry`, `regimes`, `edge_covariate` and `vary`
  passed on to it.

- root:

  Prior on the state at the root, used both in the fit and in the
  simulations (`"both"` draws the state of the root of each simulated
  community with equal probabilities).

- shared:

  When the rates are estimated, `TRUE` fits one pair shared by all
  communities (see
  [`fit_damocles()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md))
  and `FALSE` (the default) one pair per community.

- sympatry:

  Probability that a present ancestor leaves both daughters present,
  used in the simulations (see
  [`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md));
  `NULL` estimates it with the rates.

- regimes:

  Optional regime of each edge of `tree`, in the order of `tree$edge` (a
  factor or character vector), such as from
  [`clade_regimes()`](https://provetelab.org/comdynamite/reference/clade_regimes.md)
  or from a configuration of rate shifts estimated by another method:
  the rates in `vary` are estimated separately for each regime.

- edge_covariate:

  Optional positive number for each edge of `tree`, in the order of
  `tree$edge`, such as branch-specific rates of speciation estimated by
  another method: the rates in `vary` change with it as \\\mathrm{rate}
  \cdot (x / \bar{x}\_g)^\beta\\, where \\x\\ is the covariate and
  \\\bar{x}\_g\\ its geometric mean, with `rate` and \\\beta\\
  estimated. Give `regimes` or `edge_covariate`, not both.

- vary:

  Which rates change among `regimes` or with `edge_covariate`: `"both"`,
  `"gamma"` or `"mu"`; the other one is a single rate for the whole
  tree.

## Value

A data frame of class `mc_damocles_test`, one row per community, with
the richness, the parameters used (`gamma`, `mu`, or the columns of the
fit with `regimes` or `edge_covariate`, and `sympatry`; `NA` for rates
given per edge) and their log-likelihood when estimated, and, for `mpd`
and `mntd`, the observed value, the standardised effect size and p-value
against the DAMOCLES null (`ses_mpd_damocles`, `p_mpd_damocles`, ...)
and against random draws (`ses_mpd_random`, `p_mpd_random`, ...), plus
`n_valid`.

## Details

Allopatric speciation alone makes present species less related than
random draws expect, so communities that look overdispersed against
random draws may be consistent with the dynamic model (Pigot & Etienne
2015). Richness is not fixed in the communities simulated under
DAMOCLES; simulated communities with fewer than two species have no MPD
or MNTD and are left out of the null distribution, and `n_valid` counts
the ones kept. Distances are the cophenetic distances of `tree`, which
should hold every species of the regional pool.

Standardised effect sizes are \\(x\_{obs} - \bar{x}\_{null}) /
s\_{null}\\, as in
[`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md):
negative values mean the species present are more closely related than
expected (clustering) and positive values that they are less related
(overdispersion). Note that `DAMOCLES::DAMOCLES_bootstrap()` reports the
opposite sign. The p-values are lower-tail ranks, \\(r + 1) / (n + 1)\\,
with \\r\\ the number of null values at or below the observed one and
\\n\\ the number of valid null values; values near one point to
overdispersion.

## References

Pigot, A. L., & Etienne, R. S. (2015). A new dynamic null model for
phylogenetic community structure. *Ecology Letters*, 18, 153–163.
[doi:10.1111/ele.12395](https://doi.org/10.1111/ele.12395)

## See also

[`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md),
[`fit_damocles()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md).

Other dynamic null model:
[`clade_regimes()`](https://provetelab.org/comdynamite/reference/clade_regimes.md),
[`confint.mc_damocles_fit()`](https://provetelab.org/comdynamite/reference/confint.mc_damocles_fit.md),
[`damocles_loglik()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md),
[`damocles_lrt()`](https://provetelab.org/comdynamite/reference/damocles_lrt.md),
[`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md)

## Examples

``` r
set.seed(3)
tr <- treeshapesim::sim_tree_shape(40)
comm <- simulate_damocles(tr, gamma = 1, mu = 1, n = 2)
damocles_test(comm, tr, n_sim = 99)
#> Phylogenetic structure against the DAMOCLES null model (Pigot & Etienne 2015)
#> 99 null communities per model; SES < 0: clustering, SES > 0: overdispersion
#>      richness gamma    mu  mpd ses_mpd_damocles p_mpd_damocles ses_mpd_random
#> sim1        5  5.85 37.27 1.64          0.54759          0.810          0.730
#> sim2        7  2.03  5.65 1.54         -0.00356          0.414          0.264
#>      p_mpd_random  mntd ses_mntd_damocles p_mntd_damocles ses_mntd_random
#> sim1         0.78 0.878           -0.0181           0.550          -0.171
#> sim2         0.53 0.652           -0.7662           0.222          -0.458
#>      p_mntd_random
#> sim1          0.45
#> sim2          0.36
```
