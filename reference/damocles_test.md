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
  shared = FALSE
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
  communities, one value for all communities or one per community;
  `NULL` (the default) estimates both with
  [`fit_damocles()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md).

- root:

  Prior on the state at the root, used both in the fit and in the
  simulations (`"both"` draws the state of the root of each simulated
  community with equal probabilities).

- shared:

  When the rates are estimated, `TRUE` fits one pair shared by all
  communities (see
  [`fit_damocles()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md))
  and `FALSE` (the default) one pair per community.

## Value

A data frame of class `mc_damocles_test`, one row per community, with
the richness, the rates used (`gamma`, `mu`) and their log-likelihood
when estimated, and, for `mpd` and `mntd`, the observed value, the
standardised effect size and p-value against the DAMOCLES null
(`ses_mpd_damocles`, `p_mpd_damocles`, ...) and against random draws
(`ses_mpd_random`, `p_mpd_random`, ...), plus `n_valid`.

## Details

Allopatric speciation alone makes present species less related than
random draws expect, so communities that look overdispersed against
random draws may be consistent with the dynamic model (Pigot & Etienne
2015). Richness is not fixed in the communities simulated under
DAMOCLES; simulated communities with fewer than two species have no MPD
or MNTD and are left out of the null distribution, and `n_valid` counts
the ones kept. Distances are the cophenetic distances of `tree`, which
should hold every species of the regional pool.

Standardised effect sizes are
`(observed - mean of the null) / sd of the null`, as in
[`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md):
negative values mean the species present are more closely related than
expected (clustering) and positive values that they are less related
(overdispersion). Note that `DAMOCLES::DAMOCLES_bootstrap()` reports the
opposite sign. The p-values are lower-tail ranks,
`(number of null values <= observed + 1) / (n_valid + 1)`; values near
one point to overdispersion.

## References

Pigot, A. L., & Etienne, R. S. (2015). A new dynamic null model for
phylogenetic community structure. *Ecology Letters*, 18, 153–163.
[doi:10.1111/ele.12395](https://doi.org/10.1111/ele.12395)

## See also

[`simulate_damocles()`](https://provetelab.org/comdynamite/reference/simulate_damocles.md),
[`fit_damocles()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md).

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
