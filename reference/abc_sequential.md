# Sequential ABC with random forests

Narrows a wide prior towards the region that fits one observation, by
rounds of simulation and ABC with random forests: the ABC sequential
Monte Carlo with random forests (ABC-SMC-RF) of Dinh et al. (2025). The
first round is a reference table drawn from the prior, as in
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
and one forest per parameter as in
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md).
Each later round draws new parameters from the weighted parameters of
the round before, perturbs them with a normal kernel of twice their
weighted variance (Beaumont et al. 2009), draws again those outside the
prior, simulates and fits new forests on the new simulations only. The
variance of the kernel is that of the parameters of the round before
weighted by their weights. The weights of a round are those the forests
give to its simulations for the observation, times the prior and divided
by the density of the proposal they were drawn from, so that the last
round estimates the posterior under the prior and not under the
proposal.

## Usage

``` r
abc_sequential(
  prior,
  simulate,
  observed,
  n,
  rounds = 4,
  setup = NULL,
  summarise = function(sim, shared) summary_stats(sim, shared$tree),
  table = NULL,
  ntree = 500,
  min_node_size = 5,
  ncores = 1,
  seed = 1
)

# S3 method for class 'mc_abc_sequential'
predict(object, level = 0.9, round = length(object$rounds), ...)

# S3 method for class 'mc_abc_sequential'
plot(x, bins = 30, ...)
```

## Arguments

- prior:

  Named list with one numeric vector of length 2 per parameter, the
  lower and upper limits of its uniform prior, as in
  [`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md).

- simulate, setup, summarise:

  Functions to simulate and summarise, as in
  [`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md).
  `setup` is called once, with the seed that
  [`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md)
  would use, and what it returns is shared by all rounds.

- observed:

  Named numeric vector of observed statistics, as from
  [`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md),
  or a data frame with one row.

- n:

  Number of simulations per round. One whole number, 10 or more.

- rounds:

  Number of rounds, the first from the prior. One whole number, 1 or
  more.

- table:

  Optional reference table from the prior, as from
  [`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md)
  with the same `prior`, `simulate`, `setup` and `summarise`, used as
  the first round instead of new simulations. Its rows need not number
  `n`.

- ntree, min_node_size, ncores:

  Passed to
  [`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md).

- seed:

  Integer that sets all the seeds. One finite whole number.

- object:

  An `mc_abc_sequential` from `abc_sequential()`.

- level:

  Probability of the credible interval. One number greater than 0 and
  less than 1.

- round:

  Round whose weighted parameters give the estimates. By default the
  last.

- ...:

  Not used.

- x:

  An `mc_abc_sequential` from `abc_sequential()`.

- bins:

  Number of bins of the histograms. One whole number, 1 or more.

## Value

An `mc_abc_sequential`: a list with `rounds`, one element per round with
its `table` (an `mc_abc_table`), `fit` (an `mc_abc_fit`) and `weights`
(a data frame with the parameters of the simulations kept by the fit and
one column of normalised weights per parameter, named `w_<parameter>`);
and `parameters`, `statistics`, `prior` and `observed`.

A data frame with one row per parameter: `round`, `parameter`,
`estimate` (weighted median), `low` and `high` (the limits of the
central credible interval), `mean` (weighted mean) and `ess`, the
effective number of simulations behind the weights, 1 / sum(weights^2).

## Details

The gain is in precision for the same number of simulations when the
prior is wide and the posterior narrow, since later rounds simulate
where the posterior is. The result belongs to one observation: a new
observation needs new rounds,
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md)
and
[`abc_calibrate()`](https://provetelab.org/comdynamite/reference/abc_calibrate.md)
do not apply, and the coverage of the intervals can only be measured by
repeating the whole procedure on pseudo-observations
(`data-raw/abc_sequential.R`). As in
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
each parameter has its own forest and weights, so the proposal draws
each parameter on its own and the correlation among parameters is not
kept. The number of rounds is fixed: compare the posteriors of the last
rounds, with predict() and its `round`, to see whether more would help.

The rounds can make the intervals too narrow. The posterior of a forest
rests on few simulations, its effective number (`ess` of predict()), and
its error of location follows the spread of the reference table rather
than the width of the posterior. As the rounds concentrate the table,
this error grows relative to the interval and builds up from round to
round. In a toy model with a known posterior and one informative
statistic (`data-raw/abc_sequential.R`), about 14 effective simulations
per round brought the coverage of 90% intervals from 0.99 in the first
round to 0.90 in the third and 0.85 in the fourth, while the same rounds
with exact weights kept it at 0.92; wider kernels restored the coverage
only by giving back the narrower intervals
(`data-raw/abc_sequential_diagnostics.R`). In the community model of
`data-raw/abc_recovery.R`, with 117 to 161 effective simulations, the
coverage held. A warning is given when the last round has fewer than 50
effective simulations for some parameter; there, prefer fewer rounds and
compare them.

## References

Beaumont, M. A., Cornuet, J.-M., Marin, J.-M. & Robert, C. P. (2009).
Adaptive approximate Bayesian computation. *Biometrika* 96: 983-990.
[doi:10.1093/biomet/asp052](https://doi.org/10.1093/biomet/asp052)

Dinh, K. N., Liu, C., Xiang, Z., Liu, Z. & Tavaré, S. (2025).
Approximate Bayesian computation sequential Monte Carlo via random
forests. *Statistics and Computing* 35: 219.
[doi:10.1007/s11222-025-10748-x](https://doi.org/10.1007/s11222-025-10748-x)

## See also

Other inference:
[`abc_calibrate()`](https://provetelab.org/comdynamite/reference/abc_calibrate.md),
[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md),
[`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md),
[`abc_diagnose()`](https://provetelab.org/comdynamite/reference/abc_diagnose.md),
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
[`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md),
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md),
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
[`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md),
[`nobs.mc_node_test()`](https://provetelab.org/comdynamite/reference/nobs.mc_node_test.md),
[`observe_detection()`](https://provetelab.org/comdynamite/reference/observe_detection.md),
[`plot.mc_abc_choice()`](https://provetelab.org/comdynamite/reference/plot.mc_abc_choice.md),
[`plot.mc_abc_fit()`](https://provetelab.org/comdynamite/reference/plot.mc_abc_fit.md),
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)

## Examples

``` r
# one statistic, the parameter plus noise, and a prior much wider than the posterior; with
# so few simulations the forests rest on a handful of them, hence the warning
toy <- function(scenario, shared) scenario$a + rnorm(1, 0, 0.1)
prior <- list(a = c(-10, 10))
seq_fit <- abc_sequential(prior, toy, observed = c(s = 2), n = 200, rounds = 3,
                          summarise = function(sim, shared) c(s = sim), ntree = 100)
#> abc_sequential() implements the method of:
#>   Dinh, K. N., Liu, C., Xiang, Z., Liu, Z., & Tavaré, S. (2025). Approximate Bayesian computation sequential Monte Carlo via random forests. Statistics and Computing, 35, 219. https://doi.org/10.1007/s11222-025-10748-x
#>   Beaumont, M. A., Cornuet, J.-M., Marin, J.-M., & Robert, C. P. (2009). Adaptive approximate Bayesian computation. Biometrika, 96, 983-990. https://doi.org/10.1093/biomet/asp052
#>   Raynal, L., Marin, J.-M., Pudlo, P., Ribatet, M., Robert, C. P., & Estoup, A. (2019). ABC random forests for Bayesian parameter inference. Bioinformatics, 35, 1720-1728. https://doi.org/10.1093/bioinformatics/bty867
#> Please cite these papers together with comdynamite (citation("comdynamite")).
#> This message is shown once per session; options(comdynamite.cite = FALSE) turns it off.
#> Warning: abc_sequential(): Fewer than 50 effective simulations in the last round for a; the intervals may be too narrow. Compare the rounds with predict(round = ).
seq_fit
#> Sequential ABC with random forests: 3 rounds, 600 simulations used
#> Parameters: a 
#> Statistics: 1 
#>  parameter estimate  low high  ess
#>          a     1.94 1.88 2.05 5.21
#> Effective simulations per round:
#>           a
#> round 1 5.9
#> round 2 3.7
#> round 3 5.2
predict(seq_fit, round = 1)
#>   round parameter estimate     low     high     mean      ess
#> 1     1         a 1.975219 1.91424 2.106069 1.992312 5.921633
```
