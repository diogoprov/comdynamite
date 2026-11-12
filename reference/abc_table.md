# Reference table for approximate Bayesian computation

Draws parameters from uniform priors, simulates one data set per draw
and summarises it, as the first step of inference by approximate
Bayesian computation (ABC). What `setup` returns is drawn once and
shared by every simulation, so it should hold what is fixed by the data,
such as the observed tree and landscape. A parameter that spans orders
of magnitude, such as dispersal, is better drawn on the log scale and
transformed back inside `simulate`.

## Usage

``` r
abc_table(
  prior,
  simulate,
  n,
  setup = NULL,
  summarise = function(sim, shared) summary_stats(sim, shared$tree),
  seed = 1
)
```

## Arguments

- prior:

  Named list with one numeric vector of length 2 per parameter, the
  lower and upper limits of its uniform prior.

- simulate:

  Function of `scenario`, a list with one value per parameter, and
  `shared`, what `setup` returned, which gives a simulation. As in
  [`run_scenarios()`](https://provetelab.org/comdynamite/reference/run_scenarios.md).

- n:

  Number of simulations.

- setup:

  Optional function without arguments that gives what all simulations
  share.

- summarise:

  Function of `sim` and `shared` that gives a named numeric vector. The
  default is
  [`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md),
  with `shared$tree` if there is one.

- seed:

  Integer that sets all the seeds.

## Value

An `mc_abc_table`: a data frame with one row per simulation, the
parameters and the statistics, with the attributes `prior`, `parameters`
and `statistics`.

## See also

[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md)
to estimate the parameters,
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md)
to check whether they can be recovered.

Other inference:
[`abc_choose()`](https://provetelab.org/comdynamite/reference/abc_choose.md),
[`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md),
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
[`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md),
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md),
[`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md),
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)

## Examples

``` r
shared_pool <- function() {
  list(tree = treeshapesim::sim_tree_shape(10), landscape = landscape_implicit(8))
}
one_run <- function(scenario, shared) {
  simulate_metacommunity(shared$landscape, species_traits(shared$tree, scenario$lambda),
                         competition_matrix(shared$tree, "none"),
                         dispersal = 10^scenario$log_dispersal, n_steps = 30,
                         seed_steps = 1:3)
}
prior <- list(log_dispersal = c(-3, -0.3), lambda = c(0, 1))
table <- abc_table(prior, one_run, n = 20, setup = shared_pool)
head(table[, 1:5])
#>   log_dispersal    lambda alpha gamma     beta
#> 1    -2.2831266 0.9347052 2.750    10 3.636364
#> 2    -1.9952655 0.2121425 4.500     8 1.777778
#> 3    -1.4532959 0.6516738 6.250     8 1.280000
#> 4    -0.5478390 0.1255551 7.250     8 1.103448
#> 5    -2.4554588 0.2672207 3.625     8 2.206897
#> 6    -0.5743479 0.3861141 4.875     5 1.025641
```
