# Run a grid of simulation scenarios

Runs every row of `scenarios` in every replicate and returns one table.
Each replicate draws what `setup` returns once, so scenarios are
compared on the same tree and landscape. Each cell of the design has its
own seed, which comes from `seed`, its replicate and its row, so results
do not depend on the order of the cells nor on the number of workers.
Cells run in parallel under a plan of the future package, if you choose
one. The random number generator of the session is left as it was.

## Usage

``` r
run_scenarios(
  scenarios,
  simulate,
  setup = NULL,
  summarise = function(sim, shared) summary_stats(sim, shared$tree),
  n_rep = 1,
  seed = 1
)
```

## Arguments

- scenarios:

  Data frame with one scenario per row and one parameter per column.

- simulate:

  Function of `scenario`, a list with one row of `scenarios`, and
  `shared`, what `setup` returned, which gives a simulation.

- setup:

  Optional function without arguments that draws what the scenarios of
  one replicate share, such as a list with `tree` and `landscape`.

- summarise:

  Function of `sim` and `shared` that gives a named numeric vector. The
  default is
  [`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md),
  with `shared$tree` if there is one.

- n_rep:

  Number of replicates.

- seed:

  Integer that sets all the seeds.

## Value

Data frame with the columns of `scenarios`, `replicate` and the
statistics.

## See also

Other simulation studies:
[`as.data.frame.mc_sim()`](https://provetelab.org/comdynamite/reference/as.data.frame.mc_sim.md),
[`equilibrium_check()`](https://provetelab.org/comdynamite/reference/equilibrium_check.md),
[`summary_stats()`](https://provetelab.org/comdynamite/reference/summary_stats.md)

## Examples

``` r
scenarios <- expand.grid(dispersal = c(0.001, 0.1), lambda = c(0, 1))
shared_pool <- function() {
  list(tree = treeshapesim::sim_tree_shape(10), landscape = landscape_implicit(8))
}
one_run <- function(scenario, shared) {
  simulate_metacommunity(shared$landscape, species_traits(shared$tree, scenario$lambda),
                         competition_matrix(shared$tree, "none"),
                         dispersal = scenario$dispersal, n_steps = 30, seed_steps = 1:3)
}
out <- run_scenarios(scenarios, one_run, shared_pool, n_rep = 2)
aggregate(alpha ~ dispersal + lambda, out, mean)
#>   dispersal lambda  alpha
#> 1     0.001      0 2.6250
#> 2     0.100      0 5.5625
#> 3     0.001      1 3.6875
#> 4     0.100      1 5.8125
```
