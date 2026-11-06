# Environment that changes over time

Adds to the environment of each site a deviation that follows an
independent first-order autoregressive process with lag-one
autocorrelation `rho` and stationary standard deviation `sd`, starting
from its stationary distribution. This is a simple generator of temporal
autocorrelation, not the one used by Thompson et al. (2020).

## Usage

``` r
env_over_time(env, n_steps, rho = 0.9, sd = 5)
```

## Arguments

- env:

  Environment of each site, e.g. the `env` of an `mc_landscape`.

- n_steps:

  Number of time steps.

- rho:

  Lag-one autocorrelation of the deviations, in \[0, 1).

- sd:

  Stationary standard deviation of the deviations.

## Value

Matrix with one row per site and one column per time step, for the
`env_time` argument of
[`simulate_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_metacommunity.md).

## Examples

``` r
set.seed(1)
e <- env_over_time(c(20, 50, 80), n_steps = 100, rho = 0.95, sd = 5)
matplot(t(e), type = "l", lty = 1, xlab = "step", ylab = "environment")
```
