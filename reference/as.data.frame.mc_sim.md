# Simulation results as a long table

Simulation results as a long table

## Usage

``` r
# S3 method for class 'mc_sim'
as.data.frame(x, row.names = NULL, optional = FALSE, zeros = FALSE, ...)
```

## Arguments

- x:

  An `mc_sim` object.

- row.names, optional:

  Not used.

- zeros:

  Keep absences?

- ...:

  Not used.

## Value

Data frame with `site`, `species`, `step` and `abundance`, from
`history`, or from `abundance` with `step` as `NA` if no steps were
recorded.
[`tibble::as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
turns it into a tibble.

## Examples

``` r
set.seed(1)
sim <- simulate_metacommunity(landscape_implicit(4), species_traits(5),
                              competition_matrix(5, "none"), n_steps = 10,
                              seed_steps = 1:2, record_steps = c(5, 10))
head(as.data.frame(sim))
#>    site species step abundance
#> 1 site4     sp2    5        87
#> 2 site1     sp3    5         1
#> 3 site2     sp3    5        12
#> 4 site3     sp3    5        31
#> 5 site1     sp4    5         6
#> 6 site2     sp5    5        22
```
