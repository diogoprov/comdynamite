# Species traits, with or without phylogenetic signal

Environmental optima evolve by Brownian motion with Pagel's lambda and
are rescaled to `optima_range` (or, with `uniform_optima = TRUE`,
replaced by ranks with evenly spaced values on that range, so that only
their arrangement on the tree depends on the tree). Niche breadths are
drawn from a normal distribution truncated at zero, and maximum growth
rates from a uniform distribution. Either can carry phylogenetic signal:
the values are then assigned to species by the ranks of a
Brownian-motion trait, which keeps their marginal distribution.

## Usage

``` r
species_traits(
  species,
  lambda = 1,
  optima_range = c(0, 100),
  uniform_optima = FALSE,
  sigma_mean = 10,
  sigma_sd = 10,
  rmax_range = c(2, 8),
  sigma_signal = FALSE,
  rmax_signal = FALSE
)
```

## Arguments

- species:

  A `phylo` tree, or the number of species for traits without
  phylogenetic structure.

- lambda:

  Pagel's lambda of the optima (0 = no signal, 1 = Brownian motion).
  Ignored when `species` is a number.

- optima_range:

  Range of the optima.

- uniform_optima:

  Give the optima evenly spaced values by rank?

- sigma_mean, sigma_sd:

  Mean and standard deviation of the normal distribution (truncated at
  zero) of niche breadths.

- rmax_range:

  Range of the uniform distribution of maximum growth rates.

- sigma_signal, rmax_signal:

  Should niche breadth or maximum growth rate carry phylogenetic signal
  (Brownian motion)?

## Value

Data frame with columns `optimum`, `sigma` and `r_max`, one row per
species, with species names as row names.

## Examples

``` r
set.seed(1)
tr <- ape::rcoal(10)
species_traits(tr, lambda = 1, sigma_signal = TRUE)
#>       optimum     sigma    r_max
#> t6   65.34866  6.057100 4.140361
#> t4    2.75498  9.406866 4.588842
#> t9    0.00000  8.354764 2.889269
#> t8   35.42721  3.112443 2.078465
#> t2   67.23493  7.466383 6.293396
#> t7   46.93207  2.925048 2.619105
#> t1  100.00000 17.631757 4.677706
#> t3   80.94664 15.566632 5.840606
#> t10  37.44939 21.000254 7.951032
#> t5   30.26917 16.969634 4.973561
species_traits(5)
#>        optimum     sigma    r_max
#> sp1  22.027020  8.648214 7.318706
#> sp2  88.968198 21.780870 5.020037
#> sp3 100.000000 15.939462 7.262345
#> sp4   0.000000 13.329504 3.135162
#> sp5   3.141045 20.630998 6.548618
```
