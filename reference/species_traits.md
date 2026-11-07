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
tr <- treeshapesim::sim_tree_shape(10)
species_traits(tr, lambda = 1, sigma_signal = TRUE)
#>        optimum     sigma    r_max
#> t1   81.014935  5.218499 4.731647
#> t2   87.432733  9.461950 4.460504
#> t3   99.189329 14.179416 6.865221
#> t4  100.000000  8.972123 5.629600
#> t5   54.609559  5.850054 5.928344
#> t6    0.000000  9.406866 4.119184
#> t7   29.662989  6.057100 3.621561
#> t8   38.305851 23.586796 7.956104
#> t9   35.050894 13.876716 5.800960
#> t10   5.656532 21.000254 3.279249
species_traits(5)
#>       optimum      sigma    r_max
#> sp1  29.94562  0.6590237 2.382851
#> sp2 100.00000 12.9144624 4.012925
#> sp3  41.12514  5.5670813 6.342356
#> sp4   0.00000 10.0110535 4.025692
#> sp5  45.81600 10.7434132 5.782485
```
