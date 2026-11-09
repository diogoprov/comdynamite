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
  rmax_signal = FALSE,
  optima = NULL,
  sigma = NULL,
  r_max = NULL,
  rescale_optima = TRUE,
  nodes = FALSE
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

- optima, sigma, r_max:

  Optional values of each trait, named by species or in the order of the
  species, as a vector or a one-column matrix or data frame. They
  replace the values drawn for that trait.

- rescale_optima:

  Rescale external optima to `optima_range`?

- nodes:

  Also give the optima at the internal nodes of the tree, as
  [`trait_bm()`](https://provetelab.org/comdynamite/reference/trait_bm.md)
  does? They are drawn after all other values, so the traits are the
  same as with `nodes = FALSE` and the same seed, and rescaled as the
  optima of the tips, so they may fall outside `optima_range`. Not
  available with external optima or `uniform_optima`.

## Value

Data frame with columns `optimum`, `sigma` and `r_max`, one row per
species, with species names as row names. With a tree, the attribute
`tree` holds it, and with `nodes = TRUE` the attribute `nodes` holds the
optima at its internal nodes.

## Details

Any of the three traits can come from elsewhere instead, such as optima
simulated by mvMORPH or measured in the field, through `optima`, `sigma`
and `r_max`. The other traits are still drawn, and with the same seed
they are the same as without the external ones. External optima are
rescaled to `optima_range`, or ranked with `uniform_optima`, unless
`rescale_optima = FALSE`.

## See also

Other species traits:
[`arrival_ages()`](https://provetelab.org/comdynamite/reference/arrival_ages.md),
[`as_herodotools()`](https://provetelab.org/comdynamite/reference/as_herodotools.md),
[`food_web_rates()`](https://provetelab.org/comdynamite/reference/food_web_rates.md),
[`history_rates()`](https://provetelab.org/comdynamite/reference/history_rates.md),
[`paint_regimes()`](https://provetelab.org/comdynamite/reference/paint_regimes.md),
[`partner_turnover()`](https://provetelab.org/comdynamite/reference/partner_turnover.md),
[`plot.rate_matrix()`](https://provetelab.org/comdynamite/reference/plot.rate_matrix.md),
[`plot_rates()`](https://provetelab.org/comdynamite/reference/plot_rates.md),
[`plot_traits()`](https://provetelab.org/comdynamite/reference/plot_traits.md),
[`rate_matrix()`](https://provetelab.org/comdynamite/reference/rate_matrix.md),
[`region_pool()`](https://provetelab.org/comdynamite/reference/region_pool.md),
[`tip_history_metrics()`](https://provetelab.org/comdynamite/reference/tip_history_metrics.md),
[`trait_bm()`](https://provetelab.org/comdynamite/reference/trait_bm.md),
[`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md),
[`trait_ou()`](https://provetelab.org/comdynamite/reference/trait_ou.md),
[`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md),
[`turnover_rates()`](https://provetelab.org/comdynamite/reference/turnover_rates.md)

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
# optima from another model, here an Ornstein-Uhlenbeck process
ou <- trait_ou(tr, alpha = 1, theta = 0)
species_traits(tr, optima = ou[, 1])
#>       optimum     sigma    r_max
#> t1   95.73952 20.537509 6.528234
#> t2  100.00000 13.356172 5.722460
#> t3   53.57144 14.947958 3.017461
#> t4   61.74475 11.380527 2.373284
#> t5   29.88675  8.812080 2.654176
#> t6   55.81391 11.976843 4.290298
#> t7   26.41955  1.967868 3.015865
#> t8   55.07890 25.800917 3.791915
#> t9    0.00000 24.978188 3.153257
#> t10  16.88222 12.626455 3.543020
```
