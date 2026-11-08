# Local adaptation and character displacement from population means

Estimates how far the populations of a species follow the environment of
their site and how far co-occurring species push them apart, by fitting
both terms of the model of
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md)
at once to the means of the populations:

## Usage

``` r
population_shifts(
  individuals,
  trait = "trait",
  site = "site",
  species = "species",
  env = "env",
  abundance = "abundance",
  width = NULL
)
```

## Arguments

- individuals:

  Data frame with one row per individual, such as from
  [`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md).

- trait, site, species, env, abundance:

  Names of the columns with the trait value, the site, the species, the
  environment of the site and the abundance of the population.

- width:

  Width `w` of the competition kernel; `NULL` (the default) chooses it
  from a grid of multiples (1/8 to 4) of the median standard deviation
  of population means within sites, refined by
  [`stats::optimize()`](https://rdrr.io/r/stats/optimize.html) between
  the neighbours of the best value of the grid. `c` is estimated at that
  width, and the two trade off when data are few.

## Value

A data frame with one row and the columns `local_adaptation` (`h`),
`displacement` (`c`), `width` (`w`), `r_squared` (of the fit within
species), `n_populations` and `n_species` (that entered as responses).
The attribute `profile` holds the residual sum of squares for each width
tried.

## Details

\$\$m\_{ij} = a_j + h E_i + c \sum_k p\_{ik} (m\_{ij} - m\_{ik})
\exp\left(-\frac{(m\_{ij} - m\_{ik})^2}{2 w^2}\right),\$\$

where \\m\_{ij}\\ is the mean of species \\j\\ at site \\i\\, \\a_j\\ a
constant of the species, \\E_i\\ the environment of the site and
\\p\_{ik}\\ the relative abundance of species \\k\\ among the species
measured at the site. In
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md),
\\a_j = (1 - h) \theta_j\\, `h` is `local_adaptation`, `c` is
`displacement` and `w` is `displacement_width`. The constants of the
species are removed by centring each variable within species, and `h`
and `c` are the coefficients of a regression through the origin; `w`,
when not given, is the value of a grid that leaves the smallest residual
sum of squares.

Local adaptation alone moves the populations of a species along the
environment, and displacement alone moves them away from the species
they meet. The two leave different signatures, but each term alone picks
up part of the other: the slope of population means on the environment
falls when displacement pulls co-occurring populations apart, and the
spread of a site rises with displacement and falls with local
adaptation. Fitting both terms together separates them when populations
of the same species meet different environments and different
competitors. The approach follows the advice of Goldberg & Lande (2006)
and Adams & Collyer (2007) to judge divergence among sympatric species
after the variation of the trait along the environment is accounted for,
here with the strength of competition at each site in place of a
contrast between sympatry and allopatry.

## Sampling error

The mean of a population enters both sides of the equation, so its
sampling error would correlate the push of competitors with the response
and inflate `c`. The individuals of each population are therefore split
at random in two halves (so the estimates change slightly with the seed
of the random number generator): the response is the mean of one half,
and the own mean in the push is that of the other half (competitors
enter with their full means). Both assignments are stacked in one fit.
Only populations with at least two individuals, of species with at least
two such populations, enter as responses; all populations enter as
competitors. Sampling error left in the push still biases `c` toward
zero, more so with few individuals per population and with strong
displacement: strongly displaced populations sit in the tail of the
kernel of their competitors, where the push is small but its sampling
error is not. Competitors are the populations present in `individuals`,
so species that occur at a site but were not measured (as with
`fraction` in
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md))
are missing from the push and bias `c` as well. The estimates are best
read as summary statistics, with their uncertainty and bias found by
simulating the design, as in approximate Bayesian computation.

Without displacement the width is not identified: the residual sum of
squares changes little along the grid, and a narrow kernel, which acts
on few pairs of populations, can give a value of `c` far from zero, of
either sign. Look at the attribute `profile`, or fix `width` from
independent knowledge of the niche widths, before reading `c` as
evidence of displacement.

## References

Adams, D. C., & Collyer, M. L. (2007). Analysis of character divergence
along environmental gradients and other covariates. *Evolution*, 61,
510–515.
[doi:10.1111/j.1558-5646.2007.00063.x](https://doi.org/10.1111/j.1558-5646.2007.00063.x)

Goldberg, E. E., & Lande, R. (2006). Ecological and reproductive
character displacement on an environmental gradient. *Evolution*, 60,
1344–1357. [doi:10.1554/05-696.1](https://doi.org/10.1554/05-696.1)

## See also

[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md)
with `null_ip = "species"`.

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(20)
sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
                              competition_matrix(tree, "none"), n_steps = 100,
                              seed_steps = 1:5)
ind <- sample_individuals(sim, n = 20, local_adaptation = 0.5, displacement = 2)
population_shifts(ind)
#>   local_adaptation displacement    width r_squared n_populations n_species
#> 1        0.4914459     1.990583 12.13176 0.9960785            91        14
```
