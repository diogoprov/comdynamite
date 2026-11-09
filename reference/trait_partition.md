# Partition trait or phylogenetic variation between environment and space

The phylogenetic gradient analysis of Peres-Neto et al. (2012). Trait
values are standardised with species weights, which gives each site a
mean and a standard deviation. The chosen component is related to `env`
and `space` with site weights, and the variation is split into the
fractions of Peres-Neto et al. (2006): `a`, environment alone; `b`,
shared; `c`, space alone; and `d`, unexplained.

## Usage

``` r
trait_partition(
  comm,
  traits,
  env,
  space,
  component = c("mean", "variance"),
  n_perm = 999,
  site_permutations = NULL,
  shared_null = c("none", "msr", "sites"),
  coords = NULL,
  level = 0.05
)
```

## Arguments

- comm:

  Matrix with one row per site and one column per species. Every site
  needs at least one individual. Or an
  [`mc_data()`](https://provetelab.org/comdynamite/reference/mc_data.md)
  object, whose parts fill the arguments left out.

- traits:

  Numeric vector or matrix with one row per species, which can hold
  traits or the vectors of
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md).
  Species are matched to `comm` by name, in any order, and those in only
  one of the two are dropped with a warning.

- env, space:

  Numeric vector, or matrix or data frame, with one row per site. For
  `space`, the article uses Moran's eigenvector maps, which
  [`adespatial::dbmem()`](http://adeverse.github.io/adespatial/reference/dbmem.md)
  gives.

- component:

  `"mean"` for the composition of sites or `"variance"` for their
  dispersion.

- n_perm:

  Number of permutations of sites and of species.

- site_permutations:

  Optional matrix of site indices with one permutation per row, for
  example from
  [`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md),
  to use instead of free permutations of sites. Each set of predictors
  is then moved alone in its tests (see Details). Empty sites are
  allowed: they leave the analysis but their predictors can reach
  occupied sites.

- shared_null:

  Null model of the test of the shared fraction `b`: `"none"` leaves it
  untested, `"msr"` replaces `env` by Moran spectral randomisations
  built from `coords` with
  [`adespatial::msr()`](http://adeverse.github.io/adespatial/reference/msr.md),
  and `"sites"` moves `env` alone with `site_permutations`, for a
  complete grid. `"msr"` needs the package adespatial.

- coords:

  Matrix of site coordinates with one row per site of `comm`, for
  `shared_null = "msr"`. The spatial weights are those of
  [`adespatial::dbmem()`](http://adeverse.github.io/adespatial/reference/dbmem.md).

- level:

  Significance level of the test of `env` against `space` that decides
  whether `b` is tested.

## Value

List of class `mc_trait_partition` with `fractions`, a data frame with
one row for each of `abc`, `ab`, `bc`, `a`, `b`, `c` and `d`, and the
columns `r2`, `r2_adjusted`, `r2_component`, `statistic`, the value
tested, `p_site`, `p_species` and `p_max`, which are `NA` for `d`, for
`b` without `shared_null`, and for a fraction that cannot be estimated,
as when a set of predictors has as many columns as sites; `component`;
`component_share`, the share of the total variation in the component;
`semipartial`, a matrix of predictors by traits; `p_partial_max`,
`p_partial_site` and `p_partial_species`, with the larger p-value and
those of sites and of species, all `NA` if predictors are collinear;
`shared`, with `null`, `draws`, `p_env_space`, the p-value of `env`
against `space`, and `tested`; and `permutations`, whose `moved` is
`"together"` with free permutations and `"each set alone"` with
`site_permutations`.

## Details

The total variation includes both the mean and the dispersion within
sites, so `r2` answers how much of all the variation each fraction
holds, and `r2_component` how much of the component. Each fraction is
tested by permutations of sites, which move the predictors of all sites
together, and by permutations of the trait values among species, and
`p_max` is the larger p-value, as in
[`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md).
The statistic of `a` and `c` is the share of the component that the
other set leaves unexplained, \\a / (a + d)\\ within the component,
which is a monotone function of the partial F statistic of the set given
the other (Anderson & Legendre 1999). The fraction itself, which the
original code tests, is not, and its test loses most of its power when
the other set explains much: see the section on simulations. `d` cannot
be tested.

With `site_permutations`, the moves keep the spatial structure of what
they move, so each set is moved alone while the response and the other
set stay in place: `env` for `ab`, `a` and the semipartials of its
predictors, and `space` for `bc`, `c` and the semipartials of its
predictors; `abc` moves both. Moving both sets together, as free
permutations do, tests `c` against an environment that no longer
explains the response, and on a grid that test was liberal when the
environment was strong and had no power when it was spatially structured
(see the simulations). The result records this in `permutations$moved`.

The shared fraction `b` has no test by these permutations, because a
permutation of sites breaks the relation between environment and space,
and is tested with `shared_null` by the scheme of
[`adespatial::envspace.test()`](http://adeverse.github.io/adespatial/reference/envspace.test.md)
(Bauman et al. 2019): the response and `space` stay fixed, `env` is
replaced by surrogates that keep its spatial autocorrelation, and `b` is
computed again. Moran spectral randomisation (Wagner & Dray 2015) builds
them from `coords` for any design, and `site_permutations` from
[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md)
move the whole map of `env` on a complete grid. The p-value of species
comes from the same permutations of species as the other fractions. A
positive `b` can arise from two spatially structured variables without
relation, which free permutations do not take into account (Clappe et
al. 2018). The test asks whether `b` is larger than expected, so a
negative `b` is never significant. It is only meaningful when `env`
holds some of the structure in `space`, so `env` is first tested against
`space` with free permutations of sites, and when that test is not
significant at `level`, `b` is left untested and `shared$tested` is
`FALSE`.

## Note

Adapted from the MATLAB function `PartitionRDATwoDataMatricesPhylogeny`
in the supplement of Peres-Neto et al. (2012).

## Simulations

`data-raw/validation_fractions.R` and the article on validation measure
the tests on a grid of 8 by 8 sites, with 500 data sets per scenario and
199 permutations, at the 5% level. Where space acted and the environment
did not, the fraction `a` itself was never rejected, and the share left
by space was rejected in 5.6% to 5.8% of data sets with free
permutations of sites and in 4.2% to 4.4% with the larger p-value; where
the environment had a weak effect beside space, the fraction was
rejected in 15% and the share in 100%. Two limits: with an environment
structured at scales that `space` does not hold, free permutations
rejected `a` in 11.4% (moves of the map with
[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md):
1.4%, larger p-value: 4.0%); and with a strong, spatially structured
environment and no effect of space, free permutations rejected `c` in
10% (larger p-value: 2.4%), where the residuals of the linear fit of the
environment keep spatial structure. Up to version 0.6.1,
`site_permutations` moved both sets together in every test, and that
test of `c` rejected 19% with a strong environment without structure and
never with a structured one.

`data-raw/validation_fractions_null.R` chose what the moves of
`site_permutations` should move, with
[`torus_permutations()`](https://provetelab.org/comdynamite/reference/torus_permutations.md),
on the first generator (300 data sets) and checked the choice on the
second (200 data sets). For `a`, moving `env` alone held the level
wherever the environment had no effect, in both generators (3.7% to 5.7%
and 4.5% to 7.1% with the p-value of sites, standard errors of 0.011 to
0.018), including an environment structured at scales outside `space`,
where free permutations rejected 9.7% and 20%, and its power was the
same. When the environment acted on a trait other than the one analysed,
the p-value of sites rejected 45% and 75%, as with free permutations,
and `p_max` 2% and 4%. For `c`, moving `space` alone held the level
(3.0% to 7.0% and 3.5% to 8.0%; larger p-value at most 3%) where free
permutations rejected 8.0% and 24% with a strong, structured environment
acting alone, and moving both sets together rejected 22% with a strong
environment without structure. Its power depends on the question: it was
complete when the species followed a field in the span of `space`, and
25% to 28% when they followed a field of exponential covariance that
`space` holds only in part (free permutations: 98%), because the maps
moved around the grid explain such a field about as well as the maps in
place. The restricted test of `c` asks whether the maps in `space`
explain more than the same maps elsewhere on the grid, not whether the
remainder has any spatial structure. For that question, use free
permutations and read `p_max`.

These limits were larger with a second generator
(`data-raw/validation_fractions_check.R`: counts, a grid of 10 by 10
sites, 8 Moran's eigenvector maps, 200 data sets) where a hidden field
that `space` holds only in part drove the species: free permutations
rejected `a` in 21% (larger p-value: 2%), and with a structured
environment acting alone they rejected `c` in 25% (larger p-value:
7.5%). The p-value of sites alone should not be read when residuals can
keep spatial structure; `p_max` should. Without effects, the rates there
were 4.5% to 6.1% with free permutations, and Moran spectral
randomisation rejected `b` in 5% to 6%.

For `b`, with an environment that was spatially structured but had no
effect, free permutations of `env` rejected 48% of data sets, Moran
spectral randomisation 5.6% and moves of the map 11%. With a trait
unrelated to the composition, the randomisation alone rejected 54% and
the larger p-value 4.8%: the permutation of species is needed, as for
the other fractions. Power was 100% when the environment alone acted,
43% when space also acted on its own, and 28% with a weak environment
beside space.

## Article and original code

Raw and adjusted fractions of the total for the mean component agree
with `PartitionRDATwoDataMatricesPhylogeny` of the supplement, checked
in Octave. The code departs from Appendix A in these points, where we
follow the text:

- Fractions relative to the component: `Meana` and `Meanc` subtract the
  fractions of the total from the fraction of the component.

- Adjusted values relative to the component use the number of predictors
  of both sets in the models of one set, and the result can exceed one.

- Total variation: Appendix A multiplies it by the number of traits, and
  the code does not.

- Several vectors: the article analyses the vectors of environment and
  space together in a redundancy analysis, and the script passes only
  the first.

- Variance component: the explained sum of squares contains the
  intercept, so `b` absorbs the mean dispersion. Appendix A asks for the
  intercept in the models. Here it is removed, and `r2` is the variation
  around the weighted mean.

- Permutations of the variance component: the code takes the square root
  twice, sets small values to zero in the observed matrix instead of the
  permuted one, and rescales the permuted values to the observed mean
  and standard deviation, which the text does not mention. It also calls
  `WeightedMeanVar`, which the supplement does not include. None of this
  is done here.

- The code runs `n_perm` - 1 permutations in each null distribution.

- The text stacks two null distributions, and here we report the larger
  of the two p-values, as in
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md).

- Tests of `a` and `c`: the code permutes environment and space with the
  same order of sites, as here, and compares the fractions themselves,
  which depend on how much the other set explains. It computes the
  partial F statistics but does not use them. Here the statistic is the
  share of the component left by the other set, which is monotone in
  that F.

Moran's eigenvector maps of the function `MEM` of the supplement are the
same as those of
[`adespatial::dbmem()`](http://adeverse.github.io/adespatial/reference/dbmem.md),
checked in Octave. Appendix A recommends centring them with site weights
and the code does not, but orthogonality is not needed to partition
variation, and here all predictors are standardised with site weights.

Semipartial correlations follow Appendix A: the square root of the loss
of R2 when a predictor leaves the model of all predictors, with the sign
of its slope when it is alone. The supplement has no code for them, so
they were checked against separate regressions and not in Octave.
Appendix A calculates the model without the predictor with the slopes of
the full model, which it had just fitted again without it, and here we
use the second slopes. With several traits it gives no rule, so here
each trait has its own value, relative to the variation of that trait.
They are tested as `a` and `c` are, with the loss of R2 relative to what
the model of all predictors leaves unexplained in that trait as the
statistic, and with many predictors they need a correction.

This function partitions variation and does not test each trait against
each predictor. For that test, ter Braak (2019) proposed multilevel and
weighted averaging regressions combined by the max test, in the package
[TraitEnvMLMWA](https://github.com/CajoterBraak/TraitEnvMLMWA).

## References

Anderson, M. J., & Legendre, P. (1999). An empirical comparison of
permutation methods for tests of partial regression coefficients in a
linear model. *Journal of Statistical Computation and Simulation*, 62,
271–303.
[doi:10.1080/00949659908811936](https://doi.org/10.1080/00949659908811936)

Bauman, D., Vleminckx, J., Hardy, O. J., & Drouet, T. (2019). Testing
and interpreting the shared space-environment fraction in variation
partitioning analyses of ecological data. *Oikos*, 128, 274–285.
[doi:10.1111/oik.05496](https://doi.org/10.1111/oik.05496)

Clappe, S., Dray, S., & Peres-Neto, P. R. (2018). Beyond neutrality:
disentangling the effects of species sorting and spurious correlations
in community analysis. *Ecology*, 99, 1737–1747.
[doi:10.1002/ecy.2376](https://doi.org/10.1002/ecy.2376)

Peres-Neto, P. R., Legendre, P., Dray, S., & Borcard, D. (2006).
Variation partitioning of species data matrices: estimation and
comparison of fractions. *Ecology*, 87, 2614–2625.
[doi:10.1890/0012-9658(2006)87\[2614:VPOSDM\]2.0.CO;2](https://doi.org/10.1890/0012-9658%282006%2987%5B2614%3AVPOSDM%5D2.0.CO%3B2)

Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the
effects of spatial contingency and environmental filtering on
metacommunity phylogenetics. *Ecology*, 93, S14–S30.
[doi:10.1890/11-0494.1](https://doi.org/10.1890/11-0494.1)

ter Braak, C. J. F. (2019). New robust weighted averaging- and
model-based methods for assessing trait-environment relationships.
*Methods in Ecology and Evolution*, 10, 1962–1971.
[doi:10.1111/2041-210X.13278](https://doi.org/10.1111/2041-210X.13278)

Wagner, H. H., & Dray, S. (2015). Generating spatially constrained null
models for irregularly spaced data using Moran spectral randomization
methods. *Methods in Ecology and Evolution*, 6, 1169–1178.
[doi:10.1111/2041-210X.12407](https://doi.org/10.1111/2041-210X.12407)

## See also

[`plot_trait_partition()`](https://provetelab.org/comdynamite/reference/plot_trait_partition.md)

## Examples

``` r
set.seed(1)
coords <- cbind(1:40)
env <- runif(40, 0, 100)
space <- as.matrix(adespatial::dbmem(coords, silent = TRUE))[, 1:3]
trait <- setNames(runif(15, 0, 100), paste0("sp", 1:15))
comm <- outer(env, trait, function(e, t) exp(-(e - t)^2 / 800))
colnames(comm) <- names(trait)
out <- trait_partition(comm, trait, env, space, n_perm = 49)
out
#> <mc_trait_partition> mean component, 0.511 of the total variation
#> permutations: 49 free of sites, 49 of species; smallest possible p: 0.02
#>   fraction       r2 r2_adjusted r2_component statistic p_site p_species p_max
#> 1      abc 0.490557      0.4323      0.96059    0.4906   0.02      0.02  0.02
#> 2       ab 0.489712      0.4763      0.95894    0.4897   0.02      0.02  0.02
#> 3       bc 0.020911     -0.0607      0.04095    0.0209   0.58      0.02  0.58
#> 4        a 0.469647      0.4930      0.91965    0.9589   0.02      0.14  0.14
#> 5        b 0.020065     -0.0167      0.03929    0.0201     NA        NA    NA
#> 6        c 0.000846     -0.0439      0.00166    0.0403   0.62      0.28  0.62
#> 7        d 0.509443      0.5677      0.03941        NA     NA        NA    NA
plot_trait_partition(out)

# the shared fraction, with an environment that follows the transect
env <- 50 + 30 * sin(1:40 / 6) + rnorm(40, 0, 5)
comm <- outer(env, trait, function(e, t) exp(-(e - t)^2 / 800))
colnames(comm) <- names(trait)
trait_partition(comm, trait, env, space, n_perm = 49, shared_null = "msr", coords = coords)
#> <mc_trait_partition> mean component, 0.404 of the total variation
#> permutations: 49 free of sites, 49 of species; smallest possible p: 0.02
#> shared fraction b: 49 Moran spectral randomisations of env (env against space: p = 0.02)
#>   fraction      r2 r2_adjusted r2_component statistic p_site p_species p_max
#> 1      abc 0.39600      0.3270      0.97986     0.396   0.02      0.02  0.02
#> 2       ab 0.39484      0.3789      0.97697     0.395   0.02      0.02  0.02
#> 3       bc 0.35841      0.3049      0.88684     0.358   0.02      0.02  0.02
#> 4        a 0.03759      0.0220      0.09302     0.822   0.02      0.14  0.14
#> 5        b 0.35724      0.3569      0.88395     0.357   0.02      0.02  0.02
#> 6        c 0.00117     -0.0519      0.00289     0.125   0.12      0.70  0.70
#> 7        d 0.60400      0.6730      0.02014        NA     NA        NA    NA
```
