# T-statistics of intraspecific trait variation

Computes the three variance ratios of Violle et al. (2012) for each site
(community) and trait, from trait values of individuals such as those of
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md):

## Usage

``` r
tstats(
  individuals,
  traits = "trait",
  site = "site",
  species = "species",
  method = c("base", "cati"),
  nperm = 99,
  ...
)
```

## Arguments

- individuals:

  Data frame with one row per individual, such as from
  [`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md).

- traits:

  Names of the columns with trait values.

- site, species:

  Names of the columns with the site and the species of each individual.

- method:

  `"base"` (the default) for the observed values only, or `"cati"` for
  [`cati::Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html)
  with null models; it needs the package cati.

- nperm:

  Number of randomisations of the null models of `method = "cati"`;
  `NULL` for observed values only.

- ...:

  Further arguments to
  [`cati::Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html),
  such as `reg.pool` or `independantTraits`.

## Value

With `method = "base"`, a data frame with one row per site and trait and
the columns `site`, `trait`, `n_populations`, `T_IP.IC`, `T_IC.IR` and
`T_PC.PR`; statistics that cannot be computed, such as `T_PC.PR` for a
site with one species, are `NA`. With `method = "cati"`, the `Tstats`
object of
[`cati::Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html),
with the site names of `individuals`, for its methods
[`cati::ses.Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`summary()`](https://rdrr.io/r/base/summary.html).

## Details

- `T_IP.IC`, the mean variance among individuals within the populations
  of the site over the variance among all individuals of the site. Low
  values mean that the populations of co-occurring species overlap
  little (niche differentiation within the community).

- `T_IC.IR`, the variance among individuals of the site over that among
  all individuals of the region. Low values mean that the site holds a
  narrow part of the regional trait range (environmental filtering of
  individuals).

- `T_PC.PR`, the variance among population means within the site over
  that among all population means of the region. Low values mean
  filtering of populations.

A population is a species at a site. Variances use `n - 1`, means of the
variances within populations are not weighted, and populations of a
single individual have no variance and are left out of `T_IP.IC`, as in
[`cati::Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html).

`method = "base"` computes the observed values in base R, with no null
model. It is fast and needs no other package, which suits simulation
studies with many replicates, where the truth is known and each
statistic is compared among scenarios rather than with a null model.
`method = "cati"` calls
[`cati::Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html)
(Taudiere & Violle 2016), which adds null models: individuals permuted
within each site for `T_IP.IC`, individuals drawn from the regional pool
for `T_IC.IR`, and population means permuted for `T_PC.PR`. Both give
the same observed values.
[`cati::Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html)
finds the populations of a site by matching its name as a regular
expression, so a site called `site1` would also take the populations of
`site10`; the sites are therefore passed under labels of equal width
(`site01`, `site02`, ...) and given back their names in the result.

## Local adaptation

With `local_adaptation` in
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md),
population means move toward the environment of their site, so
co-occurring species become more alike and the region spreads along the
environment. The variance among populations within a site falls, so
`T_PC.PR` and `T_IC.IR` fall, and `T_IP.IC` rises, because the
populations that remain overlap more. A fall in `T_IP.IC` with local
adaptation would need populations to diverge from co-occurring species,
as in character displacement, which
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md)
does not simulate.

## References

Taudiere, A., & Violle, C. (2016). cati: an R package using functional
traits to detect and quantify multi-level community assembly processes.
*Ecography*, 39, 699–708.
[doi:10.1111/ecog.01433](https://doi.org/10.1111/ecog.01433)

Violle, C., Enquist, B. J., McGill, B. J., Jiang, L., Albert, C. H.,
Hulshof, C., Jung, V., & Messier, J. (2012). The return of the variance:
intraspecific variability in community ecology. *Trends in Ecology &
Evolution*, 27, 244–252.
[doi:10.1016/j.tree.2011.11.014](https://doi.org/10.1016/j.tree.2011.11.014)

## See also

[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md).

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(20)
sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
                              competition_matrix(tree, "none"), n_steps = 100,
                              seed_steps = 1:5)
none <- tstats(sample_individuals(sim, intraspecific_sd = 2))
local <- tstats(sample_individuals(sim, local_adaptation = 0.8, intraspecific_sd = 2))
colMeans(none[, c("T_IP.IC", "T_IC.IR", "T_PC.PR")], na.rm = TRUE)
#>    T_IP.IC    T_IC.IR    T_PC.PR 
#> 0.02299744 0.89106250 1.11059919 
colMeans(local[, c("T_IP.IC", "T_IC.IR", "T_PC.PR")], na.rm = TRUE)
#>    T_IP.IC    T_IC.IR    T_PC.PR 
#> 0.33914930 0.03219661 0.03057455 
# with the null models of cati; ses.Tstats() lists T_IP.IC, T_IC.IR and T_PC.PR in order
res <- tstats(sample_individuals(sim, intraspecific_sd = 2), method = "cati", nperm = 19)
round(cati::ses.Tstats(res)[[1]]$ses[, 1], 1)
#>  site1  site2  site3  site4  site5  site6  site7  site8  site9 site10 site11 
#>  -15.6  -12.9  -16.9  -13.9   -8.7  -10.7   -8.7  -16.9   -7.1  -11.7  -14.8 
#> site12 site13 site14 site15 
#>  -11.2  -15.1  -10.2  -18.5 
```
