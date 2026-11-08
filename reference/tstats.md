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
  alternative = c("less", "greater", "two.sided"),
  null_ip = c("site", "species"),
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

  `"base"` (the default), computed in base R, or `"cati"` for
  [`cati::Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html);
  it needs the package cati.

- nperm:

  Number of randomisations of the null models; 0 or `NULL` for the
  observed values only, much faster in simulation studies with many
  replicates.

- alternative:

  Direction of the test of `method = "base"`: `"less"` (the default),
  `"greater"` or `"two.sided"`.

- null_ip:

  Null model of `T_IP.IC` for `method = "base"`: `"site"` (the default),
  trait values permuted within each site as in Taudiere & Violle (2016),
  or `"species"`, the populations of each species exchanged among its
  sites (see the section Null models).

- ...:

  Further arguments to
  [`cati::Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html),
  such as `reg.pool` or `independantTraits`.

## Value

With `method = "base"`, a data frame with one row per site and trait and
the columns `site`, `trait`, `n_populations`, `T_IP.IC`, `T_IC.IR` and
`T_PC.PR`; statistics that cannot be computed, such as `T_PC.PR` for a
site with one species, are `NA`. With `nperm > 0`, also `ses_` and `p_`
columns for each statistic, and the attribute `null`, an array
(randomisation x row x statistic) with the null values. With
`method = "cati"`, the `Tstats` object of
[`cati::Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html),
with the site names of `individuals`.

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
Individuals with a missing value of a trait are left out for that trait.

## Null models

With `nperm > 0`, `method = "base"` compares each statistic with a null
distribution, from the three null models of Taudiere & Violle (2016),
implemented in base R from their description:

- `T_IP.IC`: trait values are permuted among the individuals of each
  site, which keeps the variance of the site and breaks the structure of
  populations within it;

- `T_IC.IR`: each site receives as many individuals as it has, drawn
  without replacement from all individuals of the region;

- `T_PC.PR`: each individual takes the mean of its population, and each
  site then receives as many of these values as it has individuals,
  drawn without replacement from the region; population means are
  computed again from the values drawn.

The standardised effect size is
`(observed - mean of null) / sd of null`, and the p-value is
`(1 + r) / (1 + nperm)`, where `r` is the number of null values at or
below the observed one (`alternative = "less"`), at or above it
(`"greater"`), or twice the smaller of the two (`"two.sided"`, at most
1). Low values of the three statistics are the usual hypotheses (niche
differentiation within communities for `T_IP.IC`, environmental
filtering for the others), so `"less"` tests them. Each trait is
randomised on its own.

The null model of `T_IP.IC` asks only whether populations differ:
permuting individuals within a site removes its populations, so null
values are near 1, and co-occurring species with different means give a
large negative effect size even when each population sits where its
species would be at any site. With `null_ip = "species"`, the
populations of each species, with their individuals, are instead
exchanged at random among the sites where the species occurs. Each site
keeps its species and each species its populations; only which
population of a species sits at which site is broken. Whole populations
are moved, rather than their means, because means of populations with
different numbers of individuals are not exchangeable. The standardised
effect size is near zero when the position of a population does not
depend on the site, positive when co-occurring populations are closer
than expected from their species (convergence, as under local adaptation
to a shared environment), and negative when they are further apart
(divergence, as under character displacement);
`alternative = "two.sided"` tests both. Species found at a single site
stay in place. Local adaptation and character displacement move this
statistic in opposite directions and can cancel; see
[`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md)
to estimate both.

`method = "cati"` calls
[`cati::Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html)
(Taudiere & Violle 2016) instead and returns its object. Both methods
give the same observed values. Two problems of cati 0.99.6 affect it:
[`cati::Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.html)
finds the populations of a site by matching its name as a regular
expression, so a site called `site1` would also take the populations of
`site10` (<https://github.com/adrientaudiere/cati/issues/11>); the sites
are therefore passed under labels of equal width (`site01`, `site02`,
...) and given back their names in the result. And its null model for
`T_IP.IC` puts together the permutations of different sites in the wrong
order, so that the null values of a site come from the individuals of
other sites; its standardised effect sizes for `T_IP.IC` are then wrong
(<https://github.com/adrientaudiere/cati/issues/12>). `method = "base"`
is recommended.

## Local adaptation and character displacement

With `local_adaptation` in
[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md),
population means move toward the environment of their site, so
co-occurring species become more alike and the region spreads along the
environment. The variance among populations within a site falls, so
`T_PC.PR` and `T_IC.IR` fall, and `T_IP.IC` rises, because the
populations that remain overlap more. With `displacement`, co-occurring
populations diverge from one another instead, which raises the variance
among populations within sites and lowers `T_IP.IC`.

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

[`sample_individuals()`](https://provetelab.org/comdynamite/reference/sample_individuals.md),
[`population_shifts()`](https://provetelab.org/comdynamite/reference/population_shifts.md).

## Examples

``` r
set.seed(1)
tree <- treeshapesim::sim_tree_shape(20)
sim <- simulate_metacommunity(landscape_implicit(15), species_traits(tree),
                              competition_matrix(tree, "none"), n_steps = 100,
                              seed_steps = 1:5)
ratios <- c("T_IP.IC", "T_IC.IR", "T_PC.PR")
none <- tstats(sample_individuals(sim, intraspecific_sd = 2), nperm = 0)
local <- tstats(sample_individuals(sim, local_adaptation = 0.8, intraspecific_sd = 2),
                nperm = 0)
apart <- tstats(sample_individuals(sim, displacement = 2, intraspecific_sd = 2), nperm = 0)
rbind(none = colMeans(none[ratios], na.rm = TRUE),
      local = colMeans(local[ratios], na.rm = TRUE),
      apart = colMeans(apart[ratios], na.rm = TRUE))
#>          T_IP.IC    T_IC.IR    T_PC.PR
#> none  0.02299744 0.89106250 1.11059919
#> local 0.33914930 0.03219661 0.03057455
#> apart 0.01278061 0.91639355 1.13307743
# standardised effect sizes and p-values from the null models
res <- tstats(sample_individuals(sim, intraspecific_sd = 2), nperm = 99)
head(res[c("site", "T_IP.IC", "ses_T_IP.IC", "p_T_IP.IC")])
#>    site     T_IP.IC ses_T_IP.IC p_T_IP.IC
#> 1 site1 0.016376543  -13.545376      0.01
#> 2 site2 0.079715195  -22.236577      0.01
#> 3 site3 0.012470128  -10.490671      0.01
#> 4 site4 0.019383034  -23.070752      0.01
#> 5 site5 0.056784711   -6.326309      0.01
#> 6 site6 0.005298958   -9.945998      0.01
```
