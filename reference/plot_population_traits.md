# Plot the trait distributions of populations within communities

Draws, in one panel per site, the trait probability density of each
population from
[`trait_densities()`](https://provetelab.org/comdynamite/reference/trait_densities.md),
filled with the colour of its species and with adjustable transparency,
the density of the community as a black line and the individuals as a
rug. With `weighted = TRUE` (the default) each population is scaled by
its relative abundance, so that the populations of a site add up to the
community line; with `weighted = FALSE` every population integrates to
one, which shows the shape of rare populations.

## Usage

``` r
plot_population_traits(
  individuals,
  trait = "trait",
  site = "site",
  species = "species",
  abundance = NULL,
  sites = NULL,
  weighted = TRUE,
  community = TRUE,
  rug = TRUE,
  alpha = 0.35,
  colours = NULL,
  bw = NULL,
  n = 512,
  ncol = NULL,
  scales = "free_y"
)
```

## Arguments

- individuals:

  Data frame with one row per individual.

- trait:

  Name of the column with the trait.

- site, species:

  Names of the columns with the site and the species.

- abundance:

  Name of a column with the abundance of each population (the same value
  for all its individuals), or `NULL` to weight populations by the
  number of individuals measured.

- sites:

  Sites to draw; `NULL` (the default) draws all.

- weighted:

  Scale the population densities by their relative abundances.

- community:

  Draw the density of the community.

- rug:

  Draw the individuals as a rug.

- alpha:

  Transparency of the population densities (0, transparent, to 1,
  opaque).

- colours:

  Colours of the species: unnamed, in the order of first appearance in
  `individuals`, or named by species. `NULL` (the default) uses the
  Okabe-Ito palette for up to seven species and
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)
  (`"Dark 3"`) for more.

- bw:

  Bandwidth of the kernel, or `NULL` for the rule above.

- n:

  Number of points of the grid of trait values.

- ncol, scales:

  Passed to
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).

## Value

A ggplot object.

## Details

The ggplot object can be changed further, for instance with
`+ ggplot2::facet_wrap(~ site, ncol = 2)`,
`+ ggplot2::scale_fill_viridis_d()` or `+ ggplot2::theme_classic()`.

## See also

[`trait_densities()`](https://provetelab.org/comdynamite/reference/trait_densities.md),
[`tstats()`](https://provetelab.org/comdynamite/reference/tstats.md).

## Examples

``` r
set.seed(1)
ind <- data.frame(site = rep(paste0("s", 1:4), each = 40),
                  species = sample(c("a", "b", "c"), 160, replace = TRUE))
ind$size <- c(a = 10, b = 14, c = 17)[ind$species] + rnorm(160)
plot_population_traits(ind, "size")

plot_population_traits(ind, "size", weighted = FALSE, alpha = 0.2, sites = c("s1", "s2"))
```
