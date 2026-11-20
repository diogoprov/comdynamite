# Competition matrix

Intraspecific coefficients are 1 and interspecific coefficients lie in
\[0, `inter_max`), all multiplied by `scale`. With the default `scale`
of 0.05, as in Thompson et al. (2020), a species alone at its optimum
with `r_max = 5` reaches an equilibrium of \\(5 - 1) / 0.05 = 80\\
individuals, and `inter_max` is the strength of interspecific relative
to intraspecific competition. The matrix is symmetric.

## Usage

``` r
competition_matrix(
  species,
  type = c("none", "random", "phylo", "trait"),
  inter_max = 0.1,
  scale = 0.05,
  traits = NULL
)
```

## Arguments

- species:

  A `phylo` tree, or the number of species (not for `type = "phylo"`);
  ignored for `type = "trait"`, where the species are those of `traits`.

- type:

  Type of interspecific competition.

- inter_max:

  Maximum interspecific coefficient, relative to intraspecific
  competition.

- scale:

  Multiplier of all coefficients.

- traits:

  For `type = "trait"`: a data frame from
  [`species_traits()`](https://provetelab.org/comdynamite/reference/species_traits.md),
  or any data frame with columns `optimum` and `sigma`, or a list of
  them for a niche of several axes.

## Value

Square matrix: element \[j, k\] is the effect of species k on species j.

## Details

- `"none"`: no interspecific competition.

- `"random"`: interspecific coefficients drawn from U(0, `inter_max`).

- `"phylo"`: the same uniform values, assigned to pairs of species by
  rank of phylogenetic similarity (1 - d / max(d), with d the patristic
  distance), so that close relatives compete more strongly while the
  distribution of coefficients stays the same as with `"random"`.

- `"trait"`: competition falls with the distance between the niches of
  the two species, \\\alpha\_{jk} = q \exp(-(\theta_j - \theta_k)^2 /
  (2(\sigma_j^2 + \sigma_k^2)))\\, with \\q\\ = `inter_max`, \\\theta\\
  the optima and \\\sigma\\ the niche breadths of `traits`. The last
  factor is the overlap of two Gaussian niches, the kernel that
  [`simulate_evolving_metacommunity()`](https://provetelab.org/comdynamite/reference/simulate_evolving_metacommunity.md)
  uses with `niche_overlap = TRUE`, so that with `inter_max = 1` two
  species with the same niche compete as one species. With a list of
  traits (several axes), the overlap is the product over axes. Nothing
  is random.

## Examples

``` r
set.seed(1)
tr <- treeshapesim::sim_tree_shape(5)
round(competition_matrix(tr, "phylo", inter_max = 0.5), 3)
#>       t1    t2    t3    t4    t5
#> t1 0.050 0.023 0.003 0.016 0.012
#> t2 0.023 0.050 0.007 0.005 0.010
#> t3 0.003 0.007 0.050 0.019 0.018
#> t4 0.016 0.005 0.019 0.050 0.025
#> t5 0.012 0.010 0.018 0.025 0.050
traits <- species_traits(tr)
round(competition_matrix(type = "trait", traits = traits, inter_max = 1), 3)
#>       t1    t2    t3    t4    t5
#> t1 0.050 0.000 0.002 0.036 0.017
#> t2 0.000 0.050 0.004 0.000 0.000
#> t3 0.002 0.004 0.050 0.017 0.005
#> t4 0.036 0.000 0.017 0.050 0.049
#> t5 0.017 0.000 0.005 0.049 0.050
```
