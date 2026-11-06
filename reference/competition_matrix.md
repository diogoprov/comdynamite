# Competition matrix

Intraspecific coefficients are 1 and interspecific coefficients lie in
\[0, `inter_max`), all multiplied by `scale`. With the default `scale`
of 0.05, as in Thompson et al. (2020), a species alone at its optimum
with `r_max = 5` reaches an equilibrium of `(5 - 1) / 0.05 = 80`
individuals, and `inter_max` is the strength of interspecific relative
to intraspecific competition. The matrix is symmetric.

## Usage

``` r
competition_matrix(
  species,
  type = c("none", "random", "phylo"),
  inter_max = 0.1,
  scale = 0.05
)
```

## Arguments

- species:

  A `phylo` tree, or the number of species (not for `type = "phylo"`).

- type:

  Type of interspecific competition.

- inter_max:

  Maximum interspecific coefficient, relative to intraspecific
  competition.

- scale:

  Multiplier of all coefficients.

## Value

Square matrix: element \[j, k\] is the effect of species k on species j.

## Details

- `"none"`: no interspecific competition.

- `"random"`: interspecific coefficients drawn from U(0, `inter_max`).

- `"phylo"`: the same uniform values, assigned to pairs of species by
  rank of phylogenetic similarity (1 - d / max(d), with d the patristic
  distance), so that close relatives compete more strongly while the
  distribution of coefficients stays the same as with `"random"`.

## Examples

``` r
set.seed(1)
tr <- ape::rcoal(5)
round(competition_matrix(tr, "phylo", inter_max = 0.5), 3)
#>       t2    t5    t4    t3    t1
#> t2 0.050 0.023 0.007 0.003 0.009
#> t5 0.023 0.050 0.005 0.010 0.000
#> t4 0.007 0.005 0.050 0.022 0.016
#> t3 0.003 0.010 0.022 0.050 0.010
#> t1 0.009 0.000 0.016 0.010 0.050
```
