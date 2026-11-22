# Deprecated functions

Renamed in 0.7.0 for consistency; the old names still work, with a
warning, until 1.0.0.

## Usage

``` r
traits_ou(...)

fit_damocles(...)
```

## Arguments

- ...:

  Passed to the new function.

## Value

The result of the new function.

## Details

- `traits_ou()` is now
  [`trait_ou()`](https://provetelab.org/comdynamite/reference/trait_ou.md),
  like
  [`trait_bm()`](https://provetelab.org/comdynamite/reference/trait_bm.md),
  [`trait_mk()`](https://provetelab.org/comdynamite/reference/trait_mk.md)
  and
  [`trait_rates()`](https://provetelab.org/comdynamite/reference/trait_rates.md).

- `fit_damocles()` is now
  [`damocles_fit()`](https://provetelab.org/comdynamite/reference/damocles_loglik.md),
  like the other functions of DAMOCLES.

## Examples

``` r
tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
set.seed(1)
old <- suppressWarnings(traits_ou(tree, alpha = 1))
set.seed(1)
identical(old, trait_ou(tree, alpha = 1))
#> [1] TRUE
```
