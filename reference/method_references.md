# Original references of the methods in comdynamite

Many functions of comdynamite implement methods published by other
authors. The first time one of these functions is used in a session, it
shows a message with the original papers of its method, which should be
cited together with the package (`citation("comdynamite")`).
`method_references()` gives the same references at any time.

## Usage

``` r
method_references(fun = NULL)
```

## Arguments

- fun:

  Name of a function of comdynamite, such as `"node_test"`, or `NULL`
  for all functions that implement a published method.

## Value

A character vector of references for one function. With `fun = NULL`, a
data frame with columns `fun` and `reference`, one row per function and
paper.

## Details

The message is shown once per function and session. Set
`options(comdynamite.cite = FALSE)` to turn it off, for example in a
script that runs many sessions;
[`suppressMessages()`](https://rdrr.io/r/base/message.html) also hides
it. A function called by another one (such as
[`node_contrasts()`](https://provetelab.org/comdynamite/reference/node_contrasts.md)
inside
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md))
does not show its own message.

## See also

[`citation()`](https://rdrr.io/r/utils/citation.html)

## Examples

``` r
method_references("node_test")
#> [1] "Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity phylogenetics: separating the roles of environmental filters and historical biogeography. Ecology Letters, 13, 1290-1299. https://doi.org/10.1111/j.1461-0248.2010.01523.x"
#> [2] "Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects of spatial contingency and environmental filtering on metacommunity phylogenetics. Ecology, 93, S14-S30. https://doi.org/10.1890/11-0494.1"                           
head(method_references())
#>               fun
#> 1  node_contrasts
#> 2       node_test
#> 3       node_test
#> 4    phylo_vector
#> 5 trait_partition
#> 6 trait_partition
#>                                                                                                                                                                                                                                           reference
#> 1 Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity phylogenetics: separating the roles of environmental filters and historical biogeography. Ecology Letters, 13, 1290-1299. https://doi.org/10.1111/j.1461-0248.2010.01523.x
#> 2 Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity phylogenetics: separating the roles of environmental filters and historical biogeography. Ecology Letters, 13, 1290-1299. https://doi.org/10.1111/j.1461-0248.2010.01523.x
#> 3                            Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects of spatial contingency and environmental filtering on metacommunity phylogenetics. Ecology, 93, S14-S30. https://doi.org/10.1890/11-0494.1
#> 4                            Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects of spatial contingency and environmental filtering on metacommunity phylogenetics. Ecology, 93, S14-S30. https://doi.org/10.1890/11-0494.1
#> 5                            Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012). Assessing the effects of spatial contingency and environmental filtering on metacommunity phylogenetics. Ecology, 93, S14-S30. https://doi.org/10.1890/11-0494.1
#> 6  Peres-Neto, P. R., Legendre, P., Dray, S., & Borcard, D. (2006). Variation partitioning of species data matrices: estimation and comparison of fractions. Ecology, 87, 2614-2625. https://doi.org/10.1890/0012-9658(2006)87[2614:VPOSDM]2.0.CO;2
```
