# Results as data frames, summaries and plots

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives the
main table of a result.
[`summary()`](https://rdrr.io/r/base/summary.html) gives the same table,
with the tests of
[`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md)
in order of `p_adjusted`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) calls
[`plot_node_test()`](https://provetelab.org/comdynamite/reference/plot_node_test.md)
or
[`plot_trait_partition()`](https://provetelab.org/comdynamite/reference/plot_trait_partition.md)
or
[`plot_phylo_beta()`](https://provetelab.org/comdynamite/reference/plot_phylo_beta.md).

## Usage

``` r
# S3 method for class 'mc_node_test'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'mc_node_test'
summary(object, ...)

# S3 method for class 'mc_node_test'
plot(x, tree, ...)

# S3 method for class 'mc_pcc'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'mc_pcc'
summary(object, ...)

# S3 method for class 'mc_phylo_beta'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'mc_phylo_beta'
summary(object, ...)

# S3 method for class 'mc_phylo_beta'
plot(x, comm, tree, ...)

# S3 method for class 'mc_phylo_vector'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'mc_phylo_vector'
summary(object, ...)

# S3 method for class 'mc_trait_partition'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'mc_trait_partition'
summary(object, ...)

# S3 method for class 'mc_trait_partition'
plot(x, ...)
```

## Arguments

- x, object:

  Result of
  [`node_test()`](https://provetelab.org/comdynamite/reference/node_test.md),
  [`pcc_test()`](https://provetelab.org/comdynamite/reference/pcc_test.md),
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md),
  [`phylo_vector()`](https://provetelab.org/comdynamite/reference/phylo_vector.md)
  or
  [`trait_partition()`](https://provetelab.org/comdynamite/reference/trait_partition.md).

- row.names, optional:

  Not used.

- ...:

  Passed to the plot function, and not used otherwise.

- tree:

  The `phylo` object of the analysis.

- comm:

  The community matrix given to
  [`phylo_beta_test()`](https://provetelab.org/comdynamite/reference/phylo_beta_test.md).

## Value

A data frame, or for
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) the value of
the plot function.

## Examples

``` r
nodes <- node_contrasts(daphniids$tree)
set.seed(1)
out <- node_test(daphniids$comm, nodes, daphniids$env, n_perm = 9)
head(as.data.frame(out), 3)
#>   node n_species n_sites         r2 r2_adjusted r2_weighted  statistic p_site
#> 1   14        13     239 0.18758306  0.15565401  0.16400769 0.21465241    0.2
#> 2   15        11     235 0.06581504  0.02844764  0.05076061 0.06735415    0.5
#> 3   16         9     233 0.39831212  0.37402875  0.33394192 0.51833271    0.1
#>   p_species p_max p_joint n_draws p_adjusted
#> 1       0.4   0.4     0.1       9       0.99
#> 2       1.0   1.0     0.3       9       1.00
#> 3       0.1   0.1     0.1       9       0.99
head(summary(out), 3)
#>   node n_species n_sites        r2 r2_adjusted r2_weighted statistic p_site
#> 1   16         9     233 0.3983121   0.3740288   0.3339419 0.5183327    0.1
#> 2   24         2      31 0.8146036   0.7351480   0.8146036 4.3938471    0.1
#> 3   18         5     198 0.2930576   0.2592146   0.2347089 0.3429577    0.1
#>   p_species p_max p_joint n_draws p_adjusted
#> 1       0.1   0.1     0.1       9       0.99
#> 2       0.2   0.2     0.1       9       0.99
#> 3       0.3   0.3     0.1       9       0.99
plot(out, daphniids$tree, test = "site")
```
