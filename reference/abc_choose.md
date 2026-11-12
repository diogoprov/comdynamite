# Choose among models by ABC with random forests

Fits a classification random forest that predicts which model produced a
simulation from its statistics (Pudlo et al. 2016), with the linear
discriminant axes of the statistics added as predictors. The out-of-bag
error and the confusion matrix measure, before looking at the data, how
often the models can be told apart.
[`predict()`](https://rdrr.io/r/stats/predict.html) gives the chosen
model for observed statistics, the votes of the trees and the posterior
probability of the chosen model, which is more reliable than the
proportion of votes.

## Usage

``` r
abc_choose(tables, statistics = NULL, ntree = 500, lda = TRUE, ncores = 1)

# S3 method for class 'mc_abc_choice'
predict(object, newdata, ntree = 500, ncores = 1, ...)
```

## Arguments

- tables:

  Named list of reference tables, one per model, as from
  [`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
  or a data frame with a column `model`. Only the statistics shared by
  all models are used.

- statistics:

  Names of the statistics. By default all statistics shared by the
  models.

- ntree:

  Number of trees.

- lda:

  Add the linear discriminant axes as predictors?

- ncores:

  Number of threads.

- object:

  An `mc_abc_choice`.

- newdata:

  Named numeric vector of observed statistics or a data frame with one
  observation per row.

- ...:

  Not used.

## Value

An `mc_abc_choice`, a list with `forest`, `training`, `statistics`,
`models`, `error` (out-of-bag error) and `confusion` (rows: true model;
columns: chosen model).

For [`predict()`](https://rdrr.io/r/stats/predict.html), a data frame
with `observation`, `model`, `probability` and one column of votes per
model.

## Functions

- `predict(mc_abc_choice)`: Chosen model, votes and posterior
  probability of the chosen model.

## References

Pudlo, P., Marin, J.-M., Estoup, A., Cornuet, J.-M., Gautier, M., &
Robert, C. P. (2016). Reliable ABC model choice via random forests.
*Bioinformatics*, 32, 859–866.
[doi:10.1093/bioinformatics/btv684](https://doi.org/10.1093/bioinformatics/btv684)

## See also

Other inference:
[`abc_cv()`](https://provetelab.org/comdynamite/reference/abc_cv.md),
[`abc_diagnose()`](https://provetelab.org/comdynamite/reference/abc_diagnose.md),
[`abc_fit()`](https://provetelab.org/comdynamite/reference/abc_fit.md),
[`abc_ppc()`](https://provetelab.org/comdynamite/reference/abc_ppc.md),
[`abc_recovery()`](https://provetelab.org/comdynamite/reference/abc_recovery.md),
[`abc_table()`](https://provetelab.org/comdynamite/reference/abc_table.md),
[`fold_summarise()`](https://provetelab.org/comdynamite/reference/fold_summarise.md),
[`observe_detection()`](https://provetelab.org/comdynamite/reference/observe_detection.md),
[`spatial_folds()`](https://provetelab.org/comdynamite/reference/spatial_folds.md)

## Examples

``` r
set.seed(1)
a <- data.frame(s1 = rnorm(200, 0), s2 = rnorm(200))
b <- data.frame(s1 = rnorm(200, 1.5), s2 = rnorm(200))
choice <- abc_choose(list(neutral = a, niche = b), ntree = 100)
choice$confusion
#>         neutral niche
#> neutral     140    60
#> niche        53   147
predict(choice, c(s1 = 1.2, s2 = 0))
#>   observation model probability votes_neutral votes_niche
#> 1           1 niche      0.8764            70         430
```
