# Contributing to comdynamite

Thank you for considering a contribution. This file explains how to report a problem, propose a
change and submit code.

## Development stage

comdynamite is in active development and has not reached version 1.0.0. The analysis and
simulation functions are covered by automated tests, but their interfaces can still change
between minor versions. When a function or argument is renamed, the old name
keeps working with a deprecation warning (through `.Deprecated()`) until version 1.0.0, and the
change is listed in `NEWS.md`. Each release is archived on Zenodo with its own DOI.

## Reporting a problem

Open an [issue](https://github.com/diogoprov/comdynamite/issues) with a minimal reproducible
example, ideally made with the [reprex](https://reprex.tidyverse.org) package, and the output of
`sessionInfo()`. For simulations, include the seed (`set.seed()`) so the run can be repeated.

## Proposing a change

For anything beyond a small fix (a typo, a broken link, a clearer sentence in the documentation),
open an issue first to discuss the change. A model or method that is new to the package should
come with the reference that defines it, and with a description of how its implementation can be
checked against published results.

## Submitting code

1. Fork the repository and create a branch from `main`.
2. Follow the existing style: base R where possible, lines of at most 100 characters (see
   `.lintr`), and roxygen2 for documentation.
3. Add tests with testthat for every new behaviour, including the errors that invalid inputs
   should raise. Run `devtools::test()` and `devtools::check()` before opening the pull request.
4. Add a bullet to the top section of `NEWS.md` describing the change.
5. Open a pull request against `main`. Continuous integration runs `R CMD check` on Linux, macOS
   and Windows and reports test coverage.

Code from other packages or publications can only be included when its license is compatible
with the MIT license of comdynamite. Code under GPL can be used only through `Suggests`.

## Code of conduct

Please note that this project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree to
abide by its terms.
