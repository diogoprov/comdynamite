# Tests of comdynamite

The tests use testthat (edition 3) and run with `devtools::test()` or `R CMD check`.

## Regular tests

Every file `tests/testthat/test-*.R` except `test-extended-*.R` runs on every check, on CRAN
and in continuous integration. They cover the expected output of each function on small data,
the errors that invalid inputs raise (`test-validation-*.R`), and quick checks that the results
do not depend on trivial noise in the inputs or on the random seed beyond Monte Carlo error
(`test-noise-*.R`). Tests that need Suggested packages (treeshapesim, cmdstanr, ggplot2 and
others) are skipped when those packages are missing, and slow ones are skipped on CRAN.

## Extended tests

The files `tests/testthat/test-extended-*.R` check statistical properties over many simulated
data sets, with the parameters known:

- type I error and power of the permutation and null-model tests (`test-extended-perm.R`);
- recovery of trait gradients, variance components and the parameters of the Bayesian model
  (`test-extended-traits.R`);
- recovery of the parameters of the evolutionary and diversification models and error rates of
  their tests (`test-extended-evol.R`);
- agreement of the simulators with analytic results and recovery by ABC (`test-extended-sims.R`).

Each test states, in a comment, the rOpenSci statistical software standards it covers and how
its tolerance follows from theory or from Monte Carlo error. They need no download. They take
several minutes, so they are skipped unless the environment variable
`COMDYNAMITE_EXTENDED_TESTS` is `"true"`:

```r
Sys.setenv(COMDYNAMITE_EXTENDED_TESTS = "true")
devtools::test(filter = "extended")
```

or, from a shell,

```sh
COMDYNAMITE_EXTENDED_TESTS=true Rscript -e 'devtools::test(filter = "extended")'
```

The Bayesian tests also need cmdstanr and CmdStan (`cmdstanr::install_cmdstan()`).
The workflow `.github/workflows/extended-tests.yaml` runs the extended tests on GitHub Actions
once a month and can be started by hand from the Actions tab.
