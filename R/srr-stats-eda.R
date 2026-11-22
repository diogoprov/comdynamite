#' srr_stats: Exploratory data analysis standards
#'
#' rOpenSci statistical standards (`srr`) of this category, with the parts of the
#' package that comply with each one or the reason it does not apply.
#'
#' @srrstatsTODO {EA1.0} The documentation states what the package does and the questions
#'   it answers, but no file names the target audience. Missing: one explicit sentence on
#'   the intended users in `README.md` and in the package help (`R/comdynamite-package.R`).
#' @srrstats {EA1.1} The kinds of data are described in the article "Preparing your data"
#'   (`vignettes/articles/preparing-data.Rmd`) and in the help of `mc_data()`
#'   (`R/mc_data.R`): a numeric community matrix of sites by species, a `phylo` tree,
#'   tables of species traits, numeric site predictors, site coordinates, areas or
#'   biogeographic events, and an optional array of sites by species by steps. The
#'   functions for trait variation among individuals (`tstats()`, `community_index()`,
#'   `trait_flex()`, `nested_variance()`, `trait_densities()`) take a long data frame with
#'   one row per individual.
#' @srrstats {EA1.2} The package help (`?comdynamite`, sections "Analyse data" and
#'   "Simulate data" in `R/comdynamite-package.R`) and the table "What do you want to do?"
#'   of `README.md` list the questions the software helps explore, each with its functions
#'   and, in the README, its article.
#' @srrstats {EA1.3} Every `@param` entry states the type and layout of the input it
#'   accepts, for example `comm`, `tree`, `traits`, `env`, `coords`, `areas` and `series`
#'   in `mc_data()` (`R/mc_data.R`), and the `individuals` data frame with its `site`,
#'   `species` and `traits` column arguments in `tstats()` (`R/individuals.R`) and
#'   `community_index()` (`R/trait_variation.R`).
#' @srrstats {EA2.0} The tables of an analysis are joined by the names of species and of
#'   sites, which act as the index: species by the column names of `comm`, the tip labels
#'   of `tree` and the row names of `traits` (`match_species()` in `R/phylo_beta.R`), and
#'   sites by row names (`check_sites()` in `R/checks.R`, `site_rows()` in `R/mc_data.R`).
#'   `mc_data()` does this matching once and keeps the parts aligned, including under
#'   subsetting with `[`. Tested in `tests/testthat/test-mc-data.R` and
#'   `tests/testthat/test-species-match.R`.
#' @srrstats {EA2.1} Uniqueness is affirmed for species names (`match_species()` stops with
#'   "repeated species", tested in `tests/testthat/test-species-match.R`) and for site names:
#'   `check_unique_sites()` (`R/checks.R`) stops on repeated site names in `check_sites()`, in
#'   `site_rows()` (`R/mc_data.R`) for the row names of `env`, `coords` and `areas`, and in
#'   `trait_flex()` for the site column of `env` (tests in `tests/testthat/test-checks.R`).
#' @srrstats {EA2.2} The index is identified through the class system of EA2.2a.
#' @srrstats {EA2.2a} The class `mc_data` (`R/mc_data.R`) fixes the index of each part:
#'   the row names of `comm`, `env`, `coords` and `areas` are the sites and the column
#'   names of `comm`, the tip labels of `tree` and the row names of `traits` are the
#'   species. Its constructor sets these names and its `[` method keeps them.
#' @srrstats {EA2.3} No join assumes a column name. Community, trait and predictor tables
#'   are matched by their row and column names (`match_species()`, `check_sites()`,
#'   `site_rows()`). In the long tables of individuals, the columns that identify sites,
#'   species, traits and abundances are given by the arguments `site`, `species`, `traits`
#'   (or `trait`) and `abundance` of `tstats()`, `community_index()`, `trait_flex()`,
#'   `trait_densities()` and `subsample_individuals()`, and `check_data_columns()` stops
#'   when a named column is absent. In `nested_variance()` the columns are given by
#'   `traits` and `factors`.
#' @srrstatsTODO {EA2.4} The class `mc_data` holds the multi-table input and the analyses
#'   accept it as their first argument (`node_test()`, `pcc_test()`, `pcps_test()`,
#'   `phylo_beta()`, `phylo_beta_test()`, `phylo_vector()`, `trait_partition()`,
#'   `summary_stats()`; tested in `tests/testthat/test-mc-data.R`). Missing: the class is
#'   not demanded. The same functions also accept `comm`, `tree`, `traits` and `env` as
#'   separate objects, which are matched by name inside each function. Either coerce the
#'   separate objects with `mc_data()` at the start of each analysis or document this as a
#'   justified deviation.
#' @srrstats {EA2.5} The tables inside `mc_data` are keyed by site and species names, whose
#'   uniqueness is checked as in EA2.1, so each site and species of the index matches one row of
#'   every part.
#' @srrstats {EA2.6} Extra attributes of input vectors are ignored. `check_comm()` and
#'   `check_env()` (`R/checks.R`) convert inputs with `as.matrix()`, which keeps only the
#'   values and names, and the functions for individuals use the trait columns only in
#'   numeric operations. Checked by hand on vectors with `units` and `label` attributes
#'   and an extra class, for which `pcps_test()`, `tstats()` and `trait_flex()` gave the
#'   same results as for plain vectors. There is no test of this in `tests/testthat/`.
#' @srrstats {EA3.0} Several steps are returned as tables by one call: `node_test()` tests
#'   every node of the tree and adjusts the p-values (`p_adjusted`); `tstats_summary()`
#'   gives, per trait and statistic, the share of sites that deviate from the null model
#'   and the sites concerned; `null_intervals()` gives the null intervals, standardised
#'   effect sizes and direction of deviation of `tstats()` and `community_index()`;
#'   `abc_choose()` chooses among simulation models and reports the out-of-bag error and
#'   the confusion matrix; and `run_scenarios()` runs a grid of scenarios and replicates
#'   and returns one table of statistics.
#' @srrstats {EA3.1} `summary_stats()` (`R/summary_stats.R`) returns the same named vector
#'   of statistics for simulated and observed data, with `NA` for statistics that lack
#'   their inputs, so results can be bound by rows and compared. `run_scenarios()`
#'   (`R/workflow.R`) compares scenarios on the same tree and landscape within each
#'   replicate, with one seed per cell. `trait_gradient_weights()` puts the results of the
#'   combinations of weights side by side. The `as.data.frame()` methods of the whole-tree
#'   tests (`R/methods.R`) share the columns `statistic`, `p_site`, `p_species` and
#'   `p_max`.
#' @srrstatsTODO {EA4.0} Partly met: `simulate_metacommunity()` returns integer abundances
#'   (tested in `tests/testthat/test-dynamics.R`), statistics that are means or ratios are
#'   double, and `trait_range()` of integers is integer. Missing: `check_comm()`
#'   (`R/checks.R`) multiplies by 1, so `mc_data()` and its `[` method return a double
#'   `comm` for an integer input; and `tstats()`, `community_index()` and
#'   `trait_densities()` return `site` as character when the input column is a factor.
#'   Either keep the input type in these cases or document the conversion in the return
#'   values.
#' @srrstatsTODO {EA4.1} No function has a parameter for numeric precision. The `print`
#'   methods fix it in the code (`format(., digits = 3)`, `signif(., 3)`, `sprintf()`).
#'   Missing: a `digits` argument in the `print` methods, passed to these calls.
#' @srrstats {EA4.2} Results with their own class have `print` methods (27 in the
#'   package, for example `mc_node_test`, `mc_pcc`, `mc_pcps`, `mc_phylo_beta`,
#'   `mc_trait_gradient`, `mc_trait_partition`, `mc_data`, `mc_sim`, `drift_test`), and
#'   the results of `tstats()`, `community_index()`, `nested_variance()` and
#'   `trait_rarefaction()` are data frames with an extra class, printed as data frames.
#'   `plot` methods exist for `mc_node_test`, `mc_pcps`, `mc_phylo_beta`,
#'   `mc_trait_gradient`, `mc_trait_gradient_test`, `mc_trait_partition`, `tstats`,
#'   `community_index`, `trait_flex`, `nested_variance`, `trait_rarefaction`,
#'   `population_shifts`, `rate_matrix`, `drift_test` and `food_web`, and `summary`
#'   methods are documented in `?"result-methods"` (`R/methods.R`). The result of
#'   `trait_flex()` has no `print` method and prints as a plain list of two data frames.
#'   The numeric vector of `summary_stats()` has no plot method.
#' @srrstats {EA5.0} All plots are ggplot objects (`R/plots.R`, `R/trait_plots.R`,
#'   `R/mv_plots.R`, `R/plot_rates.R`, `R/food_web_plot.R` and the `plot` methods), so
#'   text sizes, colours and themes can be changed by adding ggplot2 layers, as the help
#'   of `plot.tstats()` and `plot_population_traits()` shows. Colours follow EA5.0b; text
#'   sizes are discussed in EA5.0a.
#' @srrstatsTODO {EA5.0a} Axis and legend text use the defaults of
#'   `ggplot2::theme_minimal()` and are not enlarged. Text drawn inside the plots defaults
#'   to 2.5 to 3.2 mm (about 7 to 9 pt), set by `text_size` in `plot_node_test()`,
#'   `plot_node_heatmap()`, `plot_trait_partition()` and `plot_phylo_beta()`, and fixed at
#'   2.5 mm for the labels of `plot_traits()` and `plot_assembly()` (`R/mv_plots.R`).
#'   Missing: larger defaults for these labels and a `text_size` argument in
#'   `plot_traits()` and `plot_assembly()`.
#' @srrstats {EA5.0b} Default colours come from the Okabe-Ito palette, which is safe for
#'   colour-blind readers: blue and vermillion for two groups (`plot_node_test()`,
#'   `plot_node_heatmap()`), blue, grey and vermillion in `plot.tstats()` and
#'   `plot_assembly()`, and up to seven colours for species or regimes in
#'   `plot_population_traits()` and `plot_traits()`, with
#'   `grDevices::hcl.colors(n, "Dark 3")` for more levels. Continuous scales use
#'   `grDevices::hcl.colors()` with `"viridis"` as the default palette (`plot_rates()`).
#'   Every plot has a `colours` or `palette` argument.
#' @srrstatsTODO {EA5.2} The `print` methods of the result classes format numbers with
#'   `format(., digits = 3)`, `signif()` or `sprintf()` (for example `print.mc_pcps()`,
#'   `print.mc_phylo_beta()`, `print.mc_trait_gradient()`, `print.drift_test()`). Missing:
#'   the number of digits is fixed in the code and not set by a parameter (EA4.1), and the
#'   results that print as data frames or lists (`tstats()`, `community_index()`,
#'   `null_intervals()`, `nested_variance()`, `trait_flex()`, the `summary` methods of
#'   `R/methods.R`) and the vector of `summary_stats()` rely on the default formatting of
#'   numbers.
#' @srrstats {EA5.4} Axis breaks are those of the ggplot2 scales, which choose round
#'   values. Numbers written inside the plots are rounded: two decimals in
#'   `plot_node_heatmap()` and `plot_trait_partition()` (`R/plots.R`), whole percentages
#'   in `plot.nested_variance()` (`R/trait_variation.R`), two significant digits in the
#'   legend of `plot_rates()` and two decimals in the subtitle of `plot.drift_test()`.
#' @srrstats {EA5.5} The inputs are plain numeric objects without units, so no unit can be
#'   read from the data. Axes of traits and predictors are titled with the column names
#'   of the input (`plot_traits()`, `plot_assembly()`, `plot_population_traits()`), and
#'   the other axes with what they measure, such as "Time from the root" and
#'   "Standardised effect size". Units can be added to the returned ggplot object with
#'   `ggplot2::labs()`.
#' @srrstats {EA6.0} Every exported function is called in `tests/testthat/`, and the
#'   tests check the properties in EA6.0a to EA6.0c and EA6.0e. Column types (EA6.0d) are
#'   not yet tested.
#' @srrstats {EA6.0a} Classes and types are tested with `expect_s3_class()` (45 test
#'   files) and `expect_type()`, for example the classes of the results in
#'   `tests/testthat/test-mc-data.R`, `test-pcps.R`, `test-trait-gradient.R` and
#'   `test-trait_variation.R`, and the integer abundances of the simulation in
#'   `test-dynamics.R`.
#' @srrstats {EA6.0b} Dimensions are tested with `expect_equal(dim(.), .)`,
#'   `expect_equal(nrow(.), .)` and `expect_length()` (55 test files), for example the
#'   subsets of `mc_data` in `tests/testthat/test-mc-data.R` and the tables and null
#'   arrays of `community_index()` and `tstats_summary()` in `test-trait_variation.R`.
#' @srrstats {EA6.0c} Column and element names are tested with `expect_named()` and
#'   `expect_equal(names(.), .)` (31 test files), for example the long table of
#'   `as.data.frame.mc_data()` in `tests/testthat/test-mc-data.R`, the columns of
#'   `community_index()` in `test-trait_variation.R` and the names of `summary_stats()`
#'   in `test-summary-stats.R`.
#' @srrstatsTODO {EA6.0d} Missing: the tests check the names and dimensions of the data
#'   frames returned but not the class of each column. Add expectations such as
#'   `expect_identical(vapply(out, class, character(1)), c(site = "character", ...))` for
#'   the tables of `tstats()`, `community_index()`, `null_intervals()`, `tstats_summary()`,
#'   `trait_flex()`, `nested_variance()`, `trait_densities()`, `trait_rarefaction()` and
#'   the `as.data.frame()` methods.
#' @srrstats {EA6.0e} Numeric results are compared with analytic, reference or recomputed
#'   values through `expect_equal()` with an explicit `tolerance` (40 test files), for
#'   example the equilibria in `tests/testthat/test-analytic.R` and the sums of squares of
#'   `trait_flex()` in `test-trait_variation.R`. The extended tests state in comments how
#'   each tolerance follows from theory or from Monte Carlo error (`tests/README.md`).
#' @srrstats {EA6.1} Graphical output is tested without vdiffr, on the ggplot objects.
#'   `expect_ggplot()` (`tests/testthat/helper-plots.R`) checks the class and builds each
#'   plot with `ggplot2::ggplot_build()`, in `test-plots.R`, `test-mv-plots.R`,
#'   `test-food-web-plot.R` and `test-trait-gradient.R`. The tests also check the data,
#'   number and type of the layers with `ggplot2::layer_data()` in `test-plots.R`,
#'   `test-trait_plots.R`, `test-trait_variation.R`, `test-sampling.R`,
#'   `test-population_shifts.R` and `test-plot-rates.R`.
#' @noRd
NULL

#' NA_standards
#'
#' Standards of this category that do not apply to this package, with the reasons.
#'
#' @srrstatsNA {EA2.2b} The index is identified by the class `mc_data` (EA2.2a), the
#'   alternative to an `index` attribute.
#' @srrstatsNA {EA5.1} The package never sets a typeface: no plot function passes a font
#'   family to ggplot2, so the typeface is the default of ggplot2 and of the graphics
#'   device.
#' @srrstatsNA {EA5.3} No function prints column-wise summaries of a table with columns of
#'   mixed types. Community, predictor and trait inputs must be numeric
#'   (`check_numeric_columns()` and `check_comm()` in `R/checks.R`), and `summary_stats()`
#'   returns statistics of the whole metacommunity, not of columns.
#' @srrstatsNA {EA5.6} The package bundles no library for dynamic visualisation. `inst/`
#'   holds only `CITATION`, `WORDLIST` and one Stan model, and `animate_food_web()` writes
#'   a GIF with the Suggested package gifski.
#' @noRd
NULL
