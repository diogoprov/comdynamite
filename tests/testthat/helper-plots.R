# a plot of the package is a ggplot that can be built without error
expect_ggplot <- function(p) {
  testthat::expect_s3_class(p, "ggplot")
  ggplot2::ggplot_build(p)
  invisible(p)
}
