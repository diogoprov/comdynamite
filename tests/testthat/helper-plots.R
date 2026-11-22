# a plot of the package is a ggplot that can be built without error
expect_ggplot <- function(p) {
  testthat::expect_s3_class(p, "ggplot")
  ggplot2::ggplot_build(p)
  invisible(p)
}

# expect_warning() and expect_message() of the third edition return the condition, not the
# value of the code; these run the expectation and return the value, so that a test can keep
# it without an assignment inside the expectation
expect_warning_value <- function(object, ...) {
  env <- new.env()
  testthat::expect_warning(assign("value", object, envir = env), ...,
                           label = deparse1(substitute(object)))
  env$value
}
expect_message_value <- function(object, ...) {
  env <- new.env()
  testthat::expect_message(assign("value", object, envir = env), ...,
                           label = deparse1(substitute(object)))
  env$value
}
