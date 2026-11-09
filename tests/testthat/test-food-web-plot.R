fw <- function() {
  food_web_dynamics(seq(0, 600, by = 20), init = c(R = 50, N = 0.1, P = 0.1),
                    traits = c(R = 10, N = 4, P = 6), introductions = c(P = 300), dt = 0.5)
}

test_that("plot.food_web draws the panels and frames", {
  x <- fw()
  skip_if_not_installed("ggplot2")
  expect_ggplot(plot(x))
  frame <- expect_ggplot(plot(x, panels = "trait", time = 200, log_abundance = TRUE))
  expect_true(all(frame$data$time <= 200))
  expect_ggplot(plot(x, panels = c("abundance", "rate"), log_abundance = TRUE,
                     legend = FALSE))
  expect_warning(plot(x, col = c("red", "blue", "green")), "`col` is deprecated")
  # a level absent during the whole run is left out
  y <- food_web_dynamics(seq(0, 100, by = 20), dt = 0.5)
  expect_ggplot(plot(y))
  expect_error(plot(x, time = c(1, 2)), "one number")
  expect_error(plot(x, colours = "red"), "three colours")
  expect_error(plot(x, panels = "abundances"), "should be one of")
  z <- x
  z[c("R", "N", "P")] <- 0
  expect_error(plot(z), "every level is absent")
  expect_error(plot.food_web(as.data.frame(x)), "result of `food_web_dynamics")
})

test_that("animate_food_web writes a GIF", {
  skip_if_not_installed("gifski")
  x <- fw()
  file <- tempfile(fileext = ".gif")
  expect_identical(animate_food_web(x, file, frames = 4, width = 200, height = 200), file)
  expect_true(file.size(file) > 0)
  expect_error(animate_food_web(x, file, frames = 1000), "number of rows")
  expect_error(animate_food_web(x, file, frames = 4, delay = 0), "positive")
})

test_that("animate_food_web asks for gifski when it is missing", {
  local_mocked_bindings(has_gifski = function() FALSE)
  expect_error(animate_food_web(fw(), tempfile()), "gifski")
})
