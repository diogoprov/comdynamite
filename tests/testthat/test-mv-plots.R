tr4 <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")

draw <- function(expr) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  })
  force(expr)
}

test_that("plot_traits() draws phylomorphospaces and traitgrams", {
  set.seed(1)
  reg <- paint_regimes(tr4, nodes = 7)
  x <- traits_ou(tr4, sigma2 = rate_matrix(c(1, 2), correlation = 0.5))
  expect_null(draw(plot_traits(x, tr4, regimes = reg, labels = TRUE)))
  expect_null(draw(plot_traits(x, tr4, traits = "trait2", labels = TRUE)))
  expect_null(draw(plot_traits(x, labels = TRUE)))
  expect_null(draw(plot_traits(x, tr4, regimes = reg, colours = c(`1` = "red", `0` = "blue"))))
  h <- tr4
  h$maps <- lapply(h$edge.length, function(len) c(a = len))
  h$maps[[1]] <- c(a = 0.5, b = 0.5)
  class(h) <- c("simmap", "phylo")
  expect_null(draw(plot_traits(x, tr4, regimes = h)))
  expect_error(plot_traits(x, tr4, traits = 1:3), "one or two")
  expect_error(plot_traits(x, tr4, traits = "z"), "columns of `x`")
  expect_error(plot_traits(x[, 1, drop = FALSE], tr4), "traits = 1")
  expect_error(plot_traits(unname(x), tr4), "row names")
  expect_error(plot_traits(x, traits = 1), "traitgram needs")
  expect_error(plot_traits(x, ape::read.tree(text = "((a:1,b:1):1,c:2);")), "does not match")
  expect_error(plot_traits(x, tr4, regimes = c("0", "1")), "one value per branch")
  expect_error(plot_traits(x, tr4, regimes = tr4), "simmap")
  expect_error(plot_traits(x, tr4, regimes = reg, colours = c(`0` = "red")), "every regime")
  expect_error(plot_traits(x, tr4, regimes = reg, colours = "red"), "2 colours")
})

test_that("ellipse_points() traces the contour of the quadratic form", {
  s <- matrix(c(2, 0.8, 0.8, 1), 2)
  p <- ellipse_points(s, 1.5)
  expect_equal(rowSums((p %*% solve(s)) * p), rep(1.5^2, nrow(p)))
})

test_that("plot.rate_matrix() draws one ellipse per regime", {
  r <- rate_matrix(c(1, 2), correlation = list(a = 0.5, b = -0.5), constraint = "variance")
  expect_null(draw(plot(r)))
  expect_null(draw(plot(rate_matrix(c(1, 2, 3)), traits = c(1, 3), level = 2)))
  expect_error(plot(rate_matrix(1)), "no ellipse")
  many <- rate_matrix(c(1, 2), constraint = "proportional",
                      scale = stats::setNames(1:8, letters[1:8]))
  expect_null(draw(plot(many)))
})

test_that("plot_assembly() draws each rule", {
  set.seed(3)
  x <- matrix(rnorm(60), 30, dimnames = list(paste0("s", 1:30), NULL))
  for (rule in c("filtering", "limiting_similarity", "neutral")) {
    a <- assemble_community(x, 10, rule, optimum = c(0, 0), distance = "mahalanobis")
    expect_null(draw(plot_assembly(a, labels = TRUE)))
  }
  expect_error(plot_assembly(list()), "output of assemble_community")
  one <- assemble_community(stats::setNames(rnorm(5), letters[1:5]), 2)
  expect_error(plot_assembly(one), "one trait")
})
