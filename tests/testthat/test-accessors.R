test_that("drift_test() accessors agree with the stored regression", {
  set.seed(1)
  G <- matrix(c(4, 1, 0.5, 0, 1, 3, 0.2, 0, 0.5, 0.2, 2, 0.3, 0, 0, 0.3, 1), 4)
  means <- matrix(rnorm(25 * 4), 25) %*% chol(G * 0.5)
  grp <- rep(seq_len(25), each = 10)
  ind <- means[grp, ] + matrix(rnorm(250 * 4), 250) %*% chol(G)
  res <- drift_test(ind, grp, correct = FALSE)
  expect_equal(coef(res), c(intercept = res$intercept, slope = res$slope))
  ci <- confint(res, level = res$level)
  expect_equal(unname(ci["slope", ]), res$conf_int)
  expect_identical(rownames(confint(res, "slope")), "slope")
  expect_identical(confint(res, 2), confint(res, "slope"))
  expect_error(confint(res, "beta"), "`parm`")
  expect_equal(sqrt(vcov(res)["slope", "slope"]), res$se)
  expect_identical(nobs(res), 4L)
  expect_equal(fitted(res) + residuals(res), log(res$components$between))
  expect_length(fitted(res), 4)
})

test_that("trait_flex() keeps its model and prints a summary", {
  set.seed(1)
  ind <- data.frame(site = rep(paste0("s", 1:6), each = 20),
                    species = sample(c("a", "b", "c"), 120, replace = TRUE))
  env <- data.frame(site = paste0("s", 1:6), moisture = 1:6)
  ind$trait <- c(a = 10, b = 20, c = 30)[ind$species] +
    env$moisture[match(ind$site, env$site)] + rnorm(120)
  res <- trait_flex(ind, formula = ~moisture, env = env)
  expect_identical(formula(res), ~moisture)
  m <- res$model
  expect_identical(names(m), c("site", "trait", "component", "response", "fitted", "residual"))
  expect_equal(fitted(res) + residuals(res), m$response)
  total <- m[m$component == "total", ]
  expect_equal(sum(total$residual^2),
               res$decomposition$residual[res$decomposition$component == "total"])
  expect_identical(nobs(res), 6L)
  expect_output(print(res), "6 sites, 1 trait")
  expect_output(print(res, digits = 2), "Model: ~moisture")
  expect_error(print(res, digits = 0), "`digits`")
  ind$height <- ind$trait + rnorm(120)
  expect_output(print(trait_flex(ind, traits = c("trait", "height"), formula = ~moisture,
                                 env = env)), "6 sites, 2 traits")
  plain <- trait_flex(ind)
  expect_equal(residuals(plain)[plain$model$component == "total"],
               total$response - mean(total$response))
})

test_that("print methods take the number of digits", {
  set.seed(1)
  G <- diag(c(4, 3, 2, 1))
  grp <- rep(seq_len(25), each = 10)
  ind <- (matrix(rnorm(25 * 4), 25) %*% chol(G * 0.5))[grp, ] +
    matrix(rnorm(250 * 4), 250) %*% chol(G)
  res <- drift_test(ind, grp, correct = FALSE)
  long <- capture.output(print(res, digits = 6))
  short <- capture.output(print(res, digits = 1))
  expect_gt(sum(nchar(long)), sum(nchar(short)))
  expect_match(short[3], sprintf("%.1f", res$slope), fixed = TRUE)
  expect_error(print(res, digits = 1.5), "`digits`")
})
