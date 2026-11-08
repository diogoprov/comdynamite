test_that("food_web_params checks its values and keeps the names of the book", {
  p <- food_web_params()
  expect_s3_class(p, "food_web_params")
  expect_equal(p$c0, 1)
  expect_equal(food_web_params(c0 = 10)$c0, 10)
  expect_error(food_web_params(c0 = -1), "zero or more")
  expect_error(food_web_params(beta = 0), "positive")
  expect_error(food_web_params(d = "a"), "one finite number")
  expect_error(food_web_params(h = c(1, 2)), "one finite number")
  # optima and steepnesses may be negative
  expect_equal(food_web_params(zc_R = -3, alpha = -0.1)$alpha, -0.1)
})

test_that("food_web_dynamics returns one row per time with the expected columns", {
  x <- food_web_dynamics(c(0, 10, 20), dt = 0.5)
  expect_s3_class(x, c("food_web", "data.frame"))
  expect_equal(nrow(x), 3)
  expect_equal(x$time, c(0, 10, 20))
  expect_named(x, c("time", "R", "N", "P", "z_R", "z_N", "z_P", "grad_R", "grad_N", "grad_P",
                    "rate_R", "rate_N", "rate_P", "a_RN", "m_NP", "v_RP"))
  expect_equal(x$rate_N, abs(0.2 * x$grad_N))
  expect_equal(attr(x, "attack"), "gaussian")
  expect_s3_class(attr(x, "params"), "food_web_params")
  expect_output(print(x), "<food_web>")
  # unnamed vectors are read in the order R, N, P
  y <- food_web_dynamics(c(0, 10, 20), init = c(50, 0.1, 0), traits = c(10, 4, 1), dt = 0.5)
  expect_equal(y$z_N, x$z_N)
  z <- food_web_dynamics(c(0, 10, 20), init = c(N = 0.1, P = 0, R = 50), dt = 0.5)
  expect_equal(z$R, x$R)
})

test_that("an equilibrium stays put", {
  # the resource alone at its optimum and R* = c0 / d
  alone <- food_web_dynamics(c(0, 100, 500), init = c(R = 50, N = 0, P = 0),
                             traits = c(R = 10, N = 4, P = 1), dt = 0.5)
  expect_equal(alone$R, rep(50, 3))
  expect_equal(alone$z_R, rep(10, 3))
  expect_equal(alone$rate_R, rep(0, 3))
  # the consumer-resource equilibrium of figure 3.2, reached by a long run
  long <- food_web_dynamics(c(0, 8000), init = c(R = 50, N = 0.01, P = 0),
                            traits = c(R = 10, N = 4, P = 1), dt = 0.5)
  end <- long[2, ]
  again <- food_web_dynamics(c(0, 200), init = c(R = end$R, N = end$N, P = 0),
                             traits = c(R = end$z_R, N = end$z_N, P = 1), dt = 0.5)
  expect_equal(again$R[2], end$R, tolerance = 1e-3)
  expect_equal(again$z_R[2], end$z_R, tolerance = 1e-3)
  expect_equal(again$z_N[2], end$z_N, tolerance = 1e-3)
  expect_lt(max(again$rate_R, again$rate_N), 1e-4)
})

test_that("absent levels stay at zero with their traits held", {
  x <- food_web_dynamics(seq(0, 300, by = 100), init = c(R = 50, N = 0.1, P = 0),
                         traits = c(R = 10, N = 4, P = 7), dt = 0.5)
  expect_equal(x$P, rep(0, 4))
  expect_equal(x$z_P, rep(7, 4))
  expect_equal(x$rate_P, rep(0, 4))
  expect_equal(x$grad_P, rep(0, 4))
  expect_true(all(x$N > 0))
  # nothing at all
  none <- food_web_dynamics(c(0, 50), init = c(R = 0, N = 0, P = 0), dt = 0.5)
  expect_equal(unlist(none[2, c("R", "N", "P")]), c(R = 0, N = 0, P = 0))
  # an abundance below `extinction` is set to zero and stays there
  dying <- food_web_dynamics(c(0, 50, 100), init = c(R = 50, N = 1e-5, P = 0),
                             traits = c(R = 10, N = 40, P = 1), dt = 0.5, extinction = 1e-4)
  expect_equal(dying$N, c(1e-5, 0, 0))
  expect_equal(dying$z_N[3], dying$z_N[2])
  expect_equal(dying$rate_N[2:3], c(0, 0))
})

test_that("halving the step changes the result little", {
  times <- c(0, 300, 1000)
  coarse <- food_web_dynamics(times, init = c(R = 50, N = 0.01, P = 0),
                              traits = c(R = 10, N = 4, P = 1), dt = 1)
  fine <- food_web_dynamics(times, init = c(R = 50, N = 0.01, P = 0),
                            traits = c(R = 10, N = 4, P = 1), dt = 0.5)
  finer <- food_web_dynamics(times, init = c(R = 50, N = 0.01, P = 0),
                             traits = c(R = 10, N = 4, P = 1), dt = 0.25)
  cols <- c("R", "N", "z_R", "z_N")
  expect_equal(fine[cols], coarse[cols], tolerance = 1e-4)
  # fourth order: the error falls by about 16 when the step halves
  err1 <- abs(coarse$z_N[3] - finer$z_N[3])
  err2 <- abs(fine$z_N[3] - finer$z_N[3])
  expect_lt(err2, err1)
  # the output times need not be multiples of the step
  odd <- food_web_dynamics(c(0, 0.3, 1.7), dt = 0.5)
  expect_equal(nrow(odd), 3)
})

test_that("traits do not move when V is zero", {
  x <- food_web_dynamics(seq(0, 300, by = 100), init = c(R = 50, N = 0.1, P = 0),
                         traits = c(R = 10, N = 4, P = 1), V = 0, dt = 0.5)
  expect_equal(x$z_R, rep(10, 4))
  expect_equal(x$z_N, rep(4, 4))
  expect_equal(x$rate_R, rep(0, 4))
  expect_true(any(x$grad_N != 0))  # selection acts, but there is no variance to respond
  expect_true(all(diff(x$N) != 0))  # the ecology still runs
  one <- food_web_dynamics(seq(0, 300, by = 100), init = c(R = 50, N = 0.1, P = 0),
                           traits = c(R = 10, N = 4, P = 1), V = c(R = 0, N = 0.2, P = 0),
                           dt = 0.5)
  expect_equal(one$z_R, rep(10, 4))
  expect_true(all(diff(one$z_N) > 0))
})

test_that("introductions happen at the given time with the given state", {
  x <- food_web_dynamics(seq(0, 500, by = 50), init = c(R = 50, N = 0.01, P = 0.5),
                         traits = c(R = 10, N = 4, P = 6), introductions = c(P = 300),
                         dt = 0.5)
  expect_equal(x$P[x$time < 300], rep(0, 6))
  expect_equal(x$P[x$time == 300], 0.5)
  expect_equal(x$z_P[x$time == 300], 6)
  expect_equal(x$z_P[x$time < 300], rep(6, 6))
  expect_true(x$rate_P[x$time == 300] > 0)
  expect_equal(attr(x, "introductions"), c(P = 300))
  # before the introduction the run is the one without the predator
  without <- food_web_dynamics(seq(0, 250, by = 50), init = c(R = 50, N = 0.01, P = 0),
                               traits = c(R = 10, N = 4, P = 6), dt = 0.5)
  expect_equal(x[1:6, c("R", "N", "z_R", "z_N")], without[, c("R", "N", "z_R", "z_N")])
  # an introduction time that is not an output time still takes effect
  y <- food_web_dynamics(c(0, 250, 350), init = c(R = 50, N = 0.01, P = 0.5),
                         traits = c(R = 10, N = 4, P = 6), introductions = c(P = 300),
                         dt = 0.5)
  expect_equal(y$P[2], 0)
  expect_gt(y$P[3], 0)
  # consumer and then predator
  both <- food_web_dynamics(c(0, 100, 200, 300), init = c(R = 50, N = 0.1, P = 0.1),
                            traits = c(R = 10, N = 4, P = 6),
                            introductions = c(N = 100, P = 200), dt = 0.5)
  expect_equal(both$N[1:2], c(0, 0.1))
  expect_equal(both$P[1:3], c(0, 0, 0.1))
  expect_equal(both$R[1:2], c(50, 50))
})

test_that("the three attack forms and omnivory run", {
  for (attack in c("gaussian", "logistic", "michaelis")) {
    x <- food_web_dynamics(c(0, 50, 100), init = c(R = 50, N = 0.1, P = 0.1),
                           traits = c(R = 10, N = 4, P = 6), attack = attack, dt = 0.5)
    expect_true(all(is.finite(unlist(x))))
    expect_true(all(x$a_RN > 0))
    expect_true(all(x$m_NP > 0))
    expect_equal(x$v_RP, rep(0, 3))
  }
  omni <- food_web_dynamics(c(0, 50, 100), init = c(R = 50, N = 0.1, P = 0.1),
                            traits = c(R = 10, N = 4, P = 6), omnivory = TRUE, dt = 0.5)
  expect_true(all(omni$v_RP > 0))
  expect_true(attr(omni, "omnivory"))
  # the attack coefficients follow the book's forms at the means
  p <- food_web_params()
  g <- food_web_dynamics(c(0, 1), traits = c(R = 10, N = 4, P = 1), dt = 0.5)
  expect_equal(g$a_RN[1], p$a0 * exp(-((4 - 10) / p$beta)^2))
  l <- food_web_dynamics(c(0, 1), traits = c(R = 10, N = 4, P = 1), attack = "logistic",
                         dt = 0.5)
  expect_equal(l$a_RN[1], p$a0 / (1 + exp(-p$alpha * (4 - 10))))
  m <- food_web_dynamics(c(0, 1), traits = c(R = 10, N = 4, P = 1), attack = "michaelis",
                         dt = 0.5)
  expect_equal(m$a_RN[1], p$a0 * 10 / (p$eps_R + 10) * 4 / (p$eps_N + 4))
})

test_that("selection gradients match the derivatives of equation 3.11", {
  # consumer-resource system with Gaussian attack (equation 3.17): compare the numerical
  # gradients with the closed forms of the book
  p <- food_web_params()
  x <- food_web_dynamics(c(0, 1), init = c(R = 40, N = 3, P = 0), traits = c(R = 12, N = 7, P = 1),
                         dt = 0.5)
  R <- 40; N <- 3; zR <- 12; zN <- 7
  delta <- zN - zR
  a <- p$a0 * exp(-(delta / p$beta)^2)
  grad_N <- -2 * p$b * p$a0 * delta * exp(-(delta / p$beta)^2) * R /
    (p$beta^2 * (1 + a * p$h * R)^2) - 2 * p$f0 * p$theta * (zN - p$zf_N)
  grad_R <- -2 * p$c0 * p$gamma * (zR - p$zc_R) -
    2 * p$a0 * delta * exp(-(delta / p$beta)^2) * N / (p$beta^2 * (1 + a * p$h * R))
  expect_equal(x$grad_N[1], grad_N, tolerance = 1e-6)
  expect_equal(x$grad_R[1], grad_R, tolerance = 1e-6)
})

test_that("input errors are reported", {
  expect_error(food_web_dynamics(5), "increasing")
  expect_error(food_web_dynamics(c(0, 10, 5)), "increasing")
  expect_error(food_web_dynamics(c(0, 10), init = c(R = 1, N = 1)), "three finite numbers")
  expect_error(food_web_dynamics(c(0, 10), init = c(R = 1, N = 1, X = 1)), "named")
  expect_error(food_web_dynamics(c(0, 10), init = c(R = -1, N = 1, P = 1)), "zero or more")
  expect_error(food_web_dynamics(c(0, 10), V = -1), "zero or more")
  expect_error(food_web_dynamics(c(0, 10), V = c(1, 1)), "three finite numbers")
  expect_error(food_web_dynamics(c(0, 10), params = list(c0 = 1)), "food_web_params")
  expect_error(food_web_dynamics(c(0, 10), omnivory = NA), "TRUE or FALSE")
  expect_error(food_web_dynamics(c(0, 10), dt = 0), "positive")
  expect_error(food_web_dynamics(c(0, 10), extinction = -1), "zero or more")
  expect_error(food_web_dynamics(c(0, 10), introductions = 5), "named")
  expect_error(food_web_dynamics(c(0, 10), introductions = c(Q = 5)), "named")
  expect_error(food_web_dynamics(c(0, 10), introductions = c(P = 5, P = 6)), "each once")
  expect_error(food_web_dynamics(c(0, 10), init = c(R = 50, N = 0.1, P = 1),
                                 introductions = c(P = 0)), "after the first time")
  expect_error(food_web_dynamics(c(0, 10), init = c(R = 50, N = 0.1, P = 1),
                                 introductions = c(P = 20)), "within `times`")
  expect_error(food_web_dynamics(c(0, 10), init = c(R = 50, N = 0.1, P = 0),
                                 introductions = c(P = 5)), "positive abundance")
  expect_error(food_web_dynamics(c(0, 10), traits = c(R = -1, N = 4, P = 1),
                                 attack = "michaelis"), "michaelis")
  expect_error(food_web_dynamics(c(0, 400), init = c(R = 50, N = 0.01, P = 0),
                                 traits = c(R = 12, N = 6, P = 2),
                                 params = food_web_params(c0 = 10, d = 0.2, a0 = 1),
                                 dt = 3), "diverged")
})

test_that("subsets of the output print as data frames", {
  x <- food_web_dynamics(seq(0, 20, 10))
  expect_output(print(x), "<food_web>")
  expect_output(print(utils::tail(x[, c("time", "R")], 1)), "time")
})
