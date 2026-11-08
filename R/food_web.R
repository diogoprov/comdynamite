#' Parameters of the three-trophic-level eco-evolutionary model
#'
#' Collects the parameters of [food_web_dynamics()], named as in chapter 3 of McPeek (2017),
#' with the values of its figure 3.2 for the resource and the consumer and of its figure 3.11
#' (bidirectional-dependent traits) for the predator. Greek letters are spelt out and the
#' optima are prefixed with `z`: `zc_R` is the trait value at which the intrinsic birth rate
#' of the resource is highest, `zf_N` and `zx_P` the values at which the intrinsic death rates
#' of the consumer and the predator are lowest.
#'
#' The intrinsic rates are quadratic in the traits (equations 3.12, 3.13 and 3.21):
#' \eqn{c(z_R) = c_0 (1 - \gamma (z_R - zc_R)^2)},
#' \eqn{f(z_N) = f_0 (1 + \theta (z_N - zf_N)^2)} and
#' \eqn{x(z_P) = x_0 (1 + \delta (z_P - zx_P)^2)}. Each feeding link has an attack coefficient
#' with a maximum (`a0` for the consumer on the resource, `m0` for the predator on the
#' consumer, `v0` for the predator on the resource) and a parameter that sets the strength of
#' selection on it, which one depends on the `attack` form chosen in [food_web_dynamics()]:
#'
#' * `"gaussian"` (bidirectional-dependent traits, equation 3.16 and box 3.1): the attack
#'   coefficient is highest when the traits match and declines as
#'   \eqn{a_0 \exp(-((z_N - z_R) / \beta)^2)}; the widths are `beta`, `phi` and `psi`.
#' * `"logistic"` (unidirectional-dependent traits, equation 3.15): the attack coefficient
#'   rises with the difference between the traits of the consumer and its prey, as
#'   \eqn{a_0 / (1 + \exp(-\alpha (z_N - z_R)))}; the steepnesses are `alpha`, `rho` and
#'   `tau`.
#' * `"michaelis"` (unidirectional-independent traits, equation 3.14): each trait raises the
#'   attack coefficient on its own, as
#'   \eqn{a_0 \, z_R / (\epsilon_R + z_R) \cdot z_N / (\epsilon_N + z_N)}; the half-saturation
#'   constants are `eps_R` and `eps_N`, `eta_N` and `eta_P`, `kappa_R` and `kappa_P`. Traits
#'   must then be positive.
#'
#' @param c0,gamma,zc_R Maximum intrinsic birth rate of the resource, the strength of
#'   stabilizing selection on it and its optimum.
#' @param d Density-dependent death rate of the resource.
#' @param a0,b,h Maximum attack coefficient, conversion efficiency and handling time of the
#'   consumer feeding on the resource.
#' @param f0,theta,zf_N Minimum intrinsic death rate of the consumer, the strength of
#'   stabilizing selection on it and its optimum.
#' @param g Density-dependent death rate of the consumer.
#' @param m0,n,l Maximum attack coefficient, conversion efficiency and handling time of the
#'   predator feeding on the consumer.
#' @param v0,w,u The same for the predator feeding on the resource, used only with
#'   `omnivory = TRUE`.
#' @param x0,delta,zx_P Minimum intrinsic death rate of the predator, the strength of
#'   stabilizing selection on it and its optimum.
#' @param y Density-dependent death rate of the predator.
#' @param beta,phi,psi Widths of the Gaussian attack coefficients of the consumer on the
#'   resource, the predator on the consumer and the predator on the resource.
#' @param alpha,rho,tau Steepnesses of the logistic attack coefficients, in the same order.
#' @param eps_R,eps_N,eta_N,eta_P,kappa_R,kappa_P Half-saturation constants of the
#'   Michaelis-Menten attack coefficients: of the resource and consumer traits for the
#'   consumer on the resource, of the consumer and predator traits for the predator on the
#'   consumer, and of the resource and predator traits for the predator on the resource.
#' @return A named list of class `food_web_params`.
#' @references McPeek, M. A. (2017). *Evolutionary Community Ecology*. Princeton University
#'   Press. \doi{10.23943/princeton/9780691088778.001.0001}
#' @seealso [food_web_dynamics()]
#' @examples
#' p <- food_web_params(c0 = 10, d = 0.2)
#' p$c0
#' @export
food_web_params <- function(c0 = 1, gamma = 0.002, zc_R = 10, d = 0.02,
                            a0 = 0.1, b = 0.1, h = 0.1,
                            f0 = 0.1, theta = 0.01, zf_N = 4, g = 0,
                            m0 = 0.05, n = 0.1, l = 0.1,
                            v0 = 0.001, w = 0.1, u = 0.1,
                            x0 = 0.01, delta = 0.01, zx_P = 1, y = 0,
                            beta = 5, phi = 5, psi = 5,
                            alpha = 0.1, rho = 0.1, tau = 0.1,
                            eps_R = 20, eps_N = 20, eta_N = 20, eta_P = 20,
                            kappa_R = 20, kappa_P = 20) {
  p <- list(c0 = c0, gamma = gamma, zc_R = zc_R, d = d, a0 = a0, b = b, h = h,
            f0 = f0, theta = theta, zf_N = zf_N, g = g, m0 = m0, n = n, l = l,
            v0 = v0, w = w, u = u, x0 = x0, delta = delta, zx_P = zx_P, y = y,
            beta = beta, phi = phi, psi = psi, alpha = alpha, rho = rho, tau = tau,
            eps_R = eps_R, eps_N = eps_N, eta_N = eta_N, eta_P = eta_P,
            kappa_R = kappa_R, kappa_P = kappa_P)
  for (nm in names(p)) {
    if (!is.numeric(p[[nm]]) || length(p[[nm]]) != 1 || !is.finite(p[[nm]])) {
      stop("`", nm, "` must be one finite number", call. = FALSE)
    }
  }
  positive <- c("beta", "phi", "psi", "eps_R", "eps_N", "eta_N", "eta_P", "kappa_R", "kappa_P")
  for (nm in positive) {
    if (p[[nm]] <= 0) stop("`", nm, "` must be positive", call. = FALSE)
  }
  any_sign <- c("zc_R", "zf_N", "zx_P", "alpha", "rho", "tau")
  for (nm in setdiff(names(p), c(positive, any_sign))) {
    if (p[[nm]] < 0) stop("`", nm, "` must be zero or more", call. = FALSE)
  }
  structure(p, class = "food_web_params")
}

#' Coevolution of a resource, a consumer and a predator
#'
#' Simulates the joint dynamics of abundances and mean traits of a basal resource (`R`), a
#' consumer (`N`) and a predator (`P`), the three-trophic-level model of chapter 3 of McPeek
#' (2017). Each level has one quantitative trait whose mean changes at a rate proportional to
#' the selection gradient, the slope of the logarithm of mean fitness against the trait at the
#' mean (equation 3.7, after Lande 1982 and Iwasa et al. 1991 as cited there):
#' \eqn{d\bar z_i / dt = V_i \, \partial \ln W_i / \partial z_i}, where \eqn{V_i} is the
#' additive genetic variance. Abundances follow the per capita growth rates of equations 3.9
#' and 3.20: logistic growth of the resource, Holling type II feeding by the consumer on the
#' resource and by the predator on the consumer and, with `omnivory = TRUE`, on the resource
#' too, with density-dependent death at every level. The intrinsic rates are quadratic in the
#' traits and the attack coefficients take one of the three forms of [food_web_params()].
#'
#' The fitness of an individual with phenotype `z` is evaluated with every other individual at
#' the mean trait (equations 3.8 and 3.20): the satiation of a consumer depends on its own
#' trait, so the trait of an individual enters both the numerator and the denominator of the
#' functional responses of the links where it feeds, but only the numerator of the links where
#' it is eaten. Selection gradients are central differences of this fitness with respect to
#' `z` at the mean, so the three attack forms share one implementation. Abundances and means
#' are integrated together by a fourth-order Runge-Kutta method with a fixed step `dt`, which
#' must be short against the fastest rate of the system (about a tenth of `1 / c0` or less):
#' halving `dt` should leave the result unchanged, and a run whose state stops being finite
#' or whose abundances overshoot below zero stops with an error.
#'
#' A level with abundance zero is absent: its trait is held and its rate of change is zero.
#' Levels can be added during the run with `introductions`, to follow the burst of evolution
#' after an invasion and its deceleration towards equilibrium (figures 3.2 and 3.3 of the
#' book), and to compare the trait of the resource with and without the predator (figure
#' 3.11). Abundances that fall below `extinction` are set to zero.
#'
#' Two details of equation 3.20 are read as follows. Its loss term of the resource to the
#' predator carries the conversion efficiency `w`, which the matching selection gradient of
#' the resource, the consumer's loss term and equation 3.9 do not; the implementation leaves
#' `w` out of the losses of prey, as in the other equations. The attack coefficient of the
#' predator on the consumer is written there with the trait of the resource; it is taken as a
#' function of the traits of the consumer and the predator, as box 3.1 defines it.
#'
#' @param times Times at which to return the state, increasing; the run starts at the first.
#' @param init Abundances at the first time, or at introduction for the levels in
#'   `introductions`, named `R`, `N` and `P` or in that order. Zero means absent.
#' @param traits Mean traits at the first time, or at introduction, in the same order.
#' @param params Parameters from [food_web_params()].
#' @param attack Form of the attack coefficients: `"gaussian"`, `"logistic"` or
#'   `"michaelis"`; see [food_web_params()].
#' @param V Additive genetic variances of the three traits, one number or three; zero means a
#'   trait that does not evolve.
#' @param omnivory Does the predator also feed on the resource (intraguild predation, with
#'   `v0`, `w` and `u` of `params`)?
#' @param introductions Optional named numbers: the times at which levels `N` or `P` (or `R`)
#'   appear, after the first time. Until then they are absent whatever `init` says.
#' @param dt Step of the integration; the intervals between the times of `times` and of
#'   `introductions` are divided into steps no longer than this.
#' @param extinction Abundance below which a level is set to zero.
#' @return A data frame of class `food_web` with one row per time in `times` and columns
#'   `time`; the abundances `R`, `N`, `P`; the mean traits `z_R`, `z_N`, `z_P`; the selection
#'   gradients `grad_R`, `grad_N`, `grad_P` (\eqn{\partial \ln W / \partial z} at the mean, so
#'   that `V * grad` is \eqn{d\bar z / dt}); the rates of trait change `rate_R`, `rate_N`,
#'   `rate_P` (\eqn{|d\bar z / dt|}); and the realized attack coefficients `a_RN`, `m_NP` and
#'   `v_RP`. The parameters, the attack form, `V`, `omnivory` and `introductions` are kept as
#'   attributes.
#' @references McPeek, M. A. (2017). *Evolutionary Community Ecology*. Princeton University
#'   Press. \doi{10.23943/princeton/9780691088778.001.0001}
#' @seealso [food_web_params()]; [trait_rates()] for rates that vary among branches of a
#'   tree.
#' @examples
#' # a consumer invades a resource at its demographic and evolutionary equilibrium
#' # (figure 3.2 of McPeek 2017)
#' x <- food_web_dynamics(seq(0, 3000, by = 50), init = c(R = 50, N = 0.1, P = 0),
#'                        traits = c(R = 10, N = 4, P = 1), dt = 0.5)
#' x[c(1, 5, 11, 61), c("time", "R", "N", "z_R", "z_N", "rate_R", "rate_N")]
#' # the predator arrives later
#' y <- food_web_dynamics(seq(0, 4000, by = 50), init = c(R = 50, N = 0.1, P = 0.1),
#'                        traits = c(R = 10, N = 4, P = 6), introductions = c(P = 2000),
#'                        dt = 0.5)
#' tail(y, 1)[, c("R", "N", "P", "z_R", "z_N", "z_P")]
#' @export
food_web_dynamics <- function(times, init = c(R = 50, N = 0.1, P = 0),
                              traits = c(R = 10, N = 4, P = 1), params = food_web_params(),
                              attack = c("gaussian", "logistic", "michaelis"),
                              V = c(R = 0.2, N = 0.2, P = 0.2), omnivory = FALSE,
                              introductions = NULL, dt = 0.1, extinction = 1e-6) {
  attack <- match.arg(attack)
  levels <- c("R", "N", "P")
  if (!is.numeric(times) || length(times) < 2 || anyNA(times) || any(diff(times) <= 0)) {
    stop("`times` must be at least two increasing numbers", call. = FALSE)
  }
  init <- fw_three(init, "init", levels)
  traits <- fw_three(traits, "traits", levels)
  if (length(V) == 1) V <- rep(V, 3)
  V <- fw_three(V, "V", levels)
  if (any(init < 0)) stop("`init` must be zero or more", call. = FALSE)
  if (any(V < 0)) stop("`V` must be zero or more", call. = FALSE)
  if (!inherits(params, "food_web_params")) {
    stop("`params` must come from `food_web_params()`", call. = FALSE)
  }
  if (!is.logical(omnivory) || length(omnivory) != 1 || is.na(omnivory)) {
    stop("`omnivory` must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(dt) || length(dt) != 1 || !is.finite(dt) || dt <= 0) {
    stop("`dt` must be one positive number", call. = FALSE)
  }
  if (!is.numeric(extinction) || length(extinction) != 1 || !is.finite(extinction) ||
      extinction < 0) {
    stop("`extinction` must be one number, zero or more", call. = FALSE)
  }
  if (!is.null(introductions)) {
    if (!is.numeric(introductions) || is.null(names(introductions)) ||
        !all(names(introductions) %in% levels) || anyDuplicated(names(introductions)) ||
        anyNA(introductions)) {
      stop("`introductions` must be numbers named `R`, `N` or `P`, each once", call. = FALSE)
    }
    if (any(introductions <= times[1]) || any(introductions > times[length(times)])) {
      stop("`introductions` must fall after the first time and within `times`", call. = FALSE)
    }
    if (any(init[names(introductions)] <= 0)) {
      stop("levels in `introductions` need a positive abundance in `init`", call. = FALSE)
    }
  }
  if (attack == "michaelis" && any(traits < 0)) {
    stop("`traits` must be zero or more with `attack = \"michaelis\"`", call. = FALSE)
  }
  p <- unclass(params)
  if (!omnivory) p$v0 <- 0

  state <- c(init, traits)
  for (nm in names(introductions)) state[nm] <- 0
  grid <- sort(unique(c(times, introductions)))
  out <- matrix(NA_real_, length(times), 15)
  colnames(out) <- c(levels, paste0("z_", levels), paste0("grad_", levels),
                     paste0("rate_", levels), "a_RN", "m_NP", "v_RP")
  row <- 0
  for (k in seq_along(grid)) {
    t_now <- grid[k]
    if (k > 1) {
      span <- t_now - grid[k - 1]
      steps <- max(1, ceiling(span / dt - 1e-8))
      h <- span / steps
      for (s in seq_len(steps)) {
        state <- fw_rk4(state, h, p, attack, V)
        if (any(!is.finite(state)) || any(state[1:3] < -extinction)) {
          stop("the integration diverged before time ", format(t_now), "; reduce `dt`",
               call. = FALSE)
        }
        low <- state[1:3] < extinction
        state[1:3][low] <- 0
      }
    }
    arrive <- names(introductions)[introductions == t_now]
    for (nm in arrive) {
      state[nm] <- init[nm]
      state[paste0("z_", nm)] <- traits[paste0("z_", nm)]
    }
    if (t_now %in% times) {
      row <- row + 1
      out[row, ] <- fw_report(state, p, attack, V)
    }
  }
  res <- data.frame(time = times, out)
  structure(res, class = c("food_web", "data.frame"), params = params, attack = attack,
            V = V, omnivory = omnivory, introductions = introductions)
}

#' @export
print.food_web <- function(x, ...) {
  # a subset (head(), tail(), selected columns) loses the attributes: print it as a data frame
  needed <- c("time", "R", "N", "P", "z_R", "z_N", "z_P")
  if (is.null(attr(x, "attack")) || !all(needed %in% names(x))) {
    return(print(as.data.frame(unclass(x)), ...))
  }
  last <- x[nrow(x), ]
  cat("<food_web>", nrow(x), "times from", format(x$time[1]), "to", format(last$time),
      "\n")
  cat("attack:", attr(x, "attack"), if (isTRUE(attr(x, "omnivory"))) "with omnivory" else "",
      "\n")
  cat("final abundances: R", format(last$R, digits = 4), " N", format(last$N, digits = 4),
      " P", format(last$P, digits = 4), "\n")
  cat("final traits:     R", format(last$z_R, digits = 4), " N", format(last$z_N, digits = 4),
      " P", format(last$z_P, digits = 4), "\n")
  invisible(x)
}

fw_three <- function(x, name, levels) {
  if (!is.numeric(x) || length(x) != 3 || anyNA(x) || any(!is.finite(x))) {
    stop("`", name, "` must be three finite numbers, for R, N and P", call. = FALSE)
  }
  if (is.null(names(x))) {
    names(x) <- levels
  } else if (!setequal(names(x), levels)) {
    stop("`", name, "` must be named `R`, `N` and `P`", call. = FALSE)
  }
  x <- x[levels]
  if (name == "traits") names(x) <- paste0("z_", levels)
  x
}

# attack coefficient of a feeder with trait z_feed on prey with trait z_prey
fw_attack <- function(z_prey, z_feed, max, shape, attack) {
  if (max == 0) return(0)
  switch(attack,
    gaussian = max * exp(-((z_feed - z_prey) / shape[1])^2),
    logistic = max / (1 + exp(-shape[1] * (z_feed - z_prey))),
    michaelis = max * z_prey / (shape[1] + z_prey) * z_feed / (shape[2] + z_feed)
  )
}

fw_shapes <- function(p, attack) {
  switch(attack,
    gaussian = list(RN = p$beta, NP = p$phi, RP = p$psi),
    logistic = list(RN = p$alpha, NP = p$rho, RP = p$tau),
    michaelis = list(RN = c(p$eps_R, p$eps_N), NP = c(p$eta_N, p$eta_P),
                     RP = c(p$kappa_R, p$kappa_P))
  )
}

# log mean fitness of an individual of `level` with phenotype z, everyone else at the mean
# (equations 3.8 and 3.20)
fw_fitness <- function(z, level, state, p, attack, sh) {
  R <- state[["R"]]; N <- state[["N"]]; P <- state[["P"]]
  zR <- state[["z_R"]]; zN <- state[["z_N"]]; zP <- state[["z_P"]]
  a_bar <- fw_attack(zR, zN, p$a0, sh$RN, attack)
  m_bar <- fw_attack(zN, zP, p$m0, sh$NP, attack)
  v_bar <- fw_attack(zR, zP, p$v0, sh$RP, attack)
  den_P <- 1 + v_bar * p$u * R + m_bar * p$l * N
  switch(level,
    R = {
      a_z <- fw_attack(z, zN, p$a0, sh$RN, attack)
      v_z <- fw_attack(z, zP, p$v0, sh$RP, attack)
      p$c0 * (1 - p$gamma * (z - p$zc_R)^2) - p$d * R -
        a_z * N / (1 + a_bar * p$h * R) - v_z * P / den_P
    },
    N = {
      a_z <- fw_attack(zR, z, p$a0, sh$RN, attack)
      m_z <- fw_attack(z, zP, p$m0, sh$NP, attack)
      p$b * a_z * R / (1 + a_z * p$h * R) - m_z * P / den_P -
        p$f0 * (1 + p$theta * (z - p$zf_N)^2) - p$g * N
    },
    P = {
      m_z <- fw_attack(zN, z, p$m0, sh$NP, attack)
      v_z <- fw_attack(zR, z, p$v0, sh$RP, attack)
      (p$w * v_z * R + p$n * m_z * N) / (1 + v_z * p$u * R + m_z * p$l * N) -
        p$x0 * (1 + p$delta * (z - p$zx_P)^2) - p$y * P
    }
  )
}

# per capita growth rates and selection gradients at the current state
fw_gradients <- function(state, p, attack, sh) {
  levels <- c("R", "N", "P")
  growth <- numeric(3)
  grad <- numeric(3)
  for (i in 1:3) {
    if (state[i] <= 0) next
    z <- state[[3 + i]]
    growth[i] <- fw_fitness(z, levels[i], state, p, attack, sh)
    step <- 1e-4 * max(1, abs(z))
    grad[i] <- (fw_fitness(z + step, levels[i], state, p, attack, sh) -
                  fw_fitness(z - step, levels[i], state, p, attack, sh)) / (2 * step)
  }
  list(growth = growth, grad = grad)
}

fw_deriv <- function(state, p, attack, V) {
  g <- fw_gradients(state, p, attack, fw_shapes(p, attack))
  c(state[1:3] * g$growth, V * g$grad)
}

fw_rk4 <- function(state, h, p, attack, V) {
  present <- state[1:3] > 0
  k1 <- fw_deriv(state, p, attack, V)
  k2 <- fw_deriv(state + h / 2 * k1, p, attack, V)
  k3 <- fw_deriv(state + h / 2 * k2, p, attack, V)
  k4 <- fw_deriv(state + h * k3, p, attack, V)
  new <- state + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
  new[1:3][!present] <- 0
  new[4:6][!present] <- state[4:6][!present]
  new
}

fw_report <- function(state, p, attack, V) {
  sh <- fw_shapes(p, attack)
  g <- fw_gradients(state, p, attack, sh)
  c(state, g$grad, abs(V * g$grad),
    fw_attack(state[["z_R"]], state[["z_N"]], p$a0, sh$RN, attack),
    fw_attack(state[["z_N"]], state[["z_P"]], p$m0, sh$NP, attack),
    fw_attack(state[["z_R"]], state[["z_P"]], p$v0, sh$RP, attack))
}
