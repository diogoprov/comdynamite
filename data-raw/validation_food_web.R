# food_web_dynamics() against chapter 3 of McPeek (2017), Evolutionary Community Ecology.
# Run from the root of the package; takes a few minutes.
# 1. Figure 3.2: a consumer invades a resource at its demographic and evolutionary
#    equilibrium (R* = c0 / d = 50 at z_R = 10). The book shows the consumer peaking near 10
#    at about iteration 200 while the resource falls to about 34 and recovers, then both
#    traits rising at a decelerating rate to about 17.5 (resource) and 15.5 (consumer) by
#    iteration 3000, with the attack coefficient near 0.09, c(z_R) near 0.9 and f(z_N) near
#    0.23 at the end. One iteration of the book is one unit of time here.
# 2. Figure 3.3: at iteration 3000 the consumer sits at a maximum of its fitness surface and
#    the resource at a minimum, both with log mean fitness zero.
# 3. Figure 3.11 at c0 = 10, one column per attack form: the resource evolves a trait closer
#    to its intrinsic birth optimum (zc_R) when the predator is present. Read from the figure,
#    with and without the predator: bidirectional-dependent (panel I) about 12.4 against 15.6
#    (zc_R = 12, delta = 0.02); unidirectional-dependent (panel E) about 12.9 against 15.4
#    (zc_R = 12; the consumer and predator traits of this panel cannot be reproduced from the
#    caption's parameters: at the plotted traits the consumer cannot persist for any
#    abundance of the resource up to its carrying capacity); unidirectional-independent
#    (panel A) about 11.5 against 10.1 (zc_R = 12: the caption says 20, but the resource alone
#    sits at 12 at low c0 in every panel). Initial abundance of the invaders 0.01, as in the
#    caption of figure 2 of McPeek (2017, The American Naturalist 189, E91-E117,
#    doi:10.1086/691101); traits 6 (N) and 10 (P) are our choice, as the predator does not
#    invade from a trait far from the consumer's. In panel E the predator trait is still
#    moving at iteration 18000 (the gradient of its death rate is small), so its final value
#    depends on its initial trait. The chapter's MATLAB code (enallagma.com) was offline and
#    not archived when checked in November 2026, so these comparisons are against the book.
#    The resource-consumer code of the supplement of McPeek (2017, The American Naturalist)
#    integrates with ode45 in continuous time (one iteration = one unit of time); run with the
#    parameters of figure 3.2 (outside this package, as it has no declared licence), it gives
#    the consumer peak at t = 200 and trajectories within 1e-4 of food_web_dynamics().
suppressMessages(devtools::load_all(".", quiet = TRUE))

# 1. Figure 3.2
fig32 <- food_web_dynamics(0:3000, init = c(R = 50, N = 0.01, P = 0),
                           traits = c(R = 10, N = 4, P = 1), dt = 0.1)
p <- attr(fig32, "params")
fig32$c_R <- p$c0 * (1 - p$gamma * (fig32$z_R - p$zc_R)^2)
fig32$f_N <- p$f0 * (1 + p$theta * (fig32$z_N - p$zf_N)^2)
peak <- fig32[which.max(fig32$N), ]
rate_peak <- sapply(c("rate_R", "rate_N"), function(v) fig32$time[which.max(fig32[[v]])])
summary_32 <- list(
  consumer_peak = c(time = peak$time, N = peak$N, R = peak$R),
  resource_dip = c(time = fig32$time[which.min(fig32$R[fig32$time <= 500])],
                   R = min(fig32$R[fig32$time <= 500])),
  final = unlist(fig32[nrow(fig32), c("R", "N", "z_R", "z_N", "a_RN", "c_R", "f_N")]),
  time_of_max_rate = rate_peak,
  rate_at_max_over_rate_at_3000 = c(
    R = max(fig32$rate_R) / fig32$rate_R[nrow(fig32)],
    N = max(fig32$rate_N) / fig32$rate_N[nrow(fig32)]))
print(summary_32)

# 2. Figure 3.3: curvature of the fitness surfaces at the mean, iteration 3000
state <- unlist(fig32[nrow(fig32), c("R", "N", "P", "z_R", "z_N", "z_P")])
pp <- unclass(p)
pp$v0 <- 0
sh <- comdynamite:::fw_shapes(pp, "gaussian")
lnw <- function(level, z) comdynamite:::fw_fitness(z, level, state, pp, "gaussian", sh)
curvature <- function(level, z, h = 0.05) (lnw(level, z + h) - 2 * lnw(level, z) +
                                             lnw(level, z - h)) / h^2
summary_33 <- c(lnW_R = lnw("R", state[["z_R"]]), lnW_N = lnw("N", state[["z_N"]]),
                curvature_R = curvature("R", state[["z_R"]]),
                curvature_N = curvature("N", state[["z_N"]]))
print(summary_33)

# 3. Figure 3.11 at c0 = 10
panels <- list(
  bidirectional = list(attack = "gaussian", pars = list(
    c0 = 10, d = 0.2, a0 = 1, b = 0.1, h = 0.1, m0 = 0.05, n = 0.1, l = 0.1, f0 = 0.2, g = 0,
    x0 = 0.01, y = 0, beta = 5, phi = 5, gamma = 0.01, theta = 0.01, delta = 0.02, zc_R = 12,
    zf_N = 1, zx_P = 1)),
  unidirectional_dependent = list(attack = "logistic", pars = list(
    c0 = 10, d = 0.2, a0 = 0.5, b = 0.1, h = 0.2, m0 = 0.01, n = 0.1, l = 0.2, f0 = 0.2, g = 0,
    x0 = 0.01, y = 0, alpha = 0.1, rho = 0.1, gamma = 0.01, theta = 0.01, delta = 0.01,
    zc_R = 12, zf_N = 1, zx_P = 1)),
  unidirectional_independent = list(attack = "michaelis", pars = list(
    c0 = 10, d = 0.2, a0 = 0.5, b = 0.1, h = 0.3, m0 = 0.05, n = 0.1, l = 0.3, f0 = 0.05, g = 0,
    x0 = 0.01, y = 0, eps_R = 20, eps_N = 20, eta_N = 20, eta_P = 20, gamma = 0.01,
    theta = 0.01, delta = 0.01, zc_R = 12, zf_N = 1, zx_P = 1)))
t_pred <- 6000
t_end <- 18000
fig311 <- lapply(panels, function(panel) {
  pars <- do.call(food_web_params, panel$pars)
  init <- c(R = panel$pars$c0 / panel$pars$d, N = 0.01, P = 0.01)
  traits <- c(R = panel$pars$zc_R, N = 6, P = 10)
  times <- c(0, t_pred, t_end)
  without <- food_web_dynamics(times, init = replace(init, "P", 0), traits = traits,
                               params = pars, attack = panel$attack, dt = 0.1)
  with <- food_web_dynamics(times, init = init, traits = traits, params = pars,
                            attack = panel$attack, dt = 0.1, introductions = c(P = t_pred))
  keep <- c("R", "N", "P", "z_R", "z_N", "z_P", "rate_R", "rate_N", "rate_P")
  out <- rbind(without_predator = without[2, keep], with_predator = with[3, keep])
  out$zc_R <- panel$pars$zc_R
  out
})
print(lapply(fig311, round, 3))

saveRDS(list(fig32 = fig32[fig32$time %% 10 == 0, ], summary_32 = summary_32,
             summary_33 = summary_33, fig311 = fig311),
        "data-raw/validation_food_web.rds")
