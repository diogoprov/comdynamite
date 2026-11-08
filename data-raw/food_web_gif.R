# GIF shown on the reference page of food_web_dynamics(): a consumer invades the resource and a
# predator arrives at time 6000, with the parameters of the Gaussian (bidirectional) panel of
# figure 3.11 of McPeek (2017) at c0 = 10, as in data-raw/validation_food_web.R.
# Needs gifski. Run from the root of the package.
devtools::load_all(".", quiet = TRUE)
pars <- food_web_params(c0 = 10, d = 0.2, a0 = 1, b = 0.1, h = 0.1, m0 = 0.05, n = 0.1,
                        l = 0.1, f0 = 0.2, g = 0, x0 = 0.01, y = 0, beta = 5, phi = 5,
                        gamma = 0.01, theta = 0.01, delta = 0.02, zc_R = 12, zf_N = 1,
                        zx_P = 1)
x <- food_web_dynamics(seq(0, 12000, by = 50), init = c(R = 50, N = 0.01, P = 0.01),
                       traits = c(R = 12, N = 6, P = 10), params = pars,
                       introductions = c(P = 6000), dt = 0.1)
print(tail(x, 1)[, c("R", "N", "P", "z_R", "z_N", "z_P")])
animate_food_web(x, "man/figures/food_web_dynamics.gif", frames = 80, delay = 0.08,
                 width = 600, height = 600)
png(file.path(Sys.getenv("OUT", tempdir()), "fw_gauss.png"), 600, 600)
plot(x)
dev.off()
