# Animation of the environment of Thompson et al. (2020) for the validation article: the field
# of env_space_time() on a 30 x 30 grid of coordinates 1 to 100 (their scales of 50 and 500
# steps), next to env_over_time() with the same temporal autocorrelation per step but each
# site on its own. Both are converted to ranks over all sites and steps, so that they share
# the colour scale. Run from the package root; needs ggplot2, gganimate and gifski, and writes
# vignettes/articles/env_space_time.gif.
library(comdynamite)
library(ggplot2)
library(gganimate)

set.seed(2020)
side <- 30
xy <- expand.grid(x = seq(1, 100, length.out = side), y = seq(1, 100, length.out = side))
n_steps <- 2200
burn_in <- 200
frames <- seq(burn_in, n_steps, by = 25)

space_time <- env_space_time(xy, n_steps = n_steps, burn_in = burn_in)
# the same autocorrelation per step (temporal scale 500), without spatial structure
independent <- env_over_time(rep(0, side^2), n_steps = n_steps - burn_in, rho = exp(-1 / 500),
                             sd = 1)
independent <- cbind(matrix(independent[, 1], side^2, burn_in), independent)
independent[] <- stats::ecdf(independent)(independent)

panel <- function(env, label) {
  do.call(rbind, lapply(frames, function(t) {
    data.frame(xy, step = t, panel = label, value = env[, t])
  }))
}
data <- rbind(panel(space_time, "env_space_time()"), panel(independent, "env_over_time()"))
data$panel <- factor(data$panel, c("env_space_time()", "env_over_time()"))

p <- ggplot(data, aes(x, y, fill = value)) +
  geom_raster() +
  facet_wrap(~panel, nrow = 1) +
  coord_equal() +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.5, 1), name = "environment") +
  guides(fill = guide_colourbar(barwidth = 10)) +
  labs(title = "Step {current_frame}") +
  theme_void(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA)) +
  transition_manual(step)
animate(p, nframes = length(frames), fps = 10, width = 520, height = 320, res = 96,
        renderer = gifski_renderer("vignettes/articles/env_space_time.gif", loop = TRUE))
