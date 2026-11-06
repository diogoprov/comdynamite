# Animations of metacommunity dynamics for the article "Watching a metacommunity
# assemble". Run from the package root; needs ggplot2, gganimate and gifski, and
# writes two GIFs to vignettes/articles/.
library(comdynamite)
library(ggplot2)
library(gganimate)

set.seed(42)
side <- 20
rates <- c("0.001" = 0.001, "0.05" = 0.05, "0.5" = 0.5)
labels <- paste("dispersal", names(rates))

tree <- ape::rcoal(30)
land <- landscape_grid(side, env_bandwidth = 3, disp_scale = 1.5)
traits <- species_traits(tree, lambda = 1, sigma_mean = 8, sigma_sd = 2)
alpha <- competition_matrix(tree, "phylo", inter_max = 0.5)
xy <- expand.grid(x = seq_len(side), y = seq_len(side))

# one map per recorded step and dispersal rate
run <- function(rate, label, frames, env_time = NULL) {
  sim <- simulate_metacommunity(land, traits, alpha, dispersal = rate, n_steps = max(frames),
                                env_time = env_time, record_steps = frames)
  do.call(rbind, lapply(seq_along(frames), function(i) {
    n <- sim$history[, , i]
    total <- rowSums(n)
    data.frame(xy, step = frames[i], panel = paste("dispersal", label),
               richness = rowSums(n > 0),
               value = ifelse(total > 0, drop(n %*% traits$optimum) / total, NA))
  }))
}

map_theme <- theme_void(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA))

animate_map <- function(data, variable, label, palette, n_panels, file) {
  p <- ggplot(data, aes(x, y, fill = .data[[variable]])) +
    geom_raster() +
    facet_wrap(~panel, nrow = 1) +
    coord_equal() +
    palette +
    labs(title = "Step {current_frame}", fill = label) +
    map_theme +
    transition_manual(step)
  animate(p, nframes = length(unique(data$step)), fps = 10, width = 240 * n_panels, height = 320,
          res = 96, renderer = gifski_renderer(file, loop = TRUE))
}

dir.create("vignettes/articles", showWarnings = FALSE)

# Assembly in a constant environment: local richness
frames <- c(1:10 * 2, seq(25, 400, 5))
assembly <- do.call(rbind, Map(run, rates, names(rates), MoreArgs = list(frames = frames)))
assembly$panel <- factor(assembly$panel, labels)
animate_map(assembly, "richness", "Species per site", scale_fill_viridis_c(limits = c(0, NA)),
            3, "vignettes/articles/dynamics_richness.gif")

# A moving environment: from step 200 the map shifts one cell every 5 steps,
# returning to its start every 100 steps
n_steps <- 400
env_map <- matrix(land$env, side)
shift <- function(k) c(env_map[, (seq_len(side) + k - 1) %% side + 1])
env_time <- sapply(seq_len(n_steps), function(t) shift(max(0, t - 200) %/% 5))
frames <- seq(180, n_steps, 2)
tracking <- do.call(rbind, Map(run, rates, names(rates),
                               MoreArgs = list(frames = frames, env_time = env_time)))
environment <- do.call(rbind, lapply(frames, function(s) {
  data.frame(xy, step = s, panel = "environment", richness = NA, value = env_time[, s])
}))
tracking <- rbind(environment, tracking)
tracking$panel <- factor(tracking$panel, c("environment", labels))
animate_map(tracking, "value", "Environment or mean optimum",
            scale_fill_viridis_c(option = "cividis", na.value = "grey90"),
            4, "vignettes/articles/dynamics_tracking.gif")

# Share of sites still occupied at the last step, and how closely the mean
# optimum follows the environment there
last <- tracking[tracking$step == n_steps & tracking$panel != "environment", ]
sapply(split(last, droplevels(last$panel)), function(d) {
  occupied <- !is.na(d$value)
  c(occupied = mean(occupied), match = cor(d$value[occupied], env_time[occupied, n_steps]))
}) |> round(2) |> print()
