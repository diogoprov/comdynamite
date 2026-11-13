# Animations and figures for the article "Quantitative genetics in the metacommunity". Run
# from the package root; needs ggplot2, gganimate and gifski, and writes GIFs and PNGs to
# vignettes/articles/ and the numbers quoted in the text to data-raw/quantgen_numbers.rds.
library(comdynamite)
library(ggplot2)
library(gganimate)

out <- "vignettes/articles"
dir.create(out, showWarnings = FALSE)
numbers <- list()

map_theme <- theme_void(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA),
        strip.text = element_text(margin = margin(4, 0, 4, 0)))
plot_theme <- theme_bw(base_size = 13) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        legend.position = "bottom")

animate_plot <- function(p, n_frames, width, height, file) {
  animate(p, nframes = n_frames, fps = 10, width = width, height = height, res = 96,
          renderer = gifski_renderer(file.path(out, file), loop = TRUE))
}

save_plot <- function(p, file, width = 8, height = 4) {
  ggsave(file.path(out, file), p, width = width, height = height, dpi = 110, bg = "white")
}

# one species on an environmental gradient: maps of the load of each population
load_maps <- function(sim, frames, side, panel) {
  m <- maladaptation(sim, step = frames)
  xy <- expand.grid(x = seq_len(side), y = seq_len(side))
  grid <- merge(data.frame(site = rownames(sim$abundance), xy),
                data.frame(step = frames), by = NULL)
  grid <- merge(grid, m[, c("step", "site", "load")], all.x = TRUE)
  grid$panel <- panel
  grid
}

# 1. Local adaptation and migration load along a gradient ----------------------------------

set.seed(1)
side <- 20
land <- landscape_grid(side, env_bandwidth = 4, disp_scale = 1.5)
one <- data.frame(optimum = mean(land$env), sigma = 12, r_max = 5, row.names = "sp1")
a1 <- matrix(0.02, dimnames = list("sp1", "sp1"))
rates <- c(0.001, 0.05, 0.3)
frames <- c(1:10, seq(12, 60, 2), seq(65, 300, 5))
gradient <- lapply(rates, function(d) {
  simulate_metacommunity(land, one, a1, dispersal = d, n_steps = 300, record_steps = frames,
                         evolution = evolution_control(4, env_var = 2))
})
labels <- paste("dispersal", rates)
maps <- do.call(rbind, Map(load_maps, gradient, MoreArgs = list(frames = frames, side = side),
                           panel = labels))
env_panel <- merge(data.frame(site = rownames(gradient[[1]]$abundance),
                              expand.grid(x = seq_len(side), y = seq_len(side)),
                              load = land$env / 100, panel = "environment"),
                   data.frame(step = frames), by = NULL)
maps <- rbind(env_panel[, names(maps)], maps)
maps$panel <- factor(maps$panel, c("environment", labels))
p <- ggplot(maps, aes(x, y, fill = load)) +
  geom_raster() +
  facet_wrap(~panel, nrow = 1) +
  coord_equal() +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "grey85") +
  labs(title = "Step {current_frame}", fill = "Lag load (environment / 100 in the first panel)") +
  map_theme +
  transition_manual(step)
animate_plot(p, length(frames), 960, 330, "quantgen_gradient.gif")

final <- do.call(rbind, Map(function(sim, label) {
  m <- maladaptation(sim)
  m$panel <- label
  m
}, gradient, labels))
final$panel <- factor(final$panel, labels)
final$side <- ifelse(final$env < mean(land$env), "below", "above")
p <- ggplot(final, aes(env, load, colour = panel)) +
  geom_point(alpha = 0.5, size = 1.1) +
  geom_vline(xintercept = mean(land$env), linetype = 2, colour = "grey50") +
  geom_vline(data = do.call(rbind, Map(function(sim, label) {
    m <- maladaptation(sim)
    data.frame(panel = factor(label, labels),
               phenotype = stats::weighted.mean(m$phenotype, m$abundance))
  }, gradient, labels)), aes(xintercept = phenotype, colour = panel), linewidth = 0.8) +
  facet_wrap(~panel, nrow = 1) +
  scale_colour_viridis_d(end = 0.85, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Environment of the site", y = "Lag load") +
  plot_theme
save_plot(p, "quantgen_gradient.png", width = 10, height = 3.6)
no_evolution <- lapply(rates[2:3], function(d) {
  simulate_metacommunity(land, one, a1, dispersal = d, n_steps = 300,
                         evolution = evolution_control(0, env_var = 2))
})
fixed <- lapply(no_evolution, maladaptation)
numbers$gradient <- data.frame(
  scenario = c(labels, paste(labels[2:3], "G = 0")),
  mean_load = c(tapply(final$load, final$panel, mean), vapply(fixed, function(m) mean(m$load), 1)),
  load_below = c(tapply(final$load[final$side == "below"], final$panel[final$side == "below"],
                        mean),
                 vapply(fixed, function(m) mean(m$load[m$env < mean(land$env)]), 1)),
  load_above = c(tapply(final$load[final$side == "above"], final$panel[final$side == "above"],
                        mean),
                 vapply(fixed, function(m) mean(m$load[m$env >= mean(land$env)]), 1)),
  mean_phenotype = c(vapply(gradient, function(sim) {
    m <- maladaptation(sim)
    stats::weighted.mean(m$phenotype, m$abundance)
  }, 1), NA, NA),
  occupied = c(tapply(final$site, final$panel, length),
               vapply(fixed, nrow, 1L)))
numbers$regional_mean <- mean(land$env)
numbers$share_below <- mean(land$env < mean(land$env))

# 2. Evolutionary rescue under a warming environment -----------------------------------------

n_steps <- 500
warming <- 0.15
env_time <- sapply(seq_len(n_steps), function(t) land$env + warming * max(0, t - 100))
gs <- c(0, 1, 4, 16)
frames <- seq(90, n_steps, 5)
rescue <- lapply(gs, function(g) {
  set.seed(1)
  simulate_metacommunity(land, one, a1, dispersal = 0.02, n_steps = n_steps,
                         env_time = env_time, record_steps = frames,
                         evolution = evolution_control(g, env_var = 2))
})
labels <- paste("G =", gs)
maps <- do.call(rbind, Map(load_maps, rescue, MoreArgs = list(frames = frames, side = side),
                           panel = labels))
maps$panel <- factor(maps$panel, labels)
p <- ggplot(maps, aes(x, y, fill = load)) +
  geom_raster() +
  facet_wrap(~panel, nrow = 1) +
  coord_equal() +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "grey85") +
  labs(title = "Step {current_frame} (warming from step 100)", fill = "Lag load") +
  map_theme +
  transition_manual(step)
animate_plot(p, length(frames), 960, 330, "quantgen_rescue.gif")

series <- do.call(rbind, Map(function(sim, label) {
  m <- maladaptation(sim, step = frames)
  data.frame(step = frames, panel = label,
             occupied = as.vector(table(factor(m$step, frames))) / side^2,
             load = as.vector(tapply(m$load, factor(m$step, frames), mean)))
}, rescue, labels))
series$panel <- factor(series$panel, labels)
long <- rbind(data.frame(series[, c("step", "panel")], what = "Share of sites occupied",
                         value = series$occupied),
              data.frame(series[, c("step", "panel")], what = "Mean lag load",
                         value = series$load))
p <- ggplot(long, aes(step, value, colour = panel)) +
  geom_vline(xintercept = 100, linetype = 2, colour = "grey50") +
  geom_line(linewidth = 0.9) +
  facet_wrap(~what, scales = "free_y") +
  scale_colour_viridis_d(end = 0.85) +
  labs(x = "Step", y = NULL, colour = NULL) +
  plot_theme
save_plot(p, "quantgen_rescue.png")
numbers$rescue <- series[series$step %in% c(90, n_steps), ]
numbers$rescue_minimum <- do.call(rbind, lapply(split(series, series$panel), function(d) {
  d[which.min(d$occupied), ]
}))

# 3. Genetic constraints: the orientation of G -----------------------------------------------

set.seed(1)
side3 <- 12
land3 <- landscape_grid(side3, env_bandwidth = 3, disp_scale = 1)
size <- (land3$env - 50) / 4
shape <- 0.8 * size + stats::rnorm(length(size), 0, 3)
env3 <- cbind(size = size, shape = shape)
two <- list(size = data.frame(optimum = 0, sigma = 8, r_max = 4, row.names = "sp1"),
            shape = data.frame(optimum = 0, sigma = 8, r_max = 1, row.names = "sp1"))
correlations <- c(0.9, 0, -0.9)
labels <- c("aligned (r = 0.9)", "uncorrelated", "opposed (r = -0.9)")
frames <- c(1:10, seq(12, 40, 2), seq(45, 400, 5))
constrained <- lapply(correlations, function(r) {
  m <- genetic_matrices("sp1", P = 4 * matrix(c(1, r, r, 1), 2), heritability = 0.5)
  set.seed(2)
  simulate_metacommunity(land3, two, a1, env = env3, dispersal = 0.001, n_steps = 400,
                         record_steps = frames, seed_steps = numeric(0),
                         initial = matrix(50, side3^2, 1, dimnames = list(NULL, "sp1")),
                         evolution = evolution_control(m$genetic_var, env_var = m$env_var,
                                                       ne_ratio = Inf, founder = FALSE))
})
means <- do.call(rbind, Map(function(sim, label) {
  h <- sim$evolution$history
  do.call(rbind, lapply(seq_along(frames), function(i) {
    data.frame(step = frames[i], panel = label, env_size = env3[, "size"],
               size = h[, 1, "size", i], shape = h[, 1, "shape", i])
  }))
}, constrained, labels))
means$panel <- factor(means$panel, labels)
targets <- merge(data.frame(env_size = env3[, "size"], size = env3[, "size"],
                            shape = env3[, "shape"]),
                 data.frame(panel = factor(labels, labels)), by = NULL)
ellipse <- do.call(rbind, Map(function(sim, label) {
  g <- sim$evolution$genetic_var$sp1
  e <- eigen(g, symmetric = TRUE)
  t <- seq(0, 2 * pi, length.out = 100)
  xy <- cbind(cos(t), sin(t)) %*% diag(2 * sqrt(e$values)) %*% t(e$vectors)
  data.frame(panel = label, size = xy[, 1], shape = xy[, 2])
}, constrained, labels))
ellipse$panel <- factor(ellipse$panel, labels)
p <- ggplot(means, aes(size, shape)) +
  geom_point(data = targets, aes(colour = env_size), shape = 4, size = 1.6, alpha = 0.5) +
  geom_path(data = ellipse, colour = "grey30", linewidth = 0.8) +
  geom_point(aes(colour = env_size), size = 1.6) +
  facet_wrap(~panel, nrow = 1) +
  coord_equal() +
  scale_colour_viridis_c(option = "cividis") +
  labs(title = "Step {current_frame}", x = "Trait 1 (size)", y = "Trait 2 (shape)",
       colour = "Environment (size axis)") +
  plot_theme +
  transition_manual(step)
animate_plot(p, length(frames), 960, 420, "quantgen_constraints.gif")

loads <- do.call(rbind, Map(function(sim, label) {
  m <- maladaptation(sim, step = frames)
  data.frame(step = frames, panel = label,
             load = as.vector(tapply(m$load, factor(m$step, frames), mean)))
}, constrained, labels))
loads$panel <- factor(loads$panel, labels)
p <- ggplot(loads, aes(step, load, colour = panel)) +
  geom_line(linewidth = 0.9) +
  scale_y_log10() +
  scale_colour_viridis_d(end = 0.85) +
  labs(x = "Step", y = "Mean lag load (log scale)", colour = "G relative to the environment") +
  plot_theme
save_plot(p, "quantgen_constraints.png")
numbers$constraints <- loads[loads$step %in% c(10, 50, 100, 400), ]

# 4. Evolving priority effects and regional neutrality ---------------------------------------

set.seed(1)
side4 <- 15
land4 <- landscape_grid(side4, env_bandwidth = 3, disp_scale = 1)
sp <- paste0("sp", 1:5)
five <- data.frame(optimum = seq(15, 85, length.out = 5), sigma = 25, r_max = 5,
                   row.names = sp)
a5 <- matrix(0.02, 5, 5, dimnames = list(sp, sp))
frames <- seq(50, 5000, 50)
gs <- c(0, 16)
monopoly <- lapply(gs, function(g) {
  set.seed(1)
  simulate_metacommunity(land4, five, a5, dispersal = 0.002, extirpation = 0.02,
                         n_steps = 5000, record_steps = frames,
                         evolution = evolution_control(g, env_var = 2))
})
labels <- paste("G =", gs)
xy4 <- expand.grid(x = seq_len(side4), y = seq_len(side4))
dominant <- do.call(rbind, Map(function(sim, label) {
  do.call(rbind, lapply(seq_along(frames), function(i) {
    n <- sim$history[, , i]
    top <- ifelse(rowSums(n) > 0, sp[max.col(n, ties.method = "first")], NA)
    data.frame(xy4, step = frames[i], panel = label, species = top)
  }))
}, monopoly, labels))
# the environment panel shows, in each site, the species whose ancestral optimum is closest
nearest <- sp[apply(abs(outer(land4$env, five$optimum, "-")), 1, which.min)]
envmap <- merge(data.frame(xy4, species = nearest, panel = "best ancestral match"),
                data.frame(step = frames), by = NULL)
dominant <- rbind(envmap[, names(dominant)], dominant)
dominant$panel <- factor(dominant$panel, c("best ancestral match", labels))
colours <- stats::setNames(viridisLite::viridis(5, option = "turbo", begin = 0.1, end = 0.9),
                           sp)
p <- ggplot(dominant, aes(x, y, fill = species)) +
  geom_raster() +
  scale_fill_manual(values = colours, na.value = "grey85") +
  facet_wrap(~panel, nrow = 1) +
  coord_equal() +
  labs(title = "Step {current_frame}", fill = "Dominant species (ordered by ancestral optimum)") +
  map_theme +
  transition_manual(step)
animate_plot(p, length(frames), 900, 360, "quantgen_monopolization.gif")

trends <- do.call(rbind, Map(function(sim, label) {
  do.call(rbind, lapply(seq_along(frames), function(i) {
    n <- sim$history[, , i]
    z <- sim$evolution$history[, , i]
    w <- n
    w[is.na(z)] <- 0
    z[is.na(z)] <- 0
    occ <- rowSums(n) > 0
    regional <- colSums(w * z) / colSums(w)
    data.frame(step = frames[i], panel = label,
               sorting = stats::cor(land4$env[occ],
                                    drop(n %*% five$optimum)[occ] / rowSums(n)[occ]),
               local = stats::cor(land4$env[occ], rowSums(w * z)[occ] / rowSums(w)[occ]),
               spread = stats::sd(regional, na.rm = TRUE),
               dominance = mean(apply(n[occ, ] / rowSums(n[occ, ]), 1, max)),
               species = sum(colSums(n) > 0))
  }))
}, monopoly, labels))
trends$panel <- factor(trends$panel, labels)
long <- rbind(
  data.frame(trends[, c("step", "panel")], what = "Environment vs ancestral optimum (r)",
             value = trends$sorting),
  data.frame(trends[, c("step", "panel")], what = "Environment vs evolved optimum (r)",
             value = trends$local),
  data.frame(trends[, c("step", "panel")], what = "SD of regional mean optima",
             value = trends$spread))
p <- ggplot(long, aes(step, value, colour = panel)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~what, scales = "free_y") +
  scale_colour_manual(values = c("grey40", "#B8336A")) +
  labs(x = "Step", y = NULL, colour = NULL) +
  plot_theme
save_plot(p, "quantgen_monopolization.png", width = 10, height = 3.8)
numbers$monopoly <- trends[trends$step %in% c(250, 1000, 2500, 5000), ]
late <- trends[trends$step > 4000, ]
numbers$monopoly_late <- aggregate(cbind(sorting, local, spread, dominance) ~ panel, late, mean)
numbers$best_match <- vapply(monopoly, function(sim) {
  n <- sim$abundance
  occ <- rowSums(n) > 0
  mean(sp[max.col(n, ties.method = "first")][occ] == nearest[occ])
}, numeric(1))

# range filling: share of the sites where the ancestral optimum of each species would grow
# (growth above 1 after extirpation, ignoring competition) that the species occupies, and share
# of its occupied sites outside that potential range
filling <- do.call(rbind, Map(function(sim, label) {
  r <- growth_rates(land4$env, five)
  potential <- r * (1 - 0.02) > 1
  n <- sim$abundance
  data.frame(panel = label, species = sp,
             realised_in_potential = colSums(n > 0 & potential) / colSums(potential),
             outside_potential = colSums(n > 0 & !potential) / pmax(colSums(n > 0), 1),
             potential_sites = colSums(potential))
}, monopoly, labels))
numbers$filling <- filling

saveRDS(numbers, "data-raw/quantgen_numbers.rds")
print(numbers)
