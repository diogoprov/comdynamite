# Design figure and animation for the article "Persistence in river networks under a flood
# pulse". Run from the package root; needs ggplot2, patchwork, gganimate and gifski. Uses the
# set-up of data-raw/river_pantanal_setup.R (replicate 1, symmetric dispersal along the channel)
# and writes vignettes/articles/river_design.png and vignettes/articles/river_pulse.gif.
source("data-raw/river_pantanal_setup.R")
library(ggplot2)
library(patchwork)
library(gganimate)

s <- pantanal_setup(1)
net <- s$net
zone <- ifelse(!s$plain, "plateau", ifelse(s$channel, "main channel", "floodplain lake"))
zone <- factor(zone, c("plateau", "main channel", "floodplain lake"))
zone_cols <- c(plateau = "#B07A3B", "main channel" = "#1F5FA8", "floodplain lake" = "#5FB3C9")
outlet <- which(is.na(net$downstream))
to <- ifelse(is.na(net$downstream), outlet, net$downstream)
seg <- data.frame(x = net$x, y = net$y, xend = net$x[to], yend = net$y[to])[-outlet, ]
hull <- net[s$plain, c("x", "y")]
hull <- hull[grDevices::chull(hull), ]
base_theme <- theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(face = "bold"))

# (a) the network and its zones
p_net <- ggplot() +
  geom_polygon(data = hull, aes(x, y), fill = "#DCEFF5", colour = NA) +
  geom_segment(data = seg, aes(x, y, xend = xend, yend = yend), colour = "grey55") +
  geom_point(data = cbind(net, zone = zone), aes(x, y, colour = zone), size = 3) +
  annotate("text", x = net$x[outlet], y = net$y[outlet] - 0.6,
           label = "outlet", size = 3.5) +
  scale_colour_manual(values = zone_cols, name = NULL) +
  coord_equal() +
  labs(title = "a  Network", subtitle = "shaded: area linked overland at high water") +
  base_theme +
  theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom")

# (b) the hydrological year
months <- month.abb
cal <- rbind(
  data.frame(process = "plateau: spates (washout 0.05 per site and month)", month = rain),
  data.frame(process = "migrants: piracema (emigration 0.05 to 0.3)", month = piracema),
  data.frame(process = "floodplain: overland links; others emigrate 3x", month = high),
  data.frame(process = "floodplain lakes: drying (0.03; 0.15 in drought)", month = low))
cal_cols <- stats::setNames(c("#B07A3B", "#7B4FA0", "#1F5FA8", "#E3A33B"), unique(cal$process))
cal$process <- factor(cal$process, rev(unique(cal$process)))
p_cal <- ggplot(cal, aes(factor(month, 1:12, months), process, fill = process)) +
  geom_tile(colour = "white", height = 0.7) +
  scale_fill_manual(values = cal_cols, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  labs(title = "b  Hydrological year (monthly steps)", x = NULL, y = NULL) +
  base_theme +
  theme(panel.grid = element_blank())

# (c) the regimes after 100 years of assembly under the normal pulse
reg <- rbind(
  data.frame(regime = "normal", start = -10, end = 60, pulse = "annual pulse"),
  data.frame(regime = "drought", start = c(-10, 20, 32), end = c(20, 32, 60),
             pulse = c("annual pulse", "no pulse, harsher drying", "annual pulse")),
  data.frame(regime = "regulated", start = c(-10, 0), end = c(0, 60),
             pulse = c("annual pulse", "no pulse, no drying")))
reg$regime <- factor(reg$regime, c("regulated", "drought", "normal"))
p_reg <- ggplot(reg, aes(xmin = start, xmax = end, ymin = as.numeric(regime) - 0.3,
                         ymax = as.numeric(regime) + 0.3, fill = pulse)) +
  geom_rect() +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_y_continuous(breaks = 1:3, labels = levels(reg$regime)) +
  scale_fill_manual(values = c("annual pulse" = "#1F5FA8", "no pulse, harsher drying" = "#E3A33B",
                               "no pulse, no drying" = "grey60"), name = NULL) +
  labs(title = "c  Regimes", x = "Year (after 100 years of assembly)", y = NULL) +
  base_theme +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

design <- p_net + (p_cal / p_reg) + plot_layout(widths = c(1, 1.5))
ggsave("vignettes/articles/river_design.png", design, width = 12, height = 6.5, dpi = 110)

# animation: local richness month by month from year 17 to year 38, across the drought
frames <- burn_years * period + seq(16 * period + 1, 38 * period)
regimes <- c("normal", "drought", "regulated")
frame_data <- do.call(rbind, lapply(regimes, function(g) {
  sim <- pantanal_run(s, "symmetric", g, record_steps = frames)
  sched <- pantanal_schedule(s, g)
  do.call(rbind, lapply(seq_along(frames), function(i) {
    t <- frames[i]
    data.frame(net[, c("x", "y")], zone = zone, regime = g, frame = i,
               richness = rowSums(sim$history[, , i] > 0),
               hit = sched$disturbance[, t] > 0, flood = sched$flood_on[t],
               label = sprintf("Year %d, %s", s$year[t], month.abb[s$month[t]]))
  }))
}))
frame_data$regime <- factor(frame_data$regime, regimes,
                            c("normal pulse", "12-year drought (years 21-32)",
                              "regulated (no pulse from year 1)"))
labels <- unique(frame_data[, c("frame", "label")])$label
hull_frames <- do.call(rbind, lapply(levels(frame_data$regime), function(g) {
  on <- unique(frame_data[frame_data$regime == g, c("frame", "flood")])
  on <- on[on$flood, ]
  if (!nrow(on)) return(NULL)
  do.call(rbind, lapply(on$frame, function(f) cbind(hull, frame = f, regime = g)))
}))
hull_frames$regime <- factor(hull_frames$regime, levels(frame_data$regime))
segs <- do.call(rbind, lapply(levels(frame_data$regime), function(g) cbind(seg, regime = g)))
segs$regime <- factor(segs$regime, levels(frame_data$regime))

p_anim <- ggplot() +
  geom_polygon(data = hull_frames, aes(x, y), fill = "#CDE8F2") +
  geom_segment(data = segs, aes(x, y, xend = xend, yend = yend), colour = "grey60") +
  geom_point(data = frame_data, aes(x, y, size = richness, fill = richness), shape = 21,
             colour = "grey25") +
  geom_point(data = frame_data[frame_data$hit, ], aes(x, y), shape = 4, size = 4,
             colour = "#C0392B", stroke = 1.3) +
  facet_wrap(~regime, nrow = 1) +
  scale_size_continuous(range = c(0.8, 6), limits = c(0, 16), name = "Species per site") +
  scale_fill_viridis_c(limits = c(0, 16), name = "Species per site") +
  coord_equal() +
  labs(title = "{labels[as.integer(current_frame)]}",
       subtitle = "Blue: floodplain linked overland at high water. Red cross: site lost all populations (spate or drying).") +
  base_theme +
  theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", plot.subtitle = element_text(size = 9)) +
  guides(size = "legend", fill = "legend") +
  transition_manual(frame)
animate(p_anim, nframes = length(frames), fps = 6, width = 1050, height = 470, res = 96,
        renderer = gifski_renderer("vignettes/articles/river_pulse.gif", loop = TRUE))
