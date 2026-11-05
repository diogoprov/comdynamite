# Four hex sticker options for comdynamite (solid colours), drawn with ggplot2
# in the geometry of hexSticker (hexagon of radius 1 centred at (1, 1), 43.9 x
# 50.8 mm). Usage: Rscript data-raw/hex_options.R out_dir
# Sourced by data-raw/hex_sticker.R with FRAME <- FALSE, which drops the hexagon
# and the name so that hexSticker draws them.
source("data-raw/hex_helpers.R")
out <- if (FRAME) commandArgs(TRUE)[1] else NULL
set.seed(7)

# 1. Dynamite stick whose fuse is a phylogeny, with sparks at the tips
tr <- ape::rcoal(6)
ts <- tree_segments(tr, c(0.86, 1.14, 1.32, 1.66))
# put the root of the fuse over the middle of the stick
r0 <- which(ts$seg$y == min(ts$seg$y) & ts$seg$yend == min(ts$seg$y)); root_x <- (ts$seg$x[r0] + ts$seg$xend[r0]) / 2
ts$seg$x <- ts$seg$x - root_x + 1; ts$seg$xend <- ts$seg$xend - root_x + 1
ts$tips$x <- ts$tips$x - root_x + 1
spark <- do.call(rbind, lapply(seq_len(nrow(ts$tips)), function(i) {
  ang <- seq(0, pi, length.out = 4) + i / 3
  data.frame(x = ts$tips$x[i], y = ts$tips$y[i],
             xend = ts$tips$x[i] + 0.03 * cos(ang), yend = ts$tips$y[i] + 0.03 * sin(ang))
}))
g1 <- base("#1D3557", "#E63946") +
  annotate("rect", xmin = 0.86, xmax = 1.14, ymin = 0.72, ymax = 1.32, fill = "#E63946") +
  annotate("rect", xmin = 0.86, xmax = 1.14, ymin = c(0.82, 1.18), ymax = c(0.86, 1.22),
           fill = "#1D3557") +
  geom_segment(aes(x, y, xend = xend, yend = yend), ts$seg, colour = "#F1FAEE",
               linewidth = 0.45, lineend = "square") +
  geom_segment(aes(x, y, xend = xend, yend = yend), spark, colour = "#FFB703", linewidth = 0.35) +
  name("#F1FAEE", y = 0.5)
save(g1, "option1_dynamite.png")

# 2. Metacommunity on a grid: sites as patches, species as coloured dots,
#    dispersal arrows; one solid colour per species
g2_sites <- expand.grid(i = 0:3, j = 0:2)
g2_sites$x <- 0.55 + g2_sites$i * 0.3
g2_sites$y <- 0.9 + g2_sites$j * 0.27
cols2 <- c("#F4A261", "#2A9D8F", "#E9C46A", "#E76F51")
dots <- do.call(rbind, lapply(seq_len(nrow(g2_sites)), function(s) {
  k <- sample(1:4, sample(1:3, 1))
  data.frame(x = g2_sites$x[s] + c(-0.05, 0.05, 0)[seq_along(k)],
             y = g2_sites$y[s] + c(0.03, 0.03, -0.05)[seq_along(k)], sp = factor(k, 1:4))
}))
g2 <- base("#264653", "#E9C46A") +
  geom_tile(aes(x, y), g2_sites, width = 0.24, height = 0.22, fill = "#2F5D6E") +
  geom_point(aes(x, y, colour = sp), dots, size = 1.3) +
  annotate("segment", x = c(0.67, 0.97, 1.27), xend = c(0.73, 1.03, 1.33),
           y = 1.17, yend = 1.17, colour = "#F1FAEE", linewidth = 0.3,
           arrow = grid::arrow(length = grid::unit(0.03, "inches"), type = "closed")) +
  scale_colour_manual(values = cols2, guide = "none") +
  name("#E9C46A")
save(g2, "option2_grid.png")

# 3. Explosion burst with a phylogeny inside
burst_a <- seq(0, 2 * pi, length.out = 25)[-25]
burst <- data.frame(x = 1 + ifelse(seq_along(burst_a) %% 2 == 1, 0.55, 0.36) * cos(burst_a),
                    y = 1.18 + ifelse(seq_along(burst_a) %% 2 == 1, 0.55, 0.36) * sin(burst_a))
ts3 <- tree_segments(ape::rcoal(10), c(0.76, 1.24, 0.98, 1.38))
g3 <- base("#FFB703", "#023047") +
  geom_polygon(aes(x, y), burst, fill = "#FB8500") +
  geom_segment(aes(x, y, xend = xend, yend = yend), ts3$seg, colour = "#023047",
               linewidth = 0.5, lineend = "square") +
  name("#023047")
save(g3, "option3_burst.png")

# 4. Tree whose tips sort into patches along an environmental gradient
ts4 <- tree_segments(ape::rcoal(12), c(0.62, 1.38, 1.25, 1.6))
patch <- data.frame(x = seq(0.6, 1.4, length.out = 5), y = 0.95)
lines4 <- data.frame(x = ts4$tips$x, y = 1.25, xend = patch$x[findInterval(ts4$tips$x,
                     seq(0.5, 1.5, length.out = 6), all.inside = TRUE)], yend = 1.02)
g4 <- base("#F1FAEE", "#1D3557") +
  geom_segment(aes(x, y, xend = xend, yend = yend), lines4, colour = "#A8DADC", linewidth = 0.3) +
  geom_segment(aes(x, y, xend = xend, yend = yend), ts4$seg, colour = "#1D3557",
               linewidth = 0.45, lineend = "square") +
  geom_tile(aes(x, y), patch, width = 0.15, height = 0.13, fill = "#457B9D") +
  annotate("segment", x = 0.55, xend = 1.45, y = 0.8, yend = 0.8, colour = "#1D3557",
           linewidth = 0.35, arrow = grid::arrow(length = grid::unit(0.04, "inches"), type = "closed")) +
  annotate("text", 1, 0.72, label = "environment", colour = "#1D3557", size = 1.6) +
  name("#1D3557", y = 0.5)
save(g4, "option4_sorting.png")
