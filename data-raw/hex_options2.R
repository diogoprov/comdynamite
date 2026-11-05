# Second round of hex sticker options for comdynamite (solid colours), combining
# the metacommunity grid (option 2) and the sorting of clades along an
# environmental gradient (option 4) of data-raw/hex_options.R.
# Usage: Rscript data-raw/hex_options2.R out_dir
source("data-raw/hex_helpers.R")
out <- if (FRAME) commandArgs(TRUE)[1] else NULL
set.seed(11)

# a tree with three clades of four tips, drawn upwards in `box`
clade_tree <- function(box) {
  tr <- ape::read.tree(text = paste0("((((a:1,b:1):1,(c:1,d:1):1):1.5,((e:1,f:1):1,(g:1,h:1):1):1.5):1,",
                                     "((i:1,j:1):1,(k:1,l:1):1):2.5);"))
  ts <- tree_segments(tr, box)
  ts$tips$clade <- factor(rep(1:3, each = 4))
  ts
}
# patches (sites) on a grid, with an environment that increases left to right
patches <- function(nx, ny, x0, x1, y0, y1) {
  g <- expand.grid(i = seq_len(nx), j = seq_len(ny))
  g$x <- x0 + (g$i - 1) / (nx - 1) * (x1 - x0)
  g$y <- y0 + (g$j - 1) / max(ny - 1, 1) * (y1 - y0)
  g$env <- factor(g$i)
  g
}
# species (dots) in each patch, mostly from the clade that matches the patch
dots_in <- function(g, nclade = 3, dx = 0.045) {
  do.call(rbind, lapply(seq_len(nrow(g)), function(s) {
    home <- ceiling(g$i[s] / max(g$i) * nclade)
    k <- c(home, if (stats::runif(1) < 0.4) sample(setdiff(seq_len(nclade), home), 1))
    data.frame(x = g$x[s] + c(-dx, dx)[seq_along(k)] * (length(k) > 1),
               y = g$y[s], clade = factor(k, seq_len(nclade)))
  }))
}
arrow_seg <- function(...) annotate("segment", ..., linewidth = 0.3,
  arrow = grid::arrow(length = grid::unit(0.025, "inches"), type = "closed"))

clade_cols <- c("#E76F51", "#2A9D8F", "#E9C46A")

# A. Clades sort into a grid of patches along the gradient
tsA <- clade_tree(c(0.6, 1.4, 1.3, 1.62))
gA <- patches(3, 3, 0.66, 1.34, 0.78, 1.14)
gA_dots <- dots_in(gA)
pA <- base("#264653", "#E9C46A") +
  geom_tile(aes(x, y, fill = env), gA, width = 0.26, height = 0.15) +
  scale_fill_manual(values = c("#2F5D6E", "#3B7084", "#4A869C"), guide = "none") +
  geom_segment(aes(x, y, xend = xend, yend = yend), tsA$seg, colour = "#F1FAEE",
               linewidth = 0.45, lineend = "square") +
  geom_point(aes(x, y, colour = clade), tsA$tips, size = 1.1) +
  geom_point(aes(x, y, colour = clade), gA_dots, size = 1.3) +
  scale_colour_manual(values = clade_cols, guide = "none") +
  name("#E9C46A", y = 0.5)
save(pA, "optionA_clades_grid.png")

# B. Gradient grid with dispersal arrows, no tree
gB <- patches(4, 3, 0.58, 1.42, 0.86, 1.44)
gB_dots <- dots_in(gB)
pB <- base("#F1FAEE", "#1D3557") +
  geom_tile(aes(x, y, fill = env), gB, width = 0.24, height = 0.25) +
  scale_fill_manual(values = c("#A8DADC", "#7FB8C9", "#5795B0", "#457B9D"), guide = "none") +
  geom_point(aes(x, y, colour = clade), gB_dots, size = 1.4) +
  scale_colour_manual(values = c("#E63946", "#1D3557", "#F4A261"), guide = "none") +
  arrow_seg(x = c(0.685, 0.965, 1.245), xend = c(0.755, 1.035, 1.315), y = 1.15, yend = 1.15,
            colour = "#1D3557") +
  arrow_seg(x = 1.14, xend = 1.14, y = 1.255, yend = 1.335, colour = "#1D3557") +
  name("#1D3557", y = 0.52)
save(pB, "optionB_gradient_grid.png")

# C. Honeycomb of patches fed by the tips of a tree
hex_patch <- function(cx, cy, r = 0.12) {
  a <- pi / 2 + seq(0, 2 * pi, length.out = 7)[-7]
  data.frame(x = cx + r * cos(a), y = cy + r * sin(a))
}
centres <- data.frame(cx = c(0.79, 1, 1.21, 0.895, 1.105),
                      cy = c(0.98, 0.98, 0.98, 0.8, 0.8))
centres$env <- factor(c(1, 2, 3, 1, 3))
hexes <- do.call(rbind, lapply(seq_len(nrow(centres)), function(k)
  cbind(hex_patch(centres$cx[k], centres$cy[k]), id = k, env = centres$env[k])))
tsC <- clade_tree(c(0.64, 1.36, 1.22, 1.6))
feed <- data.frame(x = tsC$tips$x, y = 1.22,
                   xend = centres$cx[as.integer(tsC$tips$clade)], yend = 1.08)
pC <- base("#1D3557", "#F4A261") +
  geom_segment(aes(x, y, xend = xend, yend = yend), feed, colour = "#457B9D", linewidth = 0.3) +
  geom_polygon(aes(x, y, group = id, fill = env), hexes, colour = "#1D3557", linewidth = 0.4) +
  scale_fill_manual(values = clade_cols, guide = "none") +
  geom_segment(aes(x, y, xend = xend, yend = yend), tsC$seg, colour = "#F1FAEE",
               linewidth = 0.45, lineend = "square") +
  name("#F4A261", y = 0.5)
save(pC, "optionC_honeycomb.png")

# D. Tips drop into a row of patches along the gradient, with dispersal
tsD <- clade_tree(c(0.6, 1.4, 1.25, 1.6))
gD <- patches(3, 1, 0.68, 1.32, 0.98, 0.98)
dropD <- data.frame(x = tsD$tips$x, y = 1.25, xend = gD$x[as.integer(tsD$tips$clade)],
                    yend = 1.07)
gD_dots <- data.frame(x = rep(gD$x, each = 3) + c(-0.06, 0, 0.06), y = 0.96,
                      clade = factor(c(1, 1, 2, 2, 2, 3, 3, 3, 1)))
pD <- base("#2A9D8F", "#264653") +
  geom_segment(aes(x, y, xend = xend, yend = yend), dropD, colour = "#A8DADC", linewidth = 0.3) +
  geom_tile(aes(x, y), gD, width = 0.24, height = 0.17, fill = "#264653") +
  geom_point(aes(x, y, colour = clade), gD_dots, size = 1.3) +
  scale_colour_manual(values = c("#E76F51", "#F1FAEE", "#E9C46A"), guide = "none") +
  geom_segment(aes(x, y, xend = xend, yend = yend), tsD$seg, colour = "#264653",
               linewidth = 0.5, lineend = "square") +
  annotate("curve", x = 0.76, xend = 0.92, y = 0.86, yend = 0.86, curvature = 0.5,
           colour = "#F1FAEE", linewidth = 0.3,
           arrow = grid::arrow(length = grid::unit(0.025, "inches"), type = "closed")) +
  annotate("curve", x = 1.08, xend = 1.24, y = 0.86, yend = 0.86, curvature = 0.5,
           colour = "#F1FAEE", linewidth = 0.3,
           arrow = grid::arrow(length = grid::unit(0.025, "inches"), type = "closed")) +
  arrow_seg(x = 0.62, xend = 1.38, y = 0.72, yend = 0.72, colour = "#264653") +
  annotate("text", 1, 0.66, label = "environment", colour = "#264653", size = 1.5) +
  name("#F1FAEE", y = 0.5)
save(pD, "optionD_sorting_dispersal.png")

# E. Option A with the dispersal arrows of option D and a clearer environmental
#    gradient across the patches (one solid colour per column)
tsE <- clade_tree(c(0.6, 1.4, 1.36, 1.64))
gE <- patches(3, 3, 0.66, 1.34, 0.84, 1.16)
set.seed(12)
gE_dots <- dots_in(gE, dx = 0.04)
curved <- function(x, xend, y, yend, curvature = -0.5)
  annotate("curve", x = x, xend = xend, y = y, yend = yend, curvature = curvature,
           colour = "#F1FAEE", linewidth = 0.3,
           arrow = grid::arrow(length = grid::unit(0.022, "inches"), type = "closed"))
pE <- base("#264653", "#E9C46A") +
  geom_tile(aes(x, y, fill = env), gE, width = 0.22, height = 0.13) +
  scale_fill_manual(values = c("#3A6475", "#4F8AA3", "#72B3CB"), guide = "none") +
  geom_segment(aes(x, y, xend = xend, yend = yend), tsE$seg, colour = "#F1FAEE",
               linewidth = 0.45, lineend = "square") +
  geom_point(aes(x, y, colour = clade), tsE$tips, size = 1.1) +
  geom_point(aes(x, y, colour = clade), gE_dots, size = 1.25) +
  scale_colour_manual(values = c("#E76F51", "#F1FAEE", "#E9C46A"), guide = "none") +
  curved(0.74, 0.92, 1.235, 1.235) +
  curved(1.08, 1.26, 1.235, 1.235) +
  curved(1.26, 1.08, 0.765, 0.765) +
  curved(0.92, 0.74, 0.765, 0.765) +
  arrow_seg(x = 0.6, xend = 1.4, y = 0.68, yend = 0.68, colour = "#E9C46A") +
  annotate("text", 1, 0.635, label = "environment", colour = "#E9C46A", size = 1.4) +
  name("#E9C46A", y = 0.47)
save(pE, "optionE_clades_dispersal.png")
