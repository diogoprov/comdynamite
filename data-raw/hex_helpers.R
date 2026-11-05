# Shared helpers for the hex sticker scripts in data-raw/ (hexSticker geometry:
# hexagon of radius 1 centred at (1, 1), 43.9 x 50.8 mm).
library(ggplot2)
if (!exists("FRAME")) FRAME <- TRUE

a <- pi / 2 + seq(0, 2 * pi, length.out = 7)[-7]
hex <- data.frame(x = 1 + cos(a), y = 1 + sin(a))

# segments of a rectangular tree drawn upwards in box (x0, x1, y0, y1)
tree_segments <- function(tr, box) {
  tr <- ape::reorder.phylo(tr, "cladewise")
  n <- ape::Ntip(tr); e <- tr$edge
  depth <- ape::node.depth.edgelength(tr)
  x <- numeric(n + tr$Nnode)
  tips <- e[e[, 2] <= n, 2]; x[tips] <- seq_along(tips)
  for (v in unique(ape::reorder.phylo(tr, "postorder")$edge[, 1]))
    x[v] <- mean(range(x[e[e[, 1] == v, 2]]))
  X <- box[1] + (x - 1) / (n - 1) * (box[2] - box[1])
  Y <- box[3] + depth / max(depth) * (box[4] - box[3])
  p <- e[, 1]; ch <- e[, 2]
  stems <- data.frame(x = X[ch], xend = X[ch], y = Y[p], yend = Y[ch])
  bars <- do.call(rbind, lapply(unique(p), function(v) {
    k <- ch[p == v]; data.frame(x = min(X[k]), xend = max(X[k]), y = Y[v], yend = Y[v])
  }))
  list(seg = rbind(stems, bars), tips = data.frame(x = X[tips], y = box[4]))
}

base <- function(fill, border) {
  ggplot() +
    (if (FRAME) geom_polygon(aes(x, y), data = hex, fill = fill, colour = border,
                             linewidth = 1.6)) +
    coord_fixed(xlim = c(0.12, 1.88), ylim = c(-0.02, 2.02), expand = FALSE) +
    theme_void()
}
name <- function(col, y = 0.47, size = 3.9)
  if (FRAME) annotate("text", 1, y, label = "comdynamite", colour = col, size = size,
                      fontface = "bold")
save <- function(g, f) if (!is.null(out)) ggsave(file.path(out, f), g, width = 43.9, height = 50.8,
                              units = "mm", dpi = 600, bg = "transparent")

