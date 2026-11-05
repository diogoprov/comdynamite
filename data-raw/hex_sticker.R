# Hex sticker for comdynamite: a bundle of dynamite whose fuse is a phylogeny,
# with one coloured spark per species at the tips. Run from the package root;
# writes man/figures/logo.png. Needs ape, ggplot2, hexSticker, showtext and
# sysfonts (the font comes from Google Fonts, so it needs a connection).
library(ggplot2)
library(hexSticker)
library(showtext)

sysfonts::font_add_google("Righteous", "righteous")
showtext_auto()

ink <- "#1D3557"
paper <- "#F1FAEE"
stick <- "#D62828"
band <- "#1D3557"
sparks <- c("#F77F00", "#FCBF49", "#2A9D8F", "#E76F51")

# segments of a rectangular tree drawn upwards, scaled to box (x0, x1, y0, y1)
tree_segments <- function(tr, box) {
  tr <- ape::reorder.phylo(tr, "cladewise")
  n <- ape::Ntip(tr)
  e <- tr$edge
  depth <- ape::node.depth.edgelength(tr)
  x <- numeric(n + tr$Nnode)
  tips <- e[e[, 2] <= n, 2]
  x[tips] <- seq_along(tips)
  for (v in unique(ape::reorder.phylo(tr, "postorder")$edge[, 1]))
    x[v] <- mean(range(x[e[e[, 1] == v, 2]]))
  X <- box[1] + (x - 1) / (n - 1) * (box[2] - box[1])
  Y <- box[3] + depth / max(depth) * (box[4] - box[3])
  p <- e[, 1]
  ch <- e[, 2]
  stems <- data.frame(x = X[ch], xend = X[ch], y = Y[p], yend = Y[ch])
  bars <- do.call(rbind, lapply(unique(p), function(v) {
    k <- ch[p == v]
    data.frame(x = min(X[k]), xend = max(X[k]), y = Y[v], yend = Y[v])
  }))
  root <- n + 1
  list(segments = rbind(stems, bars),
       tips = data.frame(x = X[seq_len(n)], y = Y[seq_len(n)]),
       root = c(x = X[root], y = Y[root]))
}

phy <- ape::compute.brlen(ape::read.tree(text = "(((a,b),c),((d,e),(f,(g,h))));"),
                          method = "Grafen")
tree <- tree_segments(phy, c(0.5, 3.5, 2.0, 3.55))
# centre the root over the bundle
shift <- 2 - tree$root[["x"]]
tree$segments[c("x", "xend")] <- tree$segments[c("x", "xend")] + shift
tree$tips$x <- tree$tips$x + shift
tree$root[["x"]] <- 2
tree$tips$colour <- rep_len(sparks, nrow(tree$tips))

# three sticks, centred under the root of the tree
sticks <- data.frame(xmin = tree$root[["x"]] + c(-1.02, -0.32, 0.38),
                     ymin = 0, ymax = 1.55)
sticks$xmax <- sticks$xmin + 0.64
fuse <- data.frame(x = tree$root[["x"]], xend = tree$root[["x"]],
                   y = 1.55, yend = tree$root[["y"]])

p <- ggplot() +
  geom_segment(data = rbind(tree$segments, fuse),
               aes(x, y, xend = xend, yend = yend),
               colour = paper, linewidth = 0.55, lineend = "square") +
  geom_point(data = tree$tips, aes(x, y, fill = colour), shape = 21,
             colour = paper, size = 2.1, stroke = 0.3) +
  geom_rect(data = sticks, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = stick, colour = paper, linewidth = 0.3) +
  annotate("rect", xmin = min(sticks$xmin) - 0.06, xmax = max(sticks$xmax) + 0.06,
           ymin = 0.62, ymax = 0.92, fill = paper, colour = NA) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 4), ylim = c(-0.1, 3.8), expand = FALSE) +
  theme_void() +
  theme_transparent()

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
sticker(p, package = "comdynamite",
        p_family = "righteous", p_size = 15, p_color = paper, p_y = 0.56,
        s_x = 1, s_y = 1.17, s_width = 1.05, s_height = 0.95,
        h_fill = ink, h_color = stick, h_size = 1.6,
        url = "github.com/diogoprov/comdynamite", u_size = 3.3, u_color = paper,
        filename = "man/figures/logo.png", dpi = 300)
