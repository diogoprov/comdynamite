# Third round of hex sticker options for comdynamite: option E (the current logo)
# with the patches redrawn as dynamite and the dispersal arrows as lit fuses.
# Usage: Rscript data-raw/hex_options3.R out_dir
source("data-raw/hex_helpers.R")
out <- if (FRAME) commandArgs(TRUE)[1] else NULL

ink <- "#264653"; gold <- "#E9C46A"; cream <- "#F1FAEE"
clade_cols <- c("#E76F51", cream, gold)
reds <- c("#9B2F27", "#C43D2B", "#E8582F")   # environment, low to high
caps <- c("#D9907F", "#EBA48F", "#F6B9A1")   # cut end of the stick in each column

# species in the nine communities of the current logo (rows top to bottom)
comm <- list(c(1, 2), c(2, 1), c(3, 1),
             c(1, 2), 2,       c(3, 2),
             c(1, 3), 2,       c(3, 1))

# the tree of option E: three clades of four tips, drawn upwards in `box`
clade_tree <- function(box) {
  tr <- ape::read.tree(text = paste0("((((a:1,b:1):1,(c:1,d:1):1):1.5,((e:1,f:1):1,(g:1,h:1):1):1.5):1,",
                                     "((i:1,j:1):1,(k:1,l:1):1):2.5);"))
  ts <- tree_segments(tr, box)
  ts$tips$clade <- factor(rep(1:3, each = 4))
  ts
}
tree_layers <- function(ts) list(
  geom_segment(aes(x, y, xend = xend, yend = yend), ts$seg, colour = cream,
               linewidth = 0.45, lineend = "square"),
  geom_point(aes(x, y, colour = clade), ts$tips, size = 1.1),
  scale_colour_manual(values = clade_cols, guide = "none"),
  scale_fill_manual(values = clade_cols, guide = "none"))

# polygon in stick coordinates (along, across), turned by `angle` and moved to (cx, cy)
place <- function(u, v, cx, cy, angle) {
  th <- angle * pi / 180
  list(x = cx + u * cos(th) - v * sin(th), y = cy + u * sin(th) + v * cos(th))
}
poly <- function(p, fill) annotate("polygon", x = p$x, y = p$y, fill = fill, colour = NA)

# one stick of dynamite: rounded far end, oval cut end at the `angle` side, and
# paper bands given as a list of c(from, to) along the stick (-0.5 to 0.5)
stick <- function(cx, cy, len, thick, body, cap, bands = list(), band = ink,
                  angle = 0) {
  r <- thick / 2; rx <- thick * 0.2; h <- len / 2
  t1 <- seq(-pi / 2, pi / 2, length.out = 24)
  outline <- place(c(h + rx * cos(t1), -h - rx * cos(rev(t1))),
                   c(r * sin(t1), r * sin(rev(t1))), cx, cy, angle)
  up <- if (cos(angle * pi / 180) < -0.5) -1 else 1
  shine <- place(c(-h, h, h, -h), up * r * c(0.45, 0.45, 0.72, 0.72), cx, cy, angle)
  t2 <- seq(0, 2 * pi, length.out = 40)
  end <- place(h + rx * cos(t2), r * sin(t2), cx, cy, angle)
  c(list(poly(outline, body), poly(shine, cap)),
    unlist(lapply(bands, function(b) {
      b <- b * len
      strip <- function(lo, hi, fill)
        poly(place(c(hi + rx * cos(t1), lo + rx * cos(rev(t1))),
                   c(r * sin(t1), r * sin(rev(t1))), cx, cy, angle), fill)
      list(strip(b[1] - 0.012, b[2] + 0.012, cream), strip(b[1], b[2], band))
    }), recursive = FALSE),
    list(poly(end, cap)))
}

# dispersal as a fuse: a bent path from p0 to p1 that ends in an arrow and a spark
fuse <- function(p0, p1, bend = 0.2, colour = cream, size = 0.3, spark_r = 0.032,
                 head = TRUE) {
  t <- seq(0, 1, length.out = 60)
  d <- p1 - p0
  ctrl <- (p0 + p1) / 2 + bend * c(-d[2], d[1])
  x <- (1 - t)^2 * p0[1] + 2 * t * (1 - t) * ctrl[1] + t^2 * p1[1]
  y <- (1 - t)^2 * p0[2] + 2 * t * (1 - t) * ctrl[2] + t^2 * p1[2]
  c(list(geom_path(aes(x, y), data.frame(x = x, y = y), colour = colour,
                   linewidth = size, lineend = "round",
                   arrow = if (head) grid::arrow(length = grid::unit(0.024, "inches"),
                                                 type = "closed"))),
    if (spark_r > 0) spark(p1[1], p1[2], spark_r))
}
spark <- function(x, y, r) {
  star <- function(r, k = 0.42, n = 8) {
    a <- pi / 2 + seq(0, 2 * pi, length.out = 2 * n + 1)[-(2 * n + 1)]
    rr <- rep(c(r, r * k), n)
    list(x = x + rr * cos(a), y = y + rr * sin(a))
  }
  list(poly(star(r), "#F4A261"), poly(star(r * 0.55), gold))
}
gradient_axis <- function(y = 0.68) list(
  annotate("segment", x = 0.6, xend = 1.4, y = y, yend = y, colour = gold, linewidth = 0.3,
           arrow = grid::arrow(length = grid::unit(0.025, "inches"), type = "closed")),
  annotate("text", 1, y - 0.045, label = "environment", colour = gold, size = 1.4))
dots <- function(x, y, clade, size = 1.25)
  geom_point(aes(x, y, fill = clade),
             data.frame(x = x, y = y, clade = factor(clade, 1:3)), size = size,
             shape = 21, colour = ink, stroke = 0.35)

cols <- c(0.66, 1, 1.34)
rows <- c(1.16, 1, 0.84)

# F. Each patch becomes a stick, in the same 3 x 3 grid; the band holds the species
tsF <- clade_tree(c(0.6, 1.4, 1.36, 1.64))
gF <- expand.grid(i = 1:3, j = 1:3)
pF <- base(ink, gold)
for (s in seq_len(nrow(gF))) {
  i <- gF$i[s]
  pF <- pF + stick(cols[i], rows[gF$j[s]], 0.21, 0.12, reds[i], caps[i],
                   bands = list(c(-0.42, -0.4), c(0.34, 0.36)), band = cream, angle = if (gF$j[s] == 3) 180 else 0)
}
dF <- do.call(rbind, lapply(seq_len(nrow(gF)), function(s) {
  k <- comm[[s]]
  data.frame(x = cols[gF$i[s]] + if (length(k) > 1) c(-0.036, 0.036) else 0,
             y = rows[gF$j[s]], clade = k)
}))
pF <- pF + tree_layers(tsF) + dots(dF$x, dF$y, dF$clade, 1.25) +
  fuse(c(0.775, 1.17), c(0.9, 1.235), 0.55) + fuse(c(1.115, 1.17), c(1.24, 1.235), 0.55) +
  fuse(c(1.225, 0.83), c(1.1, 0.765), 0.55) + fuse(c(0.885, 0.83), c(0.76, 0.765), 0.55) +
  gradient_axis() + name(gold, y = 0.47)
save(pF, "optionF_sticks_grid.png")

# G. One bundle per position on the gradient: three upright sticks tied with a
#    strap, a fuse from the top of each bundle to the next
tsG <- clade_tree(c(0.6, 1.4, 1.38, 1.64))
pG <- base(ink, gold)
dG <- NULL
for (i in 1:3) for (j in 1:3) {
  x <- cols[i] + (j - 2) * 0.083
  pG <- pG + stick(x, 0.965, 0.3, 0.08, reds[i], caps[i], angle = 90)
  k <- comm[[(j - 1) * 3 + i]]
  dG <- rbind(dG, data.frame(x = x, y = c(1.06, 0.87)[seq_along(k)], clade = k))
}
for (i in 1:3)
  pG <- pG + annotate("rect", xmin = cols[i] - 0.133, xmax = cols[i] + 0.133,
                      ymin = 0.945, ymax = 0.985, fill = cream)
pG <- pG + tree_layers(tsG) + dots(dG$x, dG$y, dG$clade, 1.1) +
  fuse(c(0.66, 1.135), c(0.96, 1.15), 0.5) + fuse(c(1, 1.135), c(1.3, 1.15), 0.5) +
  gradient_axis(0.7) + name(gold, y = 0.47)
save(pG, "optionG_bundles.png")

# H. The root of the tree is the fuse: it lights the middle stick and dispersal
#    carries the fire to the sticks at both ends of the gradient
tsH <- clade_tree(c(0.6, 1.4, 1.4, 1.64))
rootH <- subset(tsH$seg, y == min(y) & y == yend)
rx <- (rootH$x + rootH$xend) / 2
pH <- base(ink, gold)
dH <- NULL
for (i in 1:3) {
  pH <- pH + stick(cols[i], 0.96, 0.32, 0.17, reds[i], caps[i],
                   bands = list(c(-0.4, -0.38), c(0.27, 0.29)), band = cream, angle = 90)
  k <- unique(unlist(comm[c(i, i + 3, i + 6)]))
  dH <- rbind(dH, data.frame(x = cols[i], y = 1.02 - (seq_along(k) - 1) * 0.07, clade = k))
}
pH <- pH + tree_layers(tsH) + dots(dH$x, dH$y, dH$clade, 1.6) +
  fuse(c(rx, 1.4), c(1, 1.15), 0.12, spark_r = 0, head = FALSE, size = 0.45) +
  fuse(c(1, 1.15), c(0.68, 1.16), -0.4) + fuse(c(1, 1.15), c(1.32, 1.16), 0.4) +
  gradient_axis(0.7) + name(gold, y = 0.47)
save(pH, "optionH_tree_fuse.png")

# I. Three long sticks strapped into one horizontal bundle: the columns of the
#    grid are the stretches between the straps, and the environment runs along it
tsI <- clade_tree(c(0.6, 1.4, 1.36, 1.64))
pI <- base(ink, gold)
for (y in rows) {
  pI <- pI + stick(0.97, y, 0.92, 0.15, reds[1], caps[1]) +
    stick(1.215, y, 0.43, 0.15, reds[3], caps[3]) +
    annotate("rect", xmin = 0.83, xmax = 1.17, ymin = y - 0.075, ymax = y + 0.075,
             fill = reds[2]) +
    annotate("rect", xmin = 0.83, xmax = 1.17, ymin = y + 0.034, ymax = y + 0.054,
             fill = caps[2])
}
pI <- pI + annotate("rect", xmin = c(0.815, 1.155), xmax = c(0.845, 1.185),
                    ymin = 0.77, ymax = 1.23, fill = cream)
dI <- do.call(rbind, lapply(seq_len(nrow(gF)), function(s) {
  k <- comm[[s]]
  data.frame(x = cols[gF$i[s]] + if (length(k) > 1) c(-0.04, 0.04) else 0,
             y = rows[gF$j[s]], clade = k)
}))
pI <- pI + tree_layers(tsI) + dots(dI$x, dI$y, dI$clade, 1.25) +
  fuse(c(0.72, 1.235), c(0.94, 1.24), 0.4) + fuse(c(1.06, 1.235), c(1.28, 1.24), 0.4) +
  fuse(c(1.28, 0.765), c(1.06, 0.76), 0.4) + fuse(c(0.94, 0.765), c(0.72, 0.76), 0.4) +
  gradient_axis(0.665) + name(gold, y = 0.47)
save(pI, "optionI_long_bundle.png")

# Refinements of F: clades as symbols (circle, triangle, star) at the tips and in
# the communities, sticks with more volume, and a short wick on every stick.
# environment from low to high: wine, red, orange
bodies <- c("#7D1A2E", "#CF3A27", "#F58A1F")
ends <- c("#C9808A", "#EFA08C", "#FBCB8E")
shades <- c("#560F1F", "#9C2519", "#C2620F")

# symbols as polygons in hexagon coordinates, so they scale with the sticks
symbols_df <- function(x, y, clade, r, fill) {
  do.call(rbind, lapply(seq_along(x), function(s) {
    k <- clade[s]
    if (k == 1) {
      a <- seq(0, 2 * pi, length.out = 30)[-30]; rr <- r * 0.95
    } else if (k == 2) {
      a <- pi / 2 + c(0, 2, 4) * pi / 3; rr <- r * 1.3
    } else {
      a <- pi / 2 + seq(0, 2 * pi, length.out = 11)[-11]; rr <- rep(r * c(1.45, 0.62), 5)
    }
    data.frame(x = x[s] + rr * cos(a), y = y[s] + rr * sin(a), id = s, fill = fill[k])
  }))
}
symbols_layer <- function(d, outline = NA, stroke = 0.2)
  geom_polygon(aes(x, y, group = id, fill = fill), d, colour = outline, linewidth = stroke)

wick <- function(p0, p1, bend = 0.5, size = 0.45)
  fuse(p0, p1, bend, size = size, spark_r = 0, head = FALSE)

# one variant of F; `label` is the colour of the paper label (NA for none)
option_F <- function(sym_stick, sym_tip = sym_stick, label = NA, outline = NA) {
  ts <- clade_tree(c(0.56, 1.44, 1.38, 1.64))
  cx <- c(0.62, 1, 1.38); cy <- c(1.14, 0.955, 0.77)
  L <- 0.265; Th <- 0.155; lab <- c(-0.3, 0.24)
  g <- expand.grid(i = 1:3, j = 1:3)
  p <- base(ink, gold)
  for (s in seq_len(nrow(g))) {
    i <- g$i[s]; y <- cy[g$j[s]]
    p <- p + stick(cx[i], y, L, Th, bodies[i], ends[i], angle = if (g$j[s] == 3) 180 else 0,
                   bands = if (is.na(label)) list(c(-0.52, -0.48), c(0.37, 0.41))
                           else list(lab),
                   band = if (is.na(label)) cream else label)
    # shadow along the underside
    p <- p + annotate("rect", xmin = cx[i] - L / 2, xmax = cx[i] + L / 2,
                      ymin = y - Th / 2, ymax = y - Th * 0.3, fill = shades[i], alpha = 0.55)
  }
  d <- do.call(rbind, lapply(seq_len(nrow(g)), function(s) {
    k <- comm[[s]]
    mid <- mean(lab) * L * if (g$j[s] == 3) -1 else 1   # the label flips with the stick
    data.frame(x = cx[g$i[s]] + mid + if (length(k) > 1) c(-0.04, 0.04) else 0,
               y = cy[g$j[s]], clade = k)
  }))
  e <- L / 2 + 0.015                     # from the centre of a stick to its cut end
  top <- cy[1] + Th / 2 + 0.012; low <- cy[3] - Th / 2 - 0.012
  p + scale_fill_identity() +
    geom_segment(aes(x, y, xend = xend, yend = yend), ts$seg, colour = cream,
                 linewidth = 0.45, lineend = "square") +
    symbols_layer(symbols_df(ts$tips$x, ts$tips$y, as.integer(ts$tips$clade), 0.026, sym_tip)) +
    symbols_layer(symbols_df(d$x, d$y, d$clade, 0.028, sym_stick), outline) +
    # short wicks on the sticks that do not light a neighbour
    wick(c(cx[3] + e, cy[1]), c(cx[3] + e + 0.045, cy[1] + 0.04)) +
    wick(c(cx[1] + e, cy[2]), c(cx[1] + e + 0.04, cy[2] + 0.045)) +
    wick(c(cx[2] + e, cy[2]), c(cx[2] + e + 0.04, cy[2] + 0.045)) +
    wick(c(cx[3] + e, cy[2]), c(cx[3] + e + 0.045, cy[2] + 0.04)) +
    wick(c(cx[1] - e, cy[3]), c(cx[1] - e - 0.045, cy[3] - 0.04)) +
    fuse(c(cx[1] + e, cy[1] + 0.01), c(cx[2] - 0.1, top), 0.55, size = 0.5, spark_r = 0.058, head = FALSE) +
    fuse(c(cx[2] + e, cy[1] + 0.01), c(cx[3] - 0.1, top), 0.55, size = 0.5, spark_r = 0.058, head = FALSE) +
    fuse(c(cx[3] - e, cy[3] - 0.01), c(cx[2] + 0.1, low), 0.55, size = 0.5, spark_r = 0.058, head = FALSE) +
    fuse(c(cx[2] - e, cy[3] - 0.01), c(cx[1] + 0.1, low), 0.55, size = 0.5, spark_r = 0.058, head = FALSE) +
    name(gold, y = 0.47)
}
vivid <- c("#E76F51", "#35B8A6", "#F2B233")
dark <- rep("#1F3A45", 3)
pF1 <- option_F(vivid, label = "#F6EBD2", outline = ink)
pF2 <- option_F(rep(gold, 3), outline = ink)
pF3 <- option_F(dark, sym_tip = rep(gold, 3), label = gold)
save(pF1, "optionF1_paper_label.png")
save(pF2, "optionF2_gold_symbols.png")
save(pF3, "optionF3_gold_label.png")
