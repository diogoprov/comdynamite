# Figures of the comdynamite preprint (eco-evolutionary engine, version 0.9.0). Run from the
# package root:
#   Rscript paper/figures.R
# Figure 1 is drawn here. Figures 3 and 4 assemble panels that the scripts of the articles of the
# package website draw (data-raw/figures_leibold2019.R, figures_urban2009.R, figures_hajian2026.R
# from the saved replications, and animation_quantgen.R, which reruns its single simulations with
# fixed seeds); here they are drawn again at 300 dpi into paper/figures/panels/. Figure 2 is the
# nesting map of data-raw/figure_nesting_map.R, included directly from vignettes/articles/.
library(ggplot2)
library(magick)

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7", "#999999")
dir.create("paper/figures", showWarnings = FALSE)

# Figure 1: structure of the eco-evolutionary engine ----------------------------------------
boxes <- data.frame(
  x = c(1.4, 4.2, 7.0, 1.4, 4.2, 7.0),
  y = c(3.9, 3.9, 3.9, 1.3, 1.3, 1.3),
  title = c("Landscape and time", "Ecological dynamics", "Quantitative genetics",
            "Inference", "Outputs", "Speciation"),
  body = c("sites, kernels, env axes\nenv_over_time(), env_pulse()\ndispersal_time, kernels",
           "growth by niche match\ncompetition, dispersal\ndemographic drift\n(Thompson et al. 2020)",
           "genetic mean per population\nselection, drift (Ne), gene flow\nfounders, plasticity, G or G(t)\n(Lande 1976, 1979)",
           "abc_choose(), abc_fit()\nabc_recovery(), abc_ppc()",
           "abundances, genetic means\nfull and reconstructed trees\nhistory of ranges\nmaladaptation(), persistence()",
           "point, fission, geographic,\necological modes\nprotracted completion\n(Hubert et al. 2015)"),
  fill = c("land", "eco", "gen", "inf", "out", "spe")
)
arr <- data.frame(x = c(2.6, 5.4, 7.0, 5.8, 3.0),
                  y = c(3.9, 4.0, 3.3, 1.3, 1.3),
                  xend = c(3.0, 5.8, 7.0, 5.4, 2.6),
                  yend = c(3.9, 4.0, 1.9, 1.3, 1.3))
fig1 <- ggplot() +
  geom_segment(data = arr, aes(x, y, xend = xend, yend = yend),
               arrow = arrow(length = unit(2.2, "mm"), type = "closed"), colour = "grey35") +
  geom_segment(aes(x = 5.4, y = 3.8, xend = 5.8, yend = 3.8),
               arrow = arrow(length = unit(2.2, "mm"), type = "closed", ends = "first"),
               colour = "grey35") +
  geom_curve(aes(x = 4.2, y = 3.3, xend = 4.2, yend = 1.9), curvature = 0,
             arrow = arrow(length = unit(2.2, "mm"), type = "closed"), colour = "grey35") +
  geom_tile(data = boxes, aes(x, y, fill = fill), width = 2.4, height = 1.2,
            colour = "grey30", alpha = 0.35) +
  geom_text(data = boxes, aes(x, y + 0.45, label = title), fontface = "bold", size = 3.3) +
  geom_text(data = boxes, aes(x, y - 0.12, label = body), size = 2.4, lineheight = 0.95) +
  annotate("text", x = 5.6, y = 4.75, size = 2.3, fontface = "italic", colour = okabe[4],
           label = "phenotype vs environment;\nimmigrants, founders, Ne") +
  annotate("text", x = 7.12, y = 2.6, size = 2.3, fontface = "italic", colour = okabe[4],
           hjust = 0, label = "divergence\nof populations") +
  annotate("text", x = 4.32, y = 2.6, size = 2.3, fontface = "italic", colour = okabe[4],
           hjust = 0, label = "records") +
  scale_fill_manual(values = c(land = okabe[7], eco = okabe[2], gen = okabe[3],
                               spe = okabe[1], out = okabe[5], inf = okabe[6]), guide = "none") +
  coord_cartesian(xlim = c(-0.05, 8.45), ylim = c(0.5, 5.0), expand = FALSE) +
  theme_void()
ggsave("paper/figures/fig1_model.png", fig1, width = 7.5, height = 4.3, dpi = 300, bg = "white")

# Panels drawn again at 300 dpi by the article scripts --------------------------------------
panels <- "paper/figures/panels"
dir.create(panels, showWarnings = FALSE)
env <- c(FIG_DIR = panels, FIG_DPI = "300", FIG_ONLY = "true")
for (script in c("figures_leibold2019.R", "figures_urban2009.R", "figures_hajian2026.R",
                 "animation_quantgen.R")) {
  status <- system2("Rscript", file.path("data-raw", script), env = paste0(names(env), "=", env),
                    stdout = FALSE)
  if (status != 0) stop(script, " failed")
}

# a letter in the top-left corner of each panel
panel <- function(file, letter, width = 2700) {
  img <- image_scale(image_read(file.path(panels, file)), as.character(width))
  img <- image_border(img, "white", "0x80", operator = "copy")
  image_annotate(img, letter, size = 90, weight = 700, gravity = "northwest",
                 location = "+12+0")
}

# Figure 3: replications --------------------------------------------------------------------
fig3 <- image_append(c(panel("leibold_signatures.png", "(a)"),
                       panel("urban_delay.png", "(b)"),
                       panel("hajian_comparison.png", "(c)")), stack = TRUE)
image_write(fig3, "paper/figures/fig3_replications.png", density = 300)

# Figure 4: local adaptation along a gradient and evolutionary rescue -----------------------
fig4 <- image_append(c(panel("quantgen_gradient.png", "(a)", width = 3000),
                       panel("quantgen_rescue.png", "(b)", width = 2250)), stack = TRUE)
image_write(fig4, "paper/figures/fig4_ranges.png", density = 300)
