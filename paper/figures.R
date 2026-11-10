# Figures of the comdynamite preprint. Run from the package root:
#   Rscript paper/figures.R
# Uses only results saved by the scripts in data-raw/ and vignettes/, so it runs in seconds.
library(ggplot2)
library(patchwork)

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7", "#999999")
dir.create("paper/figures", showWarnings = FALSE)

# Figure 1: structure of the package -------------------------------------------------------
boxes <- data.frame(
  id = c("pool", "sim", "obs", "data", "null", "ana", "eval"),
  x = c(1.3, 1.3, 4.1, 4.1, 6.9, 6.9, 4.1),
  y = c(4.2, 2.4, 4.2, 2.4, 4.2, 2.4, 0.6),
  title = c("Regional pool", "Simulators", "Observed data", "mc_data()",
            "Null models", "Analyses", "Simulation studies"),
  body = c("trees (treeshapesim, gen3sis)\ntrait_bm(), trait_ou(), trait_mk()\ntrait_rates(), history_rates()",
           "simulate_metacommunity()\nsimulate_static(), assemble_community()\nsimulate_damocles(), food_web_dynamics()\nsample_individuals()",
           "community, tree, traits,\nenvironment, coordinates",
           "one object with names\nchecked once",
           "torus_permutations()\nrigid_permutations()\npatch_permutations()",
           "node_test(), pcc_test(), pcps_test()\nphylo_beta_test(), trait_partition()\ntstats(), trait_gradient(), damocles_fit()",
           "run_scenarios(), summary_stats()\nequilibrium_check()\ntype I error, power, recovery"),
  fill = c("pool", "sim", "data", "data", "ana", "ana", "eval")
)
arrows <- data.frame(
  x = c(1.3, 2.6, 4.1, 5.4, 6.9, 6.9, 2.8),
  y = c(3.65, 2.4, 3.65, 2.4, 3.65, 1.85, 0.6),
  xend = c(1.3, 2.8, 4.1, 5.6, 6.9, 5.4, 1.3),
  yend = c(2.95, 2.4, 2.95, 2.4, 2.95, 0.85, 1.75)
)
fig1 <- ggplot() +
  geom_segment(data = arrows[1:6, ], aes(x, y, xend = xend, yend = yend),
               arrow = arrow(length = unit(2.2, "mm"), type = "closed"), colour = "grey35") +
  geom_curve(data = arrows[7, ], aes(x, y, xend = xend, yend = yend), curvature = -0.3,
             arrow = arrow(length = unit(2.2, "mm"), type = "closed"), colour = okabe[4],
             linetype = "dashed") +
  annotate("text", x = 1.55, y = 0.75, label = "known truth", colour = okabe[4],
           size = 3, fontface = "italic", hjust = 0) +
  geom_tile(data = boxes, aes(x, y, fill = fill), width = 2.6, height = 1.1,
            colour = "grey30", alpha = 0.35) +
  geom_text(data = boxes, aes(x, y + 0.36, label = title), fontface = "bold", size = 3.3) +
  geom_text(data = boxes, aes(x, y - 0.12, label = body), size = 2.3, lineheight = 0.95,
            family = "mono") +
  scale_fill_manual(values = c(pool = okabe[3], sim = okabe[2], data = okabe[7],
                               ana = okabe[1], eval = okabe[6]), guide = "none") +
  coord_cartesian(xlim = c(-0.05, 8.25), ylim = c(-0.05, 4.85), expand = FALSE) +
  theme_void()
ggsave("paper/figures/fig1_structure.png", fig1, width = 7.5, height = 5, dpi = 300, bg = "white")

# Figure 3: two questions answered by simulation ----------------------------------------------
# (a) phylogenetic structure along dispersal, from vignettes/phylogenetic-structure.rds
ps <- readRDS("vignettes/phylogenetic-structure.rds")
ps <- ps[is.finite(ps$ses_mpd), ]
means <- aggregate(ses_mpd ~ dispersal + lambda + competition, ps, mean)
ses <- aggregate(ses_mpd ~ dispersal + lambda + competition, ps,
                 function(z) sd(z) / sqrt(length(z)))
means$se <- ses$ses_mpd
means$scenario <- paste0("lambda = ", means$lambda, ", ",
                         ifelse(means$competition == "none", "no competition",
                                "competition among relatives"))
panel_a <- ggplot(means, aes(dispersal, ses_mpd, colour = scenario, shape = scenario)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_errorbar(aes(ymin = ses_mpd - 2 * se, ymax = ses_mpd + 2 * se), width = 0.06) +
  geom_line() + geom_point(size = 2) +
  scale_x_log10() +
  scale_colour_manual(values = okabe[c(1, 2, 3, 4)], name = NULL) +
  scale_shape_manual(values = c(15, 16, 17, 18), name = NULL) +
  labs(x = "Dispersal (emigration probability)", y = "SES MPD (mean of sites)") +
  theme_bw(base_size = 10) +
  guides(colour = guide_legend(ncol = 1)) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# (b) type I error of free and restricted permutations
edges <- readRDS("data-raw/validation_edges.rds")
edges <- edges[edges$effect == 0, ]
long <- rbind(
  data.frame(test = "node_test(), events", case = edges$kind, null = "free", p = edges$free),
  data.frame(test = "node_test(), events", case = edges$kind, null = "torus", p = edges$torus),
  data.frame(test = "node_test(), events", case = edges$kind, null = "patch", p = edges$patch)
)
long <- long[is.finite(long$p), ]
long$case <- c(blocks = "two blocks", others = "area vs. others", tails = "two ends")[long$case]
em <- readRDS("vignettes/evaluating-methods.rds")
em <- em[em$predictor == "unrelated", ]
long <- rbind(long, data.frame(test = "community ~ predictor", case = "torus landscape",
                               null = em$null_model, p = em$p))
rates <- aggregate(p ~ test + case + null, long, function(p) mean(p <= 0.05))
n <- aggregate(p ~ test + case + null, long, length)
rates$se <- sqrt(rates$p * (1 - rates$p) / n$p)
rates$null <- factor(rates$null, c("free", "torus", "patch"),
                     c("free", "torus", "patch_permutations()"))
rates$label <- paste0(rates$case, "\n(", rates$test, ")")
panel_b <- ggplot(rates, aes(p, label, colour = null, shape = null)) +
  geom_vline(xintercept = 0.05, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = pmax(p - 2 * se, 0), xmax = p + 2 * se), height = 0.2,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  scale_x_continuous(trans = "sqrt", breaks = c(0, 0.05, 0.1, 0.25, 0.5, 0.75)) +
  scale_colour_manual(values = okabe[c(4, 5, 3)], name = "Permutations") +
  scale_shape_manual(values = c(16, 17, 15), name = "Permutations") +
  labs(x = "Rejection rate of a true null hypothesis (5% level)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.direction = "vertical", panel.grid.minor = element_blank())
fig2 <- panel_a + panel_b + plot_annotation(tag_levels = "a") + plot_layout(widths = c(1, 1.1))
ggsave("paper/figures/fig3_simulation_studies.png", fig2, width = 9, height = 4.8, dpi = 300,
       bg = "white")

# Figure 2: reproduction of Thompson et al. (2020) ------------------------------------------
th <- readRDS("data-raw/validation_thompson.rds")
th$niche <- as.numeric(th$niche)
keys <- c("dispersal", "niche", "competition")
med <- aggregate(th[c("beta_spatial", "beta_temporal")], th[keys], median, na.rm = TRUE)
persisted <- aggregate(list(n = is.finite(th$alpha)), th[keys], sum)
med <- merge(med, persisted)
med[med$n < 10, c("beta_spatial", "beta_temporal")] <- NA
long <- rbind(data.frame(med[keys], measure = "spatial beta", value = med$beta_spatial),
              data.frame(med[keys], measure = "temporal beta", value = med$beta_temporal))
long$niche <- factor(signif(long$niche, 2))
fig3 <- ggplot(long, aes(dispersal, niche, fill = value)) +
  geom_tile() +
  facet_grid(competition ~ measure) +
  scale_x_log10() +
  scale_fill_viridis_c(option = "inferno", trans = "log10", na.value = "white",
                       name = "beta\n(gamma/alpha)") +
  labs(x = "Dispersal (emigration probability)", y = "Abiotic niche breadth (sigma)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank())
ggsave("paper/figures/fig2_thompson.png", fig3, width = 7, height = 6.5, dpi = 300, bg = "white")

# Numbers quoted in the text ----------------------------------------------------------------
print(rates[order(rates$test, rates$case, rates$null), c("test", "case", "null", "p", "se")])
print(means[order(means$scenario, means$dispersal), c("scenario", "dispersal", "ses_mpd", "se")])
