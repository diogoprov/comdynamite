# Figures and numbers of the replication of Hajian-Forooshani & Chase (2026), from
# data-raw/replication_hajian2026.rds (written by data-raw/replication_hajian2026.R).
# Figures go to vignettes/articles/hajian_*.png; numbers are printed for the article.

library(ggplot2)
res <- readRDS("data-raw/replication_hajian2026.rds")
s <- res$summary
s$extent <- s$n_disturbed # 100 patches, so the number of patches is the percentage
ctrl <- res$control
theme_set(theme_minimal(base_size = 11))
save <- function(p, name, w = 9, h = 5.5) {
  ggsave(file.path("vignettes/articles", name), p, width = w, height = h, dpi = 120,
         bg = "white")
}
s$regime <- factor(s$dispersal, levels = c(1e-5, 1e-4, 1e-3),
                   labels = c("Dispersal limited (d = 1e-5)", "Intermediate (d = 1e-4)",
                              "Dispersal surplus (d = 1e-3)"))

ctrl$regime <- factor(ctrl$dispersal, levels = c(1e-5, 1e-4, 1e-3), labels = levels(s$regime))

# Fig. 2: changes from step 499 to step 700 against the extent of the disturbance, n = 0.25;
# dashed grey lines: the same changes without disturbance
f2 <- s[s$niche == 0.25 & s$step == 700, ]
c2 <- ctrl[ctrl$niche == 0.25 & ctrl$step == 700, ]
scales <- rbind(
  data.frame(f2[c("regime", "extent")], scale = "alpha", value = f2$d_alpha, se = f2$d_alpha_se),
  data.frame(f2[c("regime", "extent")], scale = "gamma", value = f2$d_gamma, se = f2$d_gamma_se))
p2a <- ggplot(scales, aes(extent, value, colour = scale)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_hline(data = rbind(data.frame(regime = c2$regime, scale = "alpha", value = c2$d_alpha),
                          data.frame(regime = c2$regime, scale = "gamma", value = c2$d_gamma)),
             aes(yintercept = value, colour = scale), linetype = 2, alpha = 0.6) +
  geom_errorbar(aes(ymin = value - 1.96 * se, ymax = value + 1.96 * se), width = 0) +
  geom_line() + geom_point(size = 1.5) + facet_wrap(~regime) +
  scale_colour_manual(values = c(alpha = "#D55E00", gamma = "#0072B2"),
                      labels = c(alpha = "Δα", gamma = "Δγ"), name = NULL) +
  labs(x = NULL, y = "Log change in richness")
p2b <- ggplot(f2, aes(extent, d_beta)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_hline(data = c2, aes(yintercept = d_beta), linetype = 2, colour = "grey50") +
  geom_errorbar(aes(ymin = d_beta - 1.96 * d_beta_se, ymax = d_beta + 1.96 * d_beta_se),
                width = 0, colour = "#6A3D9A") +
  geom_line(colour = "#6A3D9A") + geom_point(size = 1.5, colour = "#6A3D9A") +
  facet_wrap(~regime) + labs(x = "% of landscape disturbed", y = "Δβ")
save(patchwork::wrap_plots(p2a, p2b, ncol = 1), "hajian_fig2.png")

# Fig. 3: changes from step 499 to steps 601-1000 against the extent, n = 0.2, 0.3 and 0.4
f3 <- s[s$niche %in% c(0.2, 0.3, 0.4) & s$step %in% c(601, seq(650, 1000, 50)), ]
f3$niche_lab <- factor(paste("n =", f3$niche))
long3 <- rbind(data.frame(f3[c("regime", "niche_lab", "extent", "step")], scale = "Δα",
                          value = f3$d_alpha),
               data.frame(f3[c("regime", "niche_lab", "extent", "step")], scale = "Δγ",
                          value = f3$d_gamma))
p3 <- ggplot(long3, aes(extent, value, colour = step, group = interaction(step, scale),
                        linetype = scale)) +
  geom_hline(yintercept = 0, linetype = 3) + geom_line() +
  facet_grid(niche_lab ~ regime) + scale_colour_viridis_c(option = "magma", end = 0.85,
                                                          name = "Step") +
  labs(x = "% of landscape disturbed", y = "Log change in richness from step 499",
       linetype = NULL)
save(p3, "hajian_fig3.png", h = 7)

# numbers for the article
cat("runs per setting:", res$n_rep, "; minutes:", round(res$minutes, 1), "\n")
print(f2[f2$extent %in% c(1, 10, 25, 50, 75, 100),
         c("regime", "extent", "d_alpha", "d_gamma", "d_beta", "d_beta_se")], digits = 2)
crossing <- by(f2, f2$regime, function(d) {
  d <- d[order(d$extent), ]
  c(min = min(d$d_beta), extent_min = d$extent[which.min(d$d_beta)], max = max(d$d_beta),
    extent_max = d$extent[which.max(d$d_beta)])
})
print(crossing)
cat("control without disturbance, steps 700 and 1000:\n")
print(ctrl[ctrl$step %in% c(700, 1000), c("niche", "dispersal", "step", "d_alpha", "d_gamma",
                                          "d_beta", "d_gamma_se")], digits = 2, row.names = FALSE)
ext <- s[s$n_extinct > 0, ]
cat("settings with extinct replicates:\n")
print(stats::aggregate(n_extinct ~ niche + dispersal + n_disturbed, ext, max))
