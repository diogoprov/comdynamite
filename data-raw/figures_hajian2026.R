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

# Comparison with their archived simulations (data-raw/hajian2026_archived_summary.rds, written
# by data-raw/compare_hajian2026.R), kernel exp(-0.1 d) as in the article and here
arch <- readRDS("data-raw/hajian2026_archived_summary.rds")
theirs <- arch$changes[arch$changes$exp_kernal == 0.1, ]
names(theirs)[names(theirs) == "niche_breadth"] <- "niche"
names(theirs)[names(theirs) == "dispersal_prob"] <- "dispersal"
names(theirs)[names(theirs) == "n.patches.pulsed"] <- "n_disturbed"
names(theirs)[names(theirs) == "post.data.ts"] <- "step"
both <- merge(s, theirs, by = c("niche", "dispersal", "n_disturbed", "step"),
              suffixes = c("_ours", "_theirs"))
agree <- do.call(rbind, lapply(c("d_alpha", "d_gamma", "d_beta"), function(v) {
  o <- both[[paste0(v, "_ours")]]
  t <- both[[paste0(v, "_theirs")]]
  se <- sqrt(both[[paste0(v, "_se_ours")]]^2 + both[[paste0(v, "_se_theirs")]]^2)
  data.frame(measure = v, n = length(o), r = stats::cor(o, t),
             mean_abs_diff = mean(abs(o - t)), mean_diff = mean(o - t),
             within_95 = mean(abs(o - t) <= 1.96 * se))
}))
cat("agreement with their archived simulations (kernel 0.1):\n")
print(agree, digits = 3)
long_cmp <- do.call(rbind, lapply(c("d_alpha", "d_gamma", "d_beta"), function(v) {
  data.frame(measure = unname(c(d_alpha = "Δα", d_gamma = "Δγ", d_beta = "Δβ")[v]),
             regime = both$regime, ours = both[[paste0(v, "_ours")]],
             theirs = both[[paste0(v, "_theirs")]])
}))
long_cmp$measure <- factor(long_cmp$measure, c("Δα", "Δγ", "Δβ"))
pc <- ggplot(long_cmp, aes(theirs, ours, colour = regime)) +
  geom_abline(colour = "grey60") + geom_point(size = 0.8, alpha = 0.5) +
  facet_wrap(~measure, scales = "free") +
  scale_colour_manual(values = c("#009E73", "#E69F00", "#56B4E9"), name = NULL) +
  labs(x = "Hajian-Forooshani & Chase (archived simulations)", y = "comdynamite") +
  theme(legend.position = "bottom")
save(pc, "hajian_comparison.png", h = 4)

# Fig. 4: log change in gamma during the pulse (step 499 to 601) and after it (601 to 1000).
# Ours: the change after the pulse is the difference of the mean log changes at steps 1000 and
# 601, equal to the mean change per replicate except where whole metacommunities went extinct
wide <- function(d) {
  a <- d[d$step == 601, c("niche", "dispersal", "n_disturbed", "d_gamma")]
  b <- d[d$step == 1000, c("niche", "dispersal", "n_disturbed", "d_gamma")]
  m <- merge(a, b, by = c("niche", "dispersal", "n_disturbed"), suffixes = c("_601", "_1000"))
  m$during <- m$d_gamma_601
  m$after <- m$d_gamma_1000 - m$d_gamma_601
  m
}
o4 <- wide(s)
c4 <- wide(ctrl)
t4 <- arch$fig4[arch$fig4$exp_kernal == 0.1 & arch$fig4$niche_breadth %in% unique(s$niche), ]
t4 <- data.frame(niche = t4$niche_breadth, dispersal = t4$dispersal_prob,
                 n_disturbed = t4$n.patches.pulsed, during = t4$during, after = t4$after)
stack4 <- function(d, source) {
  rbind(data.frame(d[c("niche", "dispersal", "n_disturbed")], phase = "During the pulse",
                   value = d$during, source = source),
        data.frame(d[c("niche", "dispersal", "n_disturbed")], phase = "After the pulse",
                   value = d$after, source = source))
}
l4 <- rbind(stack4(o4, "comdynamite"), stack4(t4, "Hajian-Forooshani & Chase"))
l4$phase <- factor(l4$phase, c("During the pulse", "After the pulse"))
l4$regime <- factor(l4$dispersal, levels = c(1e-5, 1e-4, 1e-3), labels = levels(s$regime))
l4$niche <- factor(l4$niche)
lc <- stack4(c4, "control")
lc$phase <- factor(lc$phase, levels(l4$phase))
lc$regime <- factor(lc$dispersal, levels = c(1e-5, 1e-4, 1e-3), labels = levels(s$regime))
lc$niche <- factor(lc$niche)
p4 <- ggplot(l4, aes(n_disturbed, value, colour = niche)) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_hline(data = lc, aes(yintercept = value, colour = niche), linetype = 2, alpha = 0.6) +
  geom_line(data = l4[l4$source == "comdynamite", ]) +
  geom_point(data = l4[l4$source != "comdynamite", ], size = 1, shape = 1) +
  facet_grid(phase ~ regime, scales = "free_y") +
  scale_colour_viridis_d(end = 0.85, name = "Niche breadth n") +
  labs(x = "% of landscape disturbed", y = "Log change in regional richness")
save(p4, "hajian_fig4.png", h = 5.5)
cat("Fig. 4, n = 0.25, extents 10, 50, 90:\n")
sel <- function(d) d[d$niche == 0.25 & d$n_disturbed %in% c(10, 50, 90),
                     c("dispersal", "n_disturbed", "during", "after")]
print(merge(sel(o4), sel(t4), by = c("dispersal", "n_disturbed"),
            suffixes = c("_ours", "_theirs")), digits = 2)
cat("control, after the pulse (601 to 1000):\n")
print(c4[, c("niche", "dispersal", "after")], digits = 2)
