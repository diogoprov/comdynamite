# Figures and numbers of the replication of Leibold et al. (2019), from
# data-raw/replication_leibold2019.rds (written by data-raw/replication_leibold2019.R).
# Figures go to vignettes/articles/leibold_*.png; numbers are printed for the article.

library(ggplot2)
in_file <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(in_file)) in_file <- "data-raw/replication_leibold2019.rds"
res <- readRDS(in_file)

parts <- lapply(res, function(r) {
  cbind(G = r$G, nm = r$nm, landscape_type = r$landscape_type, rep = r$rep, r$summary)
})
# runs with one patch type have their own trait and fitness columns
cols <- unique(unlist(lapply(parts, names)))
long <- do.call(rbind, lapply(parts, function(d) {
  d[setdiff(cols, names(d))] <- NA
  d[cols]
}))
long$fitdiff <- (abs(long$fitdiff_e40) + abs(long$fitdiff_e60)) / 2
two <- long[long$landscape_type == "two", ]
g_lab <- function(g) factor(paste("G =", g), paste("G =", sort(unique(g))))
two$G_f <- g_lab(two$G)
two$nm_f <- factor(paste("Nm =", two$nm), paste("Nm =", sort(unique(two$nm))))
pal <- c("#999999", "#56B4E9", "#009E73", "#E69F00", "#D55E00", "#CC79A7")
theme_set(theme_minimal(base_size = 11))
save <- function(p, name, w = 8, h = 4) {
  ggsave(file.path("vignettes/articles", name), p, width = w, height = h, dpi = 120,
         bg = "white")
}

# 1. persistence of both species (their Fig. S1)
pers <- aggregate(both ~ step + G_f + nm_f, two, mean)
p1 <- ggplot(pers, aes(step, both, colour = G_f)) +
  geom_line(linewidth = 0.8) + facet_wrap(~nm_f) + scale_colour_manual(values = pal) +
  labs(x = "Generation", y = "Runs with both species", colour = NULL)
save(p1, "leibold_persistence.png", h = 3.2)

# 2. signatures of regional neutrality while both species persist (Nm = 10)
sig <- two[two$nm == 10 & two$both, ]
sig_long <- rbind(
  data.frame(sig[, c("step", "G_f")], what = "Patch type, composition", value = sig$env_comp),
  data.frame(sig[, c("step", "G_f")], what = "Patch type, mean optimum", value = sig$env_trait),
  data.frame(sig[, c("step", "G_f")], what = "Fitness difference", value = sig$fitdiff))
sig_long <- sig_long[is.finite(sig_long$value), ]
sig_mean <- aggregate(value ~ step + G_f + what, sig_long, mean)
p2 <- ggplot(sig_mean, aes(step, value, colour = G_f)) +
  geom_line(linewidth = 0.7) + facet_wrap(~what, scales = "free_y") +
  scale_colour_manual(values = pal) +
  labs(x = "Generation", y = "Mean over runs with both species", colour = NULL)
save(p2, "leibold_signatures.png", w = 9, h = 3.2)

# 3. mean genetic optimum of each species in each patch type (their Fig. 3B), Nm = 10
zl <- do.call(rbind, lapply(c("z_sp1_e40", "z_sp1_e60", "z_sp2_e40", "z_sp2_e60"), function(v) {
  d <- two[two$nm == 10 & two$G %in% c(0, 50, 200), c("step", "G_f", "rep")]
  data.frame(d, species = sub("z_(sp[12])_.*", "\\1", v),
             patch = paste("patch type", sub(".*_e", "", v)), z = two[[v]][two$nm == 10 &
                                                                         two$G %in% c(0, 50, 200)])
}))
zl <- zl[is.finite(zl$z), ]
zm <- aggregate(z ~ step + G_f + species + patch, zl, mean)
p3 <- ggplot(zm, aes(step, z, colour = species, linetype = patch)) +
  geom_line(linewidth = 0.7) + facet_wrap(~G_f) +
  scale_colour_manual(values = c(sp1 = "#0072B2", sp2 = "#D55E00")) +
  labs(x = "Generation", y = "Mean genetic optimum", colour = NULL, linetype = NULL)
save(p3, "leibold_optima.png", h = 3.2)

# 4. classes of patches at generation 5000 (their Fig. 5)
cls <- do.call(rbind, lapply(res, function(r) {
  if (r$landscape_type != "two" || is.null(r$classes)) return(NULL)
  data.frame(G = r$G, nm = r$nm, rep = r$rep, class = names(r$classes),
             share = as.numeric(r$classes))
}))
# only runs with both species at generation 5000: after the loss of one, every patch is
# trivially dominated by the survivor
alive5000 <- two[two$step == 5000 & two$both, c("G", "nm", "rep")]
cls <- merge(cls, alive5000)
cls_n <- aggregate(rep ~ G + nm, unique(cls[, c("G", "nm", "rep")]), length)
cls_m <- aggregate(share ~ G + nm + class, cls, mean)
cls_m$class <- factor(cls_m$class, c("preadapted", "switched", "shared"))
cls_n$nm_f <- factor(paste("Nm =", cls_n$nm), paste("Nm =", sort(unique(cls_m$nm))))
cls_m$nm_f <- factor(paste("Nm =", cls_m$nm), paste("Nm =", sort(unique(cls_m$nm))))
p4 <- ggplot(cls_m, aes(factor(G), share, fill = class)) +
  geom_col() + facet_wrap(~nm_f) +
  scale_fill_manual(values = c(preadapted = "#0072B2", switched = "#E69F00",
                               shared = "#CCCCCC")) +
  geom_text(data = cls_n, aes(factor(G), 1.04, label = paste0("n=", rep)), inherit.aes = FALSE,
            size = 3) +
  labs(x = "G", y = "Share of occupied patches", fill = NULL,
       caption = "Runs with both species at generation 5000 (n of 10)")
save(p4, "leibold_classes.png", h = 3.2)

# 5. regional drift with two patch types or one (their Fig. S2), Nm = 10
dr <- long[long$nm == 10 & long$G %in% c(0, 200), ]
dr$setting <- paste(ifelse(dr$landscape_type == "two", "two patch types", "one patch type"),
                    ifelse(dr$G > 0, "evolution", "no evolution"), sep = ", ")
p5 <- ggplot(dr, aes(step, prop_sp1, group = rep)) +
  geom_line(alpha = 0.6, colour = "#0072B2") + facet_wrap(~setting, nrow = 1) +
  labs(x = "Generation", y = "Proportion of species 1")
save(p5, "leibold_drift.png", w = 10, h = 3)

# numbers for the text
cat("\nRuns with both species at generations 1000, 5000 and 10000:\n")
print(reshape(pers[pers$step %in% c(1000, 5000, 10000), ], idvar = c("G_f", "nm_f"),
              timevar = "step", direction = "wide"), row.names = FALSE, digits = 2)
cat("\nSignatures, Nm = 10, mean over generations 4000-5000 in runs with both species:\n")
late <- sig_long[sig_long$step > 4000 & sig_long$step <= 5000, ]
print(aggregate(value ~ G_f + what, late, function(v) round(mean(v), 3)))
early <- sig_long[sig_long$step <= 200, ]
cat("\nSame, first 200 generations:\n")
print(aggregate(value ~ G_f + what, early, function(v) round(mean(v), 3)))
cat("\nClasses at generation 5000:\n")
print(reshape(cls_m[, c("G", "nm", "class", "share")], idvar = c("G", "nm"),
              timevar = "class", direction = "wide"), row.names = FALSE, digits = 2)
cat("\nOne patch type, runs with both species at 10000:\n")
print(aggregate(both ~ G, long[long$landscape_type == "one" & long$step == 10000, ], mean))
cat("\nMean time to loss of a species (generations, runs that lost one):\n")
lost <- aggregate(both ~ G + nm + landscape_type + rep, long, function(b) {
  if (all(b)) NA else min(which(!b)) * 100
}, na.action = na.pass)
print(aggregate(both ~ G + nm + landscape_type, lost, function(v) round(median(v))))
