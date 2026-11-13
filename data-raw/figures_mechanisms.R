# Figures and numbers of the tests of three missing mechanisms (lottery recruitment,
# overlapping generations, variance eroded by selection) in the replications of Urban & De
# Meester (2009) and Leibold et al. (2019), from data-raw/mechanisms_urban2009.rds and
# data-raw/mechanisms_leibold2019_recruitment.rds. Figures go to vignettes/articles/.

library(ggplot2)
theme_set(theme_minimal(base_size = 11))
save <- function(p, name, w = 9, h = 4.5) {
  ggsave(file.path("vignettes/articles", name), p, width = w, height = h, dpi = 120,
         bg = "white")
}
mech_levels <- c("Beverton-Holt", "lottery, discrete", "lottery, overlapping")
var_label <- function(v) factor(ifelse(v, "eroded by selection", "mutation and drift"),
                                c("mutation and drift", "eroded by selection"))

# Urban & De Meester (2009) ---------------------------------------------------------------
u <- readRDS("data-raw/mechanisms_urban2009.rds")
f2 <- u$fig2
f2$label <- factor(f2$label, mech_levels)
f2$share1 <- f2$n1 / pmax(f2$n1 + f2$n2, 1)
f2$first_wins <- ifelse(f2$delay > 0, f2$n1 > f2$n2, f2$n2 > f2$n1)
f2$mono <- pmax(f2$n1, f2$n2) / pmax(f2$n1 + f2$n2, 1) > 0.9
cat("\nUrban, share of the new site held by species 1, delays -10, 0, 10:\n")
tab <- aggregate(share1 ~ genetic + label + variance_selection + delay,
                 f2[f2$delay %in% c(-10, 0, 10), ], function(v) round(mean(v), 2))
print(reshape(tab, idvar = c("genetic", "label", "variance_selection"), timevar = "delay",
              direction = "wide"), row.names = FALSE)
cat("\nUrban, runs won by the first arrival (leads 1, 2, 5, 10) and monopolised runs:\n")
w <- aggregate(first_wins ~ genetic + label + variance_selection + abs(delay),
               f2[f2$delay != 0, ], function(v) round(mean(v), 2))
names(w)[4] <- "lead"
mono <- aggregate(mono ~ genetic + label + variance_selection, f2,
                  function(v) round(mean(v), 2))
print(merge(reshape(w, idvar = c("genetic", "label", "variance_selection"), timevar = "lead",
                    direction = "wide"), mono), row.names = FALSE)

both <- f2[f2$genetic %in% c("mutation_both", "standing_both"), ]
m2 <- aggregate(share1 ~ genetic + label + variance_selection + delay, both, mean)
m2$variance <- var_label(m2$variance_selection)
m2$genetic <- factor(ifelse(m2$genetic == "mutation_both", "both evolve, mutation",
                            "both evolve, standing variation"))
none <- aggregate(share1 ~ label + delay, f2[f2$genetic == "none", ], mean)
p <- ggplot(m2, aes(delay, share1, colour = label, linetype = variance)) +
  geom_hline(yintercept = 0.5, colour = "grey80") +
  geom_line(data = none, aes(delay, share1, group = label), colour = "grey60",
            linetype = "dotted", inherit.aes = FALSE) +
  geom_line(linewidth = 0.7) + facet_wrap(~genetic) +
  scale_colour_manual(values = c("#0072B2", "#E69F00", "#009E73")) +
  labs(x = "Delay (generations; positive: species 1 enters first)",
       y = "Share of the new site held by species 1", colour = NULL, linetype = "Variance")
save(p, "mechanisms_urban_delay.png", w = 10, h = 3.8)

f4 <- u$fig4
f4$label <- factor(f4$label, mech_levels)
m4 <- aggregate(cbind(dom = (n1 - n2) / 1000) ~ genetic + label + variance_selection + m1 + m2,
                f4, mean)
m4$setting <- ifelse(m4$genetic == "none", "no evolution",
                     paste("mutation,", as.character(var_label(m4$variance_selection))))
m4$setting <- factor(m4$setting, c("no evolution", "mutation, mutation and drift",
                                   "mutation, eroded by selection"))
p <- ggplot(m4, aes(factor(m1), factor(m2), fill = dom)) + geom_tile() +
  facet_grid(label ~ setting) +
  scale_fill_gradient2(low = "#D55E00", mid = "white", high = "#0072B2", limits = c(-1, 1),
                       name = "(N1 - N2) / K") +
  labs(x = "Dispersal of species 1", y = "Dispersal of species 2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save(p, "mechanisms_urban_dispersal.png", w = 10, h = 8)
# summary of the grid: mean dominance of the better disperser off the diagonal
off <- m4[m4$m1 != m4$m2, ]
off$better <- ifelse(off$m1 > off$m2, off$dom, -off$dom)
cat("\nUrban, mean dominance of the better disperser off the diagonal:\n")
print(aggregate(better ~ setting + label, off, function(v) round(mean(v), 2)))

f3 <- u$fig3[u$fig3$step == 2000, ]
f3$label <- factor(f3$label, mech_levels)
f3$winner <- ifelse(f3$new1 > f3$new2, 1, 2)
f3$own <- ifelse(f3$winner == 1, f3$src1_sp1 / pmax(f3$src1_sp1 + f3$src1_sp2, 1),
                 f3$src2_sp2 / pmax(f3$src2_sp1 + f3$src2_sp2, 1))
f3$loser <- ifelse(f3$winner == 1, f3$src2_sp2 / pmax(f3$src2_sp1 + f3$src2_sp2, 1),
                   f3$src1_sp1 / pmax(f3$src1_sp1 + f3$src1_sp2, 1))
cat("\nUrban, share of winner and loser of the new site in their own source, generation 2000:\n")
print(aggregate(cbind(own, loser) ~ m + label + variance_selection, f3,
                function(v) round(mean(v), 2)))

# Leibold et al. (2019) -------------------------------------------------------------------
l <- readRDS("data-raw/mechanisms_leibold2019_recruitment.rds")
lab <- function(d) {
  paste(c("Beverton-Holt", "lottery, discrete", "lottery, overlapping")[
    match(paste(d$recruitment, d$mortality), c("beverton_holt 1", "lottery 1", "lottery 0.1"))])
}
gen_levels <- c("none", "constant", "mutation_from_G", "mutation_from_zero",
                "mutation_selection")
summ <- do.call(rbind, lapply(l, function(r) {
  cbind(label = lab(r$design), genetic = r$design$genetic, rep = r$design$rep, r$summary)
}))
summ$label <- factor(summ$label, mech_levels)
summ$genetic <- factor(summ$genetic, gen_levels)
cat("\nLeibold, runs with both species at generations 5000 and 10000:\n")
pers <- aggregate(both ~ label + genetic + step, summ[summ$step %in% c(5000, 10000), ], mean)
print(reshape(pers, idvar = c("label", "genetic"), timevar = "step", direction = "wide"),
      row.names = FALSE)
cat("\nLeibold, regional neutrality (generations 4000-5000, runs with both species):\n")
late <- summ[summ$step > 4000 & summ$step <= 5000 & summ$both, ]
late$fitdiff <- (abs(late$fitdiff_e40) + abs(late$fitdiff_e60)) / 2
print(aggregate(cbind(env_comp, env_trait, fitdiff) ~ label + genetic, late,
                function(v) round(mean(v, na.rm = TRUE), 3), na.action = na.pass))
cls <- do.call(rbind, lapply(l, function(r) {
  alive <- r$summary$both[r$summary$step == 5000]
  if (!isTRUE(alive)) return(NULL)
  data.frame(label = lab(r$design), genetic = r$design$genetic, rep = r$design$rep,
             class = names(r$classes), share = as.numeric(r$classes))
}))
cls$label <- factor(cls$label, mech_levels)
cls$genetic <- factor(cls$genetic, gen_levels)
cls_m <- aggregate(share ~ label + genetic + class, cls, mean)
cls_n <- aggregate(rep ~ label + genetic, unique(cls[, c("label", "genetic", "rep")]), length)
cat("\nLeibold, classes of patches at generation 5000 (runs with both species):\n")
print(merge(reshape(cls_m, idvar = c("label", "genetic"), timevar = "class",
                    direction = "wide"), cls_n), row.names = FALSE, digits = 2)
cls_m$class <- factor(cls_m$class, c("preadapted", "switched", "shared"))
p <- ggplot(cls_m, aes(genetic, share, fill = class)) + geom_col() + facet_wrap(~label) +
  geom_text(data = cls_n, aes(genetic, 1.05, label = paste0("n=", rep)), inherit.aes = FALSE,
            size = 3) +
  scale_fill_manual(values = c(preadapted = "#0072B2", switched = "#E69F00",
                               shared = "#CCCCCC")) +
  scale_x_discrete(labels = c("none", "constant G", "mutation\nfrom G", "mutation\nfrom 0",
                              "eroded by\nselection")) +
  theme(axis.text.x = element_text(size = 7.5)) +
  labs(x = NULL, y = "Share of occupied patches", fill = NULL,
       caption = "Nm = 10; runs with both species at generation 5000 (n of 10)")
save(p, "mechanisms_leibold_classes.png", w = 11, h = 3.8)
cat("\nLeibold, mean genetic variance at the end:\n")
print(do.call(rbind, lapply(l, function(r) {
  data.frame(label = lab(r$design), genetic = r$design$genetic,
             G = if (is.null(r$variance)) NA else mean(r$variance, na.rm = TRUE))
})) |> aggregate(G ~ label + genetic, data = _, function(v) round(mean(v), 1)))
