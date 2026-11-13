# Figures and numbers of the replication of Urban & De Meester (2009), from
# data-raw/replication_urban2009.rds (written by data-raw/replication_urban2009.R).
# Figures go to vignettes/articles/urban_*.png; numbers are printed for the article.

library(ggplot2)
res <- readRDS("data-raw/replication_urban2009.rds")
theme_set(theme_minimal(base_size = 11))
save <- function(p, name, w = 9, h = 4.5) {
  ggsave(file.path("vignettes/articles", name), p, width = w, height = h, dpi = 120,
         bg = "white")
}
labels <- c(none = "No evolution",
            mutation_sp1 = "Species 1 evolves, mutation",
            mutation_both = "Both evolve, mutation",
            standing_sp1 = "Species 1 evolves, standing variation",
            standing_both = "Both evolve, standing variation",
            constant_sp1 = "Species 1 evolves, constant G",
            constant_both = "Both evolve, constant G")

# Fig. 2: final abundance in the new site against the delay
f2 <- res$fig2
m2 <- aggregate(cbind(n1, n2) ~ genetic + delay, f2, mean)
l2 <- rbind(data.frame(m2[1:2], species = "species 1", n = m2$n1),
            data.frame(m2[1:2], species = "species 2", n = m2$n2))
l2$genetic <- factor(labels[l2$genetic], labels)
p2 <- ggplot(l2, aes(delay, n, colour = species)) +
  geom_vline(xintercept = 0, colour = "grey70") + geom_line(linewidth = 0.7) +
  facet_wrap(~genetic, nrow = 2) +
  scale_colour_manual(values = c("species 1" = "#0072B2", "species 2" = "#D55E00")) +
  labs(x = "Delay (generations; positive: species 1 enters first)",
       y = "Mean final abundance in the new site", colour = NULL)
save(p2, "urban_delay.png", w = 11, h = 5)

# share of runs won by the first species to enter (abundance > half of the total), by delay
f2$first_wins <- ifelse(f2$delay > 0, f2$n1 > f2$n2, f2$n2 > f2$n1)
win <- aggregate(first_wins ~ genetic + abs(delay), f2[f2$delay != 0, ], mean)
names(win)[2] <- "lead"
cat("\nShare of runs won by the first species, by lead (generations):\n")
print(reshape(win[win$lead %in% c(1, 2, 5, 10, 20), ], idvar = "genetic", timevar = "lead",
              direction = "wide"), row.names = FALSE, digits = 2)
cat("\nSpecies 1 evolving: share of its abundance in the new site when it enters 10 later, at",
    "the same time and 10 earlier:\n")
sh <- aggregate(cbind(n1, n2) ~ genetic + delay, f2[f2$delay %in% c(-10, 0, 10), ], mean)
sh$share1 <- round(sh$n1 / (sh$n1 + sh$n2), 2)
print(sh[order(sh$genetic, sh$delay), c("genetic", "delay", "share1")], row.names = FALSE)
# U-shaped distribution: share of runs where one species holds more than 90% of the site
f2$tot <- f2$n1 + f2$n2
f2$mono <- pmax(f2$n1, f2$n2) / pmax(f2$tot, 1) > 0.9
cat("\nRuns in which one species holds more than 90% of the new site (all delays):\n")
print(aggregate(mono ~ genetic, f2, function(v) round(mean(v), 2)))

# Fig. 4: relative proportional dominance over the dispersal grid
f4 <- res$fig4
m4 <- aggregate(cbind(dom = (n1 - n2) / 1000) ~ m1 + m2, f4, mean)
p4 <- ggplot(m4, aes(factor(m1), factor(m2), fill = dom)) +
  geom_tile() +
  scale_fill_gradient2(low = "#D55E00", mid = "white", high = "#0072B2", limits = c(-1, 1),
                       name = "(N1 - N2) / K") +
  labs(x = "Dispersal of species 1", y = "Dispersal of species 2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save(p4, "urban_dispersal.png", w = 6.5, h = 5.2)
cat("\nDominance (N1 - N2) / 1000, rows: dispersal of species 1, columns: species 2:\n")
print(round(xtabs(dom ~ m1 + m2, m4), 2))
cat("\nSource sites at the end (mean abundance of the resident species):\n")
print(round(xtabs(cbind(src1) ~ m1 + m2, aggregate(src1 ~ m1 + m2, f4, mean))))

# Fig. 3: does the winner of the new site lose its own source?
f3 <- res$fig3
last <- f3[f3$step == max(f3$step), ]
last$winner <- ifelse(last$new1 > last$new2, 1, 2)
last$own_share <- ifelse(last$winner == 1, last$src1_sp1 / pmax(last$src1_sp1 + last$src1_sp2, 1),
                         last$src2_sp2 / pmax(last$src2_sp1 + last$src2_sp2, 1))
last$loser_share <- ifelse(last$winner == 1, last$src2_sp2 / pmax(last$src2_sp1 + last$src2_sp2, 1),
                           last$src1_sp1 / pmax(last$src1_sp1 + last$src1_sp2, 1))
cat("\nAt 2000 generations: share of each species in its own source site, winner and loser of",
    "the new site:\n")
print(aggregate(cbind(own_share, loser_share) ~ m, last, function(v) round(mean(v), 2)))
f3$winner <- last$winner[match(paste(f3$m, f3$rep), paste(last$m, last$rep))]
f3$own <- ifelse(f3$winner == 1, f3$src1_sp1 / pmax(f3$src1_sp1 + f3$src1_sp2, 1),
                 f3$src2_sp2 / pmax(f3$src2_sp1 + f3$src2_sp2, 1))
f3$z_own <- ifelse(f3$winner == 1, f3$z_src1_sp1, f3$z_src2_sp2)
m3 <- aggregate(cbind(own, z_own) ~ m + step, f3, mean)
l3 <- rbind(data.frame(m3[1:2], what = "Share of the winner in its own source", value = m3$own),
            data.frame(m3[1:2], what = "Genetic mean of the winner in its own source",
                       value = m3$z_own))
p3 <- ggplot(l3, aes(step, value, colour = factor(m))) +
  geom_line(linewidth = 0.7) + facet_wrap(~what, scales = "free_y") +
  scale_colour_manual(values = c("0.01" = "#56B4E9", "0.05" = "#009E73", "0.1" = "#D55E00")) +
  labs(x = "Generation", y = NULL, colour = "Dispersal")
save(p3, "urban_source.png", w = 9, h = 3.4)
