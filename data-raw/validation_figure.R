# Draws vignettes/articles/thompson_beta.png from data-raw/validation_thompson.rds,
# with the layout of Figure 4 of Thompson et al. (2020). Cells where species
# persisted in fewer than 10 of the 30 landscapes are left blank, as in that figure.
# Uses ggplot2. Run from the package root.
library(ggplot2)

results <- readRDS("data-raw/validation_thompson.rds")
results$niche <- as.numeric(results$niche)
keys <- c("dispersal", "niche", "competition")
medians <- aggregate(results[c("beta_spatial", "beta_temporal")], results[keys], median, na.rm = TRUE)
persisted <- aggregate(list(n = is.finite(results$alpha)), results[keys], sum)
medians <- merge(medians, persisted)
medians[medians$n < 10, c("beta_spatial", "beta_temporal")] <- NA
long <- rbind(data.frame(medians[keys], measure = "spatial beta", value = medians$beta_spatial),
              data.frame(medians[keys], measure = "temporal beta", value = medians$beta_temporal))
long$niche <- factor(long$niche)

figure <- ggplot(long, aes(dispersal, niche, fill = value)) +
  geom_tile() +
  facet_grid(competition ~ measure) +
  scale_x_log10() +
  scale_fill_viridis_c(option = "inferno", trans = "log10", na.value = "white",
                       name = "species\nrichness") +
  labs(x = "dispersal", y = "abiotic niche breadth") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
ggsave("vignettes/articles/thompson_beta.png", figure, width = 7, height = 6, dpi = 120)
