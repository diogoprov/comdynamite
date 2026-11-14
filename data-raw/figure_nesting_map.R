# Figure of the nesting map (article "Nesting map"): for each mechanism of six published
# models, whether comdynamite has it, approximates it, gained it as an option for a
# replication, or lacks it. Sources and locations of each item are in dev/nesting-map.md.

library(ggplot2)

models <- c("Thompson et al. 2020", "Urban & De Meester 2009", "Leibold et al. 2019",
            "Cabral & Kreft 2012", "Cabral et al. 2019 (BioGEEM)", "Hagen et al. 2021 (gen3sis)")
# one row per mechanism; one letter per model, in the order above:
# s = in comdynamite, o = option added for a replication, a = approximated,
# m = missing from comdynamite, . = not in that model
map <- c(
  "Gaussian niche on an environmental axis"          = "sssaaa",
  "Pairwise competition coefficients"                = "s..mma",
  "Lottery recruitment into microsites"              = ".oo...",
  "Overlapping generations"                          = ".oo...",
  "Stage structure (seeds, juveniles, adults)"       = "...mm.",
  "Allee effects"                                    = "...mm.",
  "Metabolic scaling of rates"                       = "...mm.",
  "Distance-dependent dispersal"                     = "s.saam",
  "Environment varying in space and time"            = "s...ss",
  "Habitable sites changing through time"            = "....mm",
  "Loss of populations or of whole sites"            = "s.o...",
  "Trait evolution within populations"               = ".aa..a",
  "Variance built by mutation, eroded by selection"  = ".oo...",
  "Microsite heterogeneity within sites"             = ".mm...",
  "Sexual reproduction and recombination"            = "..m...",
  "Speciation"                                       = "....ss",
  "Deep-time steps (years to 170 kyr, Myr runs)"     = "....mm"
)
states <- c(s = "in comdynamite", o = "option added for the replication", a = "approximated",
            m = "missing")
d <- do.call(rbind, lapply(names(map), function(k) {
  data.frame(mechanism = k, model = models, code = strsplit(map[[k]], "")[[1]])
}))
d <- d[d$code != ".", ]
d$state <- factor(states[d$code], states)
d$mechanism <- factor(d$mechanism, rev(names(map)))
d$model <- factor(d$model, models)

p <- ggplot(d, aes(model, mechanism, fill = state)) +
  geom_tile(colour = "white", linewidth = 1) +
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#E69F00", "#D55E00"), drop = FALSE) +
  scale_x_discrete(position = "top", labels = function(x) gsub(" (\\d{4})", "\n\\1", x)) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "bottom")
ggsave("vignettes/articles/nesting_map.png", p, width = 9.5, height = 6.5, dpi = 120,
       bg = "white")
