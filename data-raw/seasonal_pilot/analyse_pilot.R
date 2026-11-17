# ABC-RF model choice on the identifiability pilot of issue #97 (run_pilot.R): out-of-bag
# error and confusion matrix of the four community models for 1, 2 and 3 field years, with
# tadpole counts alone or with calling males, and the error of the pair SS x PM, the main
# confusion between models. Output: data-raw/seasonal_pilot/pilot_results.rds and the figure
# data-raw/seasonal_pilot/pilot_error.png.

suppressMessages(devtools::load_all(quiet = TRUE))
library(ggplot2)
pilot <- readRDS("data-raw/seasonal_pilot/pilot.rds")
models <- c("SS", "NT", "PE", "PM")
set.seed(1)

table_for <- function(years, data) {
  rows <- lapply(pilot$runs, function(r) {
    s <- r$stats[[years]]
    s <- if (data == "tadpoles") s[startsWith(names(s), "tad_")] else s
    data.frame(model = r$model, t(s))
  })
  do.call(rbind, rows)
}

settings <- expand.grid(years = 1:3, data = c("tadpoles", "tadpoles + calling males"),
                        stringsAsFactors = FALSE)
results <- lapply(seq_len(nrow(settings)), function(i) {
  tab <- table_for(settings$years[i], settings$data[i])
  four <- abc_choose(split(tab[-1], tab$model)[models], ntree = 500, ncores = 8)
  pair <- abc_choose(split(tab[-1], tab$model)[c("SS", "PM")], ntree = 500, ncores = 8)
  list(years = settings$years[i], data = settings$data[i], error = four$error,
       confusion = four$confusion, error_ss_pm = pair$error, confusion_ss_pm = pair$confusion,
       n = table(four$training$model))
})

summary_table <- data.frame(
  years = vapply(results, `[[`, 0, "years"),
  data = vapply(results, `[[`, "", "data"),
  error_four = vapply(results, `[[`, 0, "error"),
  error_ss_pm = vapply(results, `[[`, 0, "error_ss_pm"))
print(summary_table, digits = 3)
for (r in results) {
  cat("\n", r$years, "year(s),", r$data, "\n")
  print(round(prop.table(r$confusion, 1), 3))
}
saveRDS(list(summary = summary_table, results = results), "data-raw/seasonal_pilot/pilot_results.rds")

long <- rbind(data.frame(summary_table[1:2], comparison = "four models",
                         error = summary_table$error_four),
              data.frame(summary_table[1:2], comparison = "SS vs PM",
                         error = summary_table$error_ss_pm))
p <- ggplot(long, aes(years, error, colour = data)) + geom_line() + geom_point() +
  facet_wrap(~comparison) + scale_x_continuous(breaks = 1:3) +
  scale_colour_manual(values = c("#0072B2", "#E69F00")) + ylim(0, NA) +
  labs(x = "Field years sampled monthly", y = "Prior error rate (out of bag)", colour = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("data-raw/seasonal_pilot/pilot_error.png", p, width = 8, height = 3.8, dpi = 120,
       bg = "white")
