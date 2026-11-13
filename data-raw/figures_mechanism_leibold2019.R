# Table of the factorial experiment of data-raw/mechanism_leibold2019.R (numbers in the article
# "Replicating Leibold et al. (2019)", section Two explanations tested).
res <- readRDS("data-raw/mechanism_leibold2019.rds")
row <- function(r) {
  s <- r$summary
  b5 <- s$both[s$step == 5000]; b10 <- s$both[s$step == 10000]
  loss <- if (all(s$both)) NA else s$step[which(!s$both)[1]]
  late <- s[s$step > 1000 & s$step <= 5000 & s$both, ]
  data.frame(G = r$G, nm = r$nm, type = r$landscape_type, sel = r$selection, ne = r$ne_ratio, rep = r$rep,
             both5000 = b5, both10000 = b10, loss = loss,
             shared = if (b5) r$classes[["shared"]] else NA,
             switched = if (b5) r$classes[["switched"]] else NA,
             fitdiff = if (nrow(late) && r$landscape_type == "two") mean((abs(late$fitdiff_e40) + abs(late$fitdiff_e60))/2) else NA,
             env_comp = if (nrow(late)) mean(late$env_comp, na.rm = TRUE) else NA)
}
d <- do.call(rbind, lapply(res, row))
agg <- function(x) aggregate(cbind(both5000, both10000, shared, switched, fitdiff, env_comp) ~ G + nm + type + sel + ne,
  x, function(v) round(mean(v, na.rm = TRUE), 2), na.action = na.pass)
a <- agg(d)
lt <- aggregate(loss ~ G + nm + type + sel + ne, d, function(v) median(v), na.action = na.omit)
a <- merge(a, lt, all.x = TRUE)
nb <- aggregate(both5000 ~ G + nm + type + sel + ne, d, sum); names(nb)[6] <- "n_both5000"
a <- merge(a, nb)
options(width = 200)
print(a[order(a$type, a$nm, a$G, a$sel, a$ne), ], row.names = FALSE)


# genetic variance that depends on population size (block iv)
res_v <- readRDS("data-raw/mechanism_leibold2019_variance.rds")
row_v <- function(r) {
  s <- r$summary
  b5 <- s$both[s$step == 5000]
  late <- s[s$step > 1000 & s$step <= 5000 & s$both, ]
  data.frame(G = r$G, nm = r$nm, type = r$landscape_type, sel = r$selection, rep = r$rep,
             both5000 = b5, both10000 = s$both[s$step == 10000],
             loss = if (all(s$both)) NA else s$step[which(!s$both)[1]],
             shared = if (b5) r$classes[["shared"]] else NA,
             switched = if (b5) r$classes[["switched"]] else NA,
             fitdiff = if (nrow(late) && r$landscape_type == "two") mean((abs(late$fitdiff_e40) + abs(late$fitdiff_e60))/2) else NA,
             env_comp = if (nrow(late)) mean(late$env_comp, na.rm = TRUE) else NA,
             v_min = r$v_minority, v_maj = r$v_majority)
}
dv <- do.call(rbind, lapply(res_v, row_v))
av <- aggregate(cbind(both5000, both10000, shared, switched, fitdiff, env_comp, v_min, v_maj) ~ G + nm + type + sel, dv,
               function(v) round(mean(v, na.rm = TRUE), 2), na.action = na.pass)
lt <- NULL
print(av[order(av$type, av$nm, av$G, av$sel), ], row.names = FALSE)

# figure: constant against size-dependent genetic variance, without immigrant selection
library(ggplot2)
pick <- function(x, label) {
  x <- x[x$type == "two" & !x$sel & x$G > 0, ]
  if ("ne" %in% names(x)) x <- x[x$ne == 0.5, ]
  data.frame(G = x$G, nm = x$nm, variance = label, both10000 = x$both10000, shared = x$shared)
}
cmp <- rbind(pick(d, "constant G"), pick(dv, "size-dependent variance"))
part <- function(v, what) {
  out <- aggregate(cmp[[v]], cmp[c("G", "nm", "variance")], mean, na.rm = TRUE)
  names(out)[4] <- "value"
  out$what <- what
  out
}
cmp_l <- rbind(part("both10000", "Both species at 10,000"),
               part("shared", "Shared patches at 5000"))
cmp_l$nm_f <- factor(paste("Nm =", cmp_l$nm), paste("Nm =", c(1, 10, 100)))
p <- ggplot(cmp_l, aes(factor(G), value, colour = variance, group = variance)) +
  geom_point(size = 2) + geom_line() + facet_grid(what ~ nm_f) +
  scale_colour_manual(values = c("constant G" = "#D55E00", "size-dependent variance" = "#0072B2")) +
  labs(x = "G", y = NULL, colour = NULL,
       caption = "Shared patches: mean over the runs with both species at generation 5000") +
  theme_minimal(base_size = 11)
ggsave("vignettes/articles/leibold_variance.png", p, width = 9, height = 4.5, dpi = 120, bg = "white")
