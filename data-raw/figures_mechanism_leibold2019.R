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

