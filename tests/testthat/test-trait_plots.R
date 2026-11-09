toy_plots <- function() {
  set.seed(1)
  ind <- data.frame(site = rep(paste0("s", 1:6), each = 20),
                    species = sample(c("a", "b", "c"), 120, replace = TRUE))
  ind$size <- c(a = 10, b = 20, c = 30)[ind$species] + rnorm(120, sd = 3)
  rbind(ind, data.frame(site = "solo", species = "a", size = rnorm(10, 10, 3)))
}

test_that("null_intervals() summarises the nulls of tstats() and community_index()", {
  ind <- toy_plots()
  res <- tstats(ind, traits = "size", n_perm = 99)
  d <- null_intervals(res, level = 0.9, inner = 0.5)
  expect_equal(nrow(d), 3 * 7)
  null <- attr(res, "null")
  row <- d[d$site == "s1" & d$statistic == "T_IC.IR", ]
  v <- null[, 1, "T_IC.IR"]
  expect_equal(row$null_low, unname(stats::quantile(v, 0.05)))
  expect_equal(row$inner_high, unname(stats::quantile(v, 0.75)))
  expect_equal(row$null_median, stats::median(v))
  # deviation agrees with tstats_summary() at the same size
  s <- tstats_summary(res, level = 0.9)$sites
  expect_equal(d$deviation, s$deviation[match(paste(d$site, d$statistic),
                                              paste(s$site, s$statistic))])
  expect_true(is.na(d$deviation[d$site == "solo" & d$statistic == "T_IP.IC"]))
  ci <- community_index(ind, list(variance = stats::var, MNND = mnnd), traits = "size",
                        null = "regional.ind", n_perm = 19)
  dc <- null_intervals(ci)
  expect_equal(unique(dc$statistic), c("variance", "MNND"))
  expect_error(null_intervals(tstats(ind, traits = "size", n_perm = 0)), "n_perm > 0")
  expect_error(null_intervals(res, level = 0.5, inner = 0.6))
})

test_that("plot() of null tests draws intervals, effect sizes and hatching", {
  skip_if_not_installed("ggplot2")
  ind <- toy_plots()
  res <- tstats(ind, traits = "size", n_perm = 99)
  p <- plot(res)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  points <- ggplot2::layer_data(p, length(p$layers) - 1)
  expect_equal(nrow(points), sum(!is.na(null_intervals(res)$deviation)))
  # hatching lines stay inside their rectangles
  rects <- ggplot2::layer_data(p, 2)
  lines <- ggplot2::layer_data(p, 3)
  expect_gt(nrow(lines), 0)
  inside <- vapply(seq_len(nrow(lines)), function(i) {
    l <- lines[i, ]
    any(rects$PANEL == l$PANEL & l$x >= rects$xmin - 1e-9 & l$xend <= rects$xmax + 1e-9 &
          l$y >= rects$ymin - 1e-9 & l$yend <= rects$ymax + 1e-9)
  }, logical(1))
  expect_true(all(inside))
  ses <- plot(res, type = "ses", order = "ses", hatch = FALSE, statistics = "T_IC.IR")
  expect_s3_class(ses, "ggplot")
  expect_equal(levels(ses$data$statistic), "T_IC.IR")
  by_observed <- levels(plot(res, order = "observed")$data$site)
  key <- tapply(null_intervals(res)$observed, null_intervals(res)$site, mean, na.rm = TRUE)
  expect_equal(by_observed, names(sort(key)))
  expect_error(plot(res, sites = "nowhere"), "not in `x`")
  expect_error(plot(res, colours = "red"), "3 valid colours")
  ci <- community_index(ind, stats::var, traits = "size", null = "regional.ind", n_perm = 19)
  expect_s3_class(plot(ci, type = "ses"), "ggplot")
  expect_equal(levels(plot(ci)$data$statistic), "stats::var")
})

test_that("trait_densities() integrates to one and sums populations into communities", {
  ind <- toy_plots()
  ind <- rbind(ind, data.frame(site = "s1", species = "z", size = 15))
  d <- trait_densities(ind, "size", n = 1024)
  step <- diff(d$trait[1:2])
  areas <- tapply(d$density, paste(d$site, d$species, d$level), sum) * step
  expect_true(all(abs(areas - 1) < 1e-3))
  s1 <- d[d$site == "s1", ]
  pops <- s1[s1$level == "population", ]
  community <- tapply(pops$density * pops$weight, pops$trait, sum)
  expect_equal(as.vector(community), s1$density[s1$level == "community"])
  expect_equal(sum(unique(pops[c("species", "weight")])$weight), 1)
  expect_gt(attr(d, "bw"), 0)
  expect_equal(attr(trait_densities(ind, "size", bw = 0.7), "bw"), 0.7)
  ind$abundance <- c(a = 1, b = 2, c = 7, z = 1)[ind$species]
  w <- trait_densities(ind, "size", abundance = "abundance")
  w <- unique(w[w$level == "population", c("site", "species", "weight")])
  expect_equal(w$weight[w$site == "s2" & w$species == "c"],
               7 / sum(c(a = 1, b = 2, c = 7)[unique(ind$species[ind$site == "s2"])]))
  expect_error(trait_densities(ind[!duplicated(paste(ind$site, ind$species)), ], "size"),
               "give `bw`")
  expect_error(trait_densities(ind, "nope"), "no column")
})

test_that("plot_population_traits() draws populations, communities and rugs", {
  skip_if_not_installed("ggplot2")
  ind <- toy_plots()
  p <- plot_population_traits(ind, "size", alpha = 0.2)
  expect_s3_class(p, "ggplot")
  areas <- ggplot2::layer_data(p, 1)
  expect_equal(unique(areas$alpha), 0.2)
  expect_length(p$layers, 3)
  bare <- plot_population_traits(ind, "size", community = FALSE, rug = FALSE,
                                 weighted = FALSE, sites = c("s1", "s2"))
  expect_length(bare$layers, 1)
  expect_equal(levels(bare$data$site), c("s1", "s2"))
  named <- plot_population_traits(ind, "size", colours = c(c = "black", a = "red", b = "blue"))
  expect_s3_class(named, "ggplot")
  expect_error(plot_population_traits(ind, "size", colours = c(a = "red")), "no colour")
  expect_error(plot_population_traits(ind, "size", sites = "nowhere"), "not in the data")
  many <- ind
  many$species <- paste0(many$species, seq_len(nrow(many)) %% 9)
  expect_s3_class(plot_population_traits(many, "size"), "ggplot")
})

test_that("the default bandwidth is not shrunk by small populations", {
  set.seed(3)
  bws <- replicate(200, {
    ind <- data.frame(site = rep(paste0("s", 1:20), each = 2), species = "a",
                      trait = rnorm(40, rep(rnorm(20, 0, 5), each = 2)))
    attr(trait_densities(ind, n = 16), "bw")
  })
  # Silverman's rule for sd 1 and n = 2
  expect_equal(mean(bws), 0.9 * 2^(-1 / 5), tolerance = 0.12)
  ind <- data.frame(site = "s1", species = rep(c("a", "b"), each = 3), trait = rep(1:2, each = 3))
  expect_error(trait_densities(ind), "identical values")
})

test_that("plots of all-NA statistics or untestable sites give no warnings", {
  skip_if_not_installed("ggplot2")
  set.seed(4)
  ind <- data.frame(site = rep(paste0("s", 1:4), each = 10), species = rep(paste0("p", 1:4), each = 10),
                    trait = rnorm(40))
  res <- suppressWarnings(tstats(ind, n_perm = 19))
  expect_no_warning(ggplot2::ggplot_build(plot(res)))
  expect_no_warning(ggplot2::ggplot_build(plot(res, type = "ses")))
})

test_that("plot options and input checks of trait plots", {
  skip_if_not_installed("ggplot2")
  set.seed(5)
  ind <- data.frame(site = rep(paste0("s", 1:3), each = 12), species = rep(c("a", "b"), 18),
                    trait = rnorm(36), ab = 1)
  res <- tstats(ind, n_perm = 19)
  expect_s3_class(plot(res, hatch = FALSE), "ggplot")
  expect_equal(nrow(hatch_lines(data.frame(xmin = numeric(0), xmax = numeric(0),
                                           ymin = numeric(0), ymax = numeric(0)), 1)), 0)
  bad <- ind
  bad$trait <- as.character(bad$trait)
  expect_error(trait_densities(bad), "numeric")
  bad <- ind
  bad$ab[1] <- 2
  expect_error(trait_densities(bad, abundance = "ab"), "one value")
})

test_that("community_index rejects malformed lists of indices", {
  ind <- data.frame(site = "s1", species = c("a", "b"), trait = 1:2)
  expect_error(community_index(ind, list(mean, var)), "distinct names")
})

test_that("trait_densities warns when the spread rises with the mean", {
  set.seed(7)
  sizes <- exp(seq(log(15), log(200), length.out = 12))
  ind <- data.frame(site = rep(c("s1", "s2"), each = 60), species = rep(paste0("sp", 1:12), 10))
  ind$crc <- sizes[as.integer(sub("sp", "", ind$species))] * exp(rnorm(120, 0, 0.1))
  expect_warning(trait_densities(ind, "crc", n = 16), "logarithm")
  ind$log_crc <- log(ind$crc)
  expect_no_warning(trait_densities(ind, "log_crc", n = 16))
  expect_no_warning(trait_densities(ind, "crc", bw = 1, n = 16))
  # equal spread on a positive scale does not warn
  ind$even <- 100 + sizes[as.integer(sub("sp", "", ind$species))] + rnorm(120)
  expect_no_warning(trait_densities(ind, "even", n = 16))
})

test_that("plot_population_traits highlights species", {
  skip_if_not_installed("ggplot2")
  set.seed(8)
  ind <- data.frame(site = rep(c("s1", "s2"), each = 40), species = rep(letters[1:8], 10),
                    trait = rnorm(80, rep(1:8, 10)))
  p <- plot_population_traits(ind, highlight = c("c", "f"), n = 32)
  b <- ggplot2::ggplot_build(p)
  fills <- unique(b$data[[1]][c("group", "fill")])
  expect_setequal(unique(fills$fill), c("#0072B2", "#E69F00", "grey75"))
  expect_equal(as.vector(b$plot$scales$get_scales("fill")$get_breaks()), c("c", "f"))
  # highlighted populations are drawn after the grey ones in each panel
  d1 <- b$data[[1]]
  expect_equal(tail(unique(d1$fill[d1$PANEL == 1]), 2), c("#0072B2", "#E69F00"))
  p2 <- plot_population_traits(ind, highlight = "a", colours = c(a = "red"), n = 32)
  expect_true("red" %in% ggplot2::ggplot_build(p2)$data[[1]]$fill)
  expect_error(plot_population_traits(ind, highlight = "z"), "not in the data")
  expect_error(plot_population_traits(ind, highlight = "a", colours = c(b = "red")),
               "no colour")
  expect_error(plot_population_traits(ind, highlight = letters[1:8]), "seven")
})

test_that("the scale check skips populations with identical spreads", {
  ind <- data.frame(site = "s1", species = rep(paste0("sp", 1:10), each = 3))
  ind$trait <- rep(seq(10, 100, by = 10), each = 3) + rep(c(-1, 0, 1), 10)
  expect_no_warning(trait_densities(ind, n = 16))
  expect_null(check_spread_scale(ind$trait, ind$species))
})
