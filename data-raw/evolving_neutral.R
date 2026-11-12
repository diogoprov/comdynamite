# Validation 1 of speciation in simulate_evolving_metacommunity(): the neutral limit. One
# community (a second site, empty and isolated, only because landscapes need two) of about J
# individuals, all species equivalent (competition equal within and between species, a niche so
# wide that growth does not depend on the environment, no genetic variance), and speciation by
# point mutation with probability nu per individual per generation. Under the Beverton-Holt
# dynamics with Poisson offspring, all individuals are born each generation, so the community
# behaves as a haploid Wright-Fisher population of size J with mutation rate nu (our
# interpretation; the census size fluctuates around J). Expected richness at equilibrium is
# then that of the Ewens sampling formula (Ewens 1972), E[S] = sum_{i=0}^{J-1} theta / (theta +
# i), with theta = 2 J nu, and the expected number of singletons is theta J / (J + theta - 1).
# Richness (extant lineages, all of them species without protracted speciation) is averaged
# over the second half of 4 J generations; census size and singletons are those at the end (the
# records of the output keep only the species extant at the end). Run from the root of the
# package; REPS changes the size.
# Result (20 replicates, 2026-11-12): richness matched the Ewens expectation within 1-4% for
# theta from 1 to 10 and J of 250 and 500 (ratio 0.99 to 1.04, standard errors 0.2 to 0.3
# species). The excess grows with theta (41.3 against 39.7 with theta = 10 and J = 500, about
# six standard errors), perhaps because the census size fluctuates or offspring numbers are
# not exactly those of Wright-Fisher; we did not pin it down. Singletons at the end, one snapshot
# per replicate, were 10% fewer than expected with theta = 1 and 15-38% more with theta of 2 or
# more.
# Ewens, W. J. (1972). The sampling theory of selectively neutral alleles. Theoretical Population
# Biology, 3, 87-112. https://doi.org/10.1016/0040-5809(72)90035-4
suppressMessages(devtools::load_all(".", quiet = TRUE))
reps <- as.integer(Sys.getenv("REPS", "20"))
ancestor <- data.frame(optimum = 50, sigma = 1e4, r_max = 5, row.names = "sp1")
grid <- expand.grid(J = c(250, 500), theta = c(1, 2, 5, 10), rep = seq_len(reps))
ewens <- function(theta, J) sum(theta / (theta + 0:(J - 1)))
run <- function(i) {
  g <- grid[i, ]
  set.seed(5000 + i)
  ls <- landscape_implicit(2)
  n_steps <- 4 * g$J
  sim <- simulate_evolving_metacommunity(
    ls, ancestor, speciation_control("point", rate = g$theta / (2 * g$J)), dispersal = 0,
    n_steps = n_steps, inter_min = 1, inter_max = 1, scale = 4 / g$J,
    initial = matrix(c(g$J, 0L), 2, 1, dimnames = list(NULL, "sp1")), max_species = 1e6)
  size <- sum(sim$abundance)
  S <- mean(sim$richness[(n_steps / 2 + 1):n_steps])
  singletons <- sum(colSums(sim$abundance) == 1)
  data.frame(g, size = size, S = S, S_expected = ewens(g$theta, round(size)),
             singletons = singletons,
             singletons_expected = g$theta * size / (size + g$theta - 1),
             lineages = nrow(sim$species))
}
started <- Sys.time()
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), run,
                                         mc.cores = as.integer(Sys.getenv("CORES", "6"))))
summary_table <- stats::aggregate(cbind(size, S, S_expected, singletons, singletons_expected,
                                        lineages) ~ J + theta, res,
                                  function(x) round(mean(x), 2))
summary_table$S_ratio <- round(summary_table$S / summary_table$S_expected, 3)
summary_table$S_se <- stats::aggregate(S ~ J + theta, res,
                                       function(x) round(stats::sd(x) / sqrt(length(x)), 2))$S
saveRDS(list(results = res, summary = summary_table, reps = reps,
             minutes = as.numeric(difftime(Sys.time(), started, units = "mins"))),
        "data-raw/evolving_neutral.rds")
print(summary_table)
