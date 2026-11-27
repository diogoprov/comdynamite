# Validation of trait_change(), item (i) of block K of dev/ROADMAP.md: the eco-evolutionary
# partition of trait change over time. Simulations with one evolving trait on a 6 x 6 grid with
# eight species, 200 steps recorded every 50, in a warming environment (+0.1 per step) or a
# static one, over a grid of genetic variance G (0, 1, 4) and plasticity b (0, 0.3, 0.6), 20
# replicates each. Criterion set before running:
#   (1) the community components add up to the change of the community-weighted mean, and the
#       genetic and plastic parts to the change within species, with error < 1e-8; the
#       components of the growth rate add up to its change with error < 1e-8;
#   (2) the genetic parts are 0 (|x| < 1e-10) with G = 0 and the plastic parts with b = 0;
#   (3) under warming, the mean genetic part (community level) increases with G at every b > 0
#       or b = 0, and the mean plastic part increases with b at every G;
#   (4) in the static environment, the plastic part and the direct effect of the environment
#       on growth are 0 (|x| < 1e-10).
# A further check, not part of the criterion: as the log growth rate is quadratic in phenotype
# and environment, the Shapley values should equal the gradient at the midpoint of the two
# steps times the change of each component, i.e. eqn 6 of Ellner et al. (2011) integrated
# along the straight path between the steps.
# Run from the root of the package: Rscript data-raw/trait_change.R. REPS, CORES and OUTPUT
# change the run.
suppressMessages(devtools::load_all(".", quiet = TRUE))
options(comdynamite.cite = FALSE)
reps <- as.integer(Sys.getenv("REPS", "20"))
cores <- as.integer(Sys.getenv("CORES", "8"))
output <- Sys.getenv("OUTPUT", "data-raw/trait_change.rds")
n_steps <- 200
record <- seq(50, 200, 50)

one_run <- function(g, b, environment, rep) {
  set.seed(1000 * rep + 7)
  ls <- landscape_grid(6)
  tr <- ape::rcoal(8)
  rate <- if (environment == "warming") 0.1 else 0
  env_time <- outer(ls$env, rate * seq_len(n_steps), "+")
  # the same community in every setting of a replicate
  traits <- species_traits(tr)
  alpha <- competition_matrix(tr, "phylo")
  set.seed(1000 * rep + 8)
  sim <- simulate_metacommunity(ls, traits, alpha, n_steps = n_steps, record_steps = record,
                                env_time = env_time,
                                evolution = evolution_control(g, plasticity = b))
  res <- trait_change(sim)
  com <- res$community
  gr <- res$growth
  ok <- !is.na(com$total)
  # gradient at the midpoint: X = c - (E - p)^2 / (2 w), so dX/dp = (E - p) / w = -dX/dE
  mid <- midpoint_terms(sim, gr, b)
  s <- summary(res)
  data.frame(
    g = g, b = b, environment = environment, rep = rep,
    closure_community = max(abs(com$total - com$sorting - com$turnover - com$within -
                                  com$covariance)[ok]),
    closure_within = max(abs(com$within - com$genetic - com$plastic)[ok]),
    closure_growth = if (nrow(gr)) max(abs(gr$total - gr$genetic - gr$plastic -
                                             gr$environment)) else 0,
    max_genetic = max(abs(c(com$genetic[ok], com$cov_genetic[ok], gr$genetic))),
    max_plastic = max(abs(c(com$plastic[ok], com$cov_plastic[ok], gr$plastic))),
    max_environment = if (nrow(gr)) max(abs(gr$environment)) else 0,
    midpoint_error = mid,
    total = mean(s$community$total), sorting = mean(s$community$sorting),
    turnover = mean(s$community$turnover), genetic = mean(s$community$genetic),
    plastic = mean(s$community$plastic),
    covariance = mean(s$community$cov_genetic + s$community$cov_plastic),
    growth_total = mean(gr$total), growth_genetic = mean(gr$genetic),
    growth_plastic = mean(gr$plastic), growth_environment = mean(gr$environment),
    populations = nrow(gr)
  )
}

midpoint_terms <- function(sim, gr, b) {
  if (!nrow(gr)) return(0)
  ev <- sim$evolution
  err <- vapply(seq_len(nrow(gr)), function(i) {
    r <- gr[i, ]
    at <- function(s) {
      lab <- paste0("step", s)
      e <- ev$environment_history[r$site, lab]
      c(g = ev$history[r$site, r$species, lab], q = b * (e - ev$reference), e = e)
    }
    a <- at(r$from)
    z <- at(r$to)
    m <- (a + z) / 2
    sigma <- sim$traits[r$species, "sigma"]
    v <- ev$settings$genetic_var + ev$settings$env_var
    width <- sigma^2 + v
    slope <- (m[["e"]] - m[["g"]] - m[["q"]]) / width
    est <- c(slope * (z[["g"]] - a[["g"]]), slope * (z[["q"]] - a[["q"]]),
             -slope * (z[["e"]] - a[["e"]]))
    max(abs(est - c(r$genetic, r$plastic, r$environment)))
  }, numeric(1))
  max(err)
}

if (sys.nframe() == 0) {
  started <- Sys.time()
  grid <- expand.grid(g = c(0, 1, 4), b = c(0, 0.3, 0.6),
                      environment = c("warming", "static"), rep = seq_len(reps),
                      stringsAsFactors = FALSE)
  runs <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
    with(grid[i, ], one_run(g, b, environment, rep))
  }, mc.cores = cores)
  failed <- vapply(runs, inherits, logical(1), "try-error")
  if (any(failed)) stop(runs[[which(failed)[1]]])
  runs <- do.call(rbind, runs)
  means <- stats::aggregate(cbind(total, sorting, turnover, genetic, plastic, covariance,
                                  growth_total, growth_genetic, growth_plastic,
                                  growth_environment) ~ environment + g + b,
                            data = runs, FUN = mean)
  warm <- means[means$environment == "warming", ]
  increasing <- function(x) all(diff(x) > 0)
  criteria <- c(
    closure = max(runs$closure_community, runs$closure_within, runs$closure_growth) < 1e-8,
    genetic_zero = max(runs$max_genetic[runs$g == 0]) < 1e-10,
    plastic_zero = max(runs$max_plastic[runs$b == 0]) < 1e-10,
    genetic_increases = all(vapply(split(warm, warm$b), function(d) {
      increasing(d$genetic[order(d$g)])
    }, logical(1))),
    plastic_increases = all(vapply(split(warm, warm$g), function(d) {
      increasing(d$plastic[order(d$b)])
    }, logical(1))),
    static_zero = max(runs$max_plastic[runs$environment == "static"],
                      runs$max_environment[runs$environment == "static"]) < 1e-10
  )
  minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  saveRDS(list(runs = runs, means = means, criteria = criteria,
               midpoint_error = max(runs$midpoint_error), minutes = minutes,
               settings = list(reps = reps, n_steps = n_steps, record = record),
               date = Sys.Date()), output, compress = "xz")
  print(means, digits = 3, row.names = FALSE)
  print(criteria)
  cat("largest difference from the midpoint gradient:", max(runs$midpoint_error), "\n")
  cat(round(minutes, 1), "min\n")
}
