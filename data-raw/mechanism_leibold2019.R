# Why does the replication of Leibold et al. (2019) (data-raw/replication_leibold2019.R)
# depart from their model, with weaker local dominance and faster regional drift?
# Two hypotheses, crossed in a factorial design with the same seeds:
# (i) immigrants are merged into the mean of the resident population before selection, so
#     maladapted immigrants are not filtered one by one as in an individual-based model:
#     `immigrant_selection = TRUE` filters the immigrants of each source by the fitness of
#     their phenotype in the site they reach (see ?evolution_control);
# (ii) drift of the population means (and founder sampling) creates fitness differences
#     between populations that selection then amplifies, which would speed up regional drift
#     even between equivalent species (the control with one patch type): `ne_ratio = Inf`
#     turns drift of the means off (founder sampling stays); a last set also turns founder
#     sampling off (`founder = FALSE`).
# Everything else as in the replication.

source("data-raw/replication_leibold2019.R")

if (sys.nframe() == 0) {
  n_rep <- 10
  grid <- expand.grid(G = c(0, 50, 100, 200, 400), nm = c(1, 10, 100), landscape_type = "two",
                      selection = c(FALSE, TRUE), ne_ratio = c(0.5, Inf), rep = seq_len(n_rep),
                      stringsAsFactors = FALSE)
  control <- expand.grid(G = c(0, 200), nm = 10, landscape_type = "one",
                         selection = c(FALSE, TRUE), ne_ratio = c(0.5, Inf),
                         rep = seq_len(n_rep), stringsAsFactors = FALSE)
  design <- rbind(grid, control)
  # G = 0: no evolution, so drift and immigrant selection do not apply
  design <- design[!(design$G == 0 & (design$selection | is.infinite(design$ne_ratio))), ]
  t0 <- Sys.time()
  results <- parallel::mclapply(seq_len(nrow(design)), function(i) {
    d <- design[i, ]
    run <- run_leibold(d$G, d$landscape_type, seed = 1000 + d$rep, n_steps = 10000,
                       nm = d$nm, immigrant_selection = d$selection, ne_ratio = d$ne_ratio)
    list(G = d$G, nm = d$nm, landscape_type = d$landscape_type, selection = d$selection,
         ne_ratio = d$ne_ratio, rep = d$rep, summary = summarise_leibold(run),
         classes = classify_leibold(run))
  }, mc.cores = max(1, parallel::detectCores() - 2), mc.preschedule = FALSE)
  bad <- vapply(results, inherits, logical(1), "try-error")
  if (any(bad)) stop("runs failed: ", results[[which(bad)[1]]])
  message("time: ", format(Sys.time() - t0))
  saveRDS(results, "data-raw/mechanism_leibold2019.rds")

  # (iii) no drift of the means and no founder sampling (founder = FALSE), G = 200, Nm = 10
  des <- expand.grid(type = c("two", "one"), rep = seq_len(n_rep), stringsAsFactors = FALSE)
  out <- parallel::mclapply(seq_len(nrow(des)), function(i) {
    r <- run_leibold(200, des$type[i], seed = 1000 + des$rep[i], nm = 10, ne_ratio = Inf,
                     founder = FALSE)
    s <- summarise_leibold(r)
    data.frame(type = des$type[i], both5000 = s$both[s$step == 5000],
               both10000 = s$both[s$step == 10000],
               loss = if (all(s$both)) NA else s$step[which(!s$both)[1]],
               shared = if (s$both[s$step == 5000]) classify_leibold(r)[["shared"]] else NA)
  }, mc.cores = max(1, parallel::detectCores() - 2))
  saveRDS(do.call(rbind, out), "data-raw/mechanism_leibold2019_founder.rds")
}
