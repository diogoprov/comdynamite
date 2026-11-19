# Summaries of the block B runs (data-raw/blockB/blockB_*.R): the numbers of the three
# replication articles recomputed with 500 replicates, and the sweeps of their central
# parameters. Reads the chunks under OUT (see run_chunks.R), prints the numbers and writes the
# aggregated tables to data-raw/blockB/blockB_summary.rds (the raw chunks stay out of git).
# Run from the root of the package: OUT=<folder> Rscript data-raw/blockB/summarise_blockB.R

source("data-raw/blockB/run_chunks.R")
get <- function(name) collect_chunks(name)
ci <- function(v) {
  v <- v[!is.na(v)]
  m <- mean(v)
  se <- if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA
  c(mean = m, lo = m - 1.96 * se, hi = m + 1.96 * se, n = length(v))
}
agg <- function(formula, data, f = ci) {
  a <- aggregate(formula, data, f, na.action = na.pass)
  y <- as.data.frame(a[[ncol(a)]])
  cbind(a[-ncol(a)], y)
}
show <- function(title, x, digits = 3) {
  cat("\n## ", title, "\n", sep = "")
  print(x, row.names = FALSE, digits = digits)
}
out <- list()

# Leibold et al. (2019) --------------------------------------------------------------------
lr <- get("leibold_replicates")
lr$fitdiff <- rowMeans(abs(lr[grep("^fitdiff_e", names(lr))]), na.rm = TRUE)
# "first 10": the replicates of the replication article (same seeds), "all 500": block B
lr2 <- rbind(transform(lr[lr$rep <= 10, ], set = "first 10"), transform(lr, set = "all 500"))
# persistence of both species
pers <- aggregate(both ~ set + landscape_type + G + nm + step, lr2, mean)
out$leibold_persistence <- pers
show("Leibold: share of runs with both species (two patch types, Nm = 10, and control)",
     reshape(pers[pers$nm == 10, ], idvar = c("set", "landscape_type", "G", "nm"),
             timevar = "step", direction = "wide"))
show("Leibold: share of runs with both species at generation 10000, all Nm",
     xtabs(both ~ G + nm, aggregate(both ~ G + nm, lr[lr$step == 10000 &
                                                       lr$landscape_type == "two", ], mean)))
# signatures of regional neutrality at generation 5000, runs with both species
b5 <- lr2[lr2$step == 5000 & lr2$both, ]
sig <- merge(merge(agg(fitdiff ~ set + landscape_type + G + nm, b5)[, c(1:4, 5, 8)],
                   agg(env_comp ~ set + landscape_type + G + nm, b5)[, 1:5],
                   by = c("set", "landscape_type", "G", "nm"), suffixes = c("", ".comp")),
             agg(env_trait ~ set + landscape_type + G + nm, b5)[, 1:5],
             by = c("set", "landscape_type", "G", "nm"))
names(sig)[5:8] <- c("fitdiff", "n_both", "r2_comp", "r2_trait")
out$leibold_signatures <- sig
show("Leibold: fitness difference and adjusted R2 of patch type at generation 5000 (runs with both)",
     sig[sig$landscape_type == "two", ])
cls <- aggregate(cbind(preadapted, switched, shared) ~ set + landscape_type + G + nm, b5, mean)
cls$n <- aggregate(both ~ set + landscape_type + G + nm, b5, length)$both
out$leibold_classes <- cls
show("Leibold: classes of patches at generation 5000 (runs with both)", cls)

sp <- get("leibold_sweep_pext")
sp$fitdiff <- (abs(sp$fitdiff_e40) + abs(sp$fitdiff_e60)) / 2
sp_pers <- aggregate(both ~ G + nm + p_ext + step, sp, mean)
out$leibold_pext_persistence <- sp_pers
show("Leibold sweep: runs with both species at generation 10000, by G (rows), Nm and p_ext",
     ftable(xtabs(both ~ G + nm + p_ext, sp_pers[sp_pers$step == 10000, ])), digits = 2)
s5 <- sp[sp$step == 5000 & sp$both, ]
sp_sig <- aggregate(cbind(fitdiff, env_comp, env_trait, preadapted, switched, shared) ~
                      G + nm + p_ext, s5, mean)
sp_sig$n <- aggregate(both ~ G + nm + p_ext, s5, length)$both
out$leibold_pext_signatures <- sp_sig
# cells with fewer than 10 runs that kept both species are left empty
masked <- function(v) {
  x <- with(sp_sig, tapply(ifelse(n >= 10, base::get(v), NA), list(G = G, nm = nm, p_ext = p_ext), c))
  ftable(round(x, 3), row.vars = 1:2)
}
show("Leibold sweep: adjusted R2 of patch type for composition at generation 5000 (runs with both)",
     masked("env_comp"))
show("Leibold sweep: shared patches at generation 5000 (runs with both)", masked("shared"))
show("Leibold sweep: switched patches at generation 5000 (runs with both)", masked("switched"))

sv <- get("leibold_sweep_vm")
sv$fitdiff <- (abs(sv$fitdiff_e40) + abs(sv$fitdiff_e60)) / 2
sv_sig <- aggregate(cbind(both, fitdiff, env_comp, env_trait, preadapted, switched, shared) ~
                      vm + nm + step, sv, mean, na.action = na.pass)
out$leibold_vm <- sv_sig
show("Leibold sweep of Vm (from no variance): generation 5000 and 10000",
     sv_sig[sv_sig$step %in% c(5000, 10000), ])

# Urban & De Meester (2009) ----------------------------------------------------------------
f2 <- get("urban_fig2")
f2$set <- "all 500"
f2_10 <- f2[f2$rep <= 100, ]
f2_10$set <- "first 100"
f2 <- rbind(f2_10, f2)
f2$first_wins <- ifelse(f2$delay > 0, f2$n1 > f2$n2, f2$n2 > f2$n1)
f2$lead <- abs(f2$delay)
win <- agg(first_wins ~ set + genetic + lead, f2[f2$delay != 0 & abs(f2$delay) %in%
                                                     c(1, 5, 10, 20), ])
out$urban_first_wins <- win
show("Urban Fig. 2: share of runs won by the first species to enter", win)
f2$share1 <- f2$n1 / pmax(f2$n1 + f2$n2, 1)
sh <- agg(share1 ~ set + genetic + delay, f2[f2$delay %in% c(-10, 0, 10), ])
out$urban_share1 <- sh
show("Urban Fig. 2: share of species 1 in the new site, delays -10, 0 and 10", sh)
f2$mono <- pmax(f2$n1, f2$n2) / pmax(f2$n1 + f2$n2, 1) > 0.9
mono <- agg(mono ~ set + genetic, f2)
out$urban_mono <- mono
show("Urban Fig. 2: runs in which one species holds more than 90% of the new site", mono)

f4 <- get("urban_fig4")
f4$dom <- (f4$n1 - f4$n2) / 1000
m4 <- agg(dom ~ m1 + m2, f4)
out$urban_dispersal <- m4
show("Urban Fig. 4: dominance (N1 - N2) / 1000, mean over 500 runs",
     round(xtabs(mean ~ m1 + m2, m4), 2))
m4$better <- ifelse(m4$m1 > m4$m2, 1, ifelse(m4$m2 > m4$m1, 2, 0))
m4$better_wins <- ifelse(m4$better == 1, m4$lo > 0, ifelse(m4$better == 2, m4$hi < 0, NA))
m4$poorer_wins <- ifelse(m4$better == 1, m4$hi < 0, ifelse(m4$better == 2, m4$lo > 0, NA))
cat("\nOff-diagonal cells (of 110) where the better disperser wins (95% CI excludes 0):",
    sum(m4$better_wins, na.rm = TRUE), "; where the poorer disperser wins:",
    sum(m4$poorer_wins, na.rm = TRUE), "\n")
show("Urban Fig. 4: cells where the better disperser wins", m4[m4$better_wins %in% TRUE, ])

f3 <- get("urban_fig3")
last <- f3[f3$step == 2000, ]
last$winner <- ifelse(last$new1 > last$new2, 1, 2)
last$own_share <- ifelse(last$winner == 1, last$src1_sp1 / pmax(last$src1_sp1 + last$src1_sp2, 1),
                         last$src2_sp2 / pmax(last$src2_sp1 + last$src2_sp2, 1))
last$loser_share <- ifelse(last$winner == 1, last$src2_sp2 / pmax(last$src2_sp1 + last$src2_sp2, 1),
                           last$src1_sp1 / pmax(last$src1_sp1 + last$src1_sp2, 1))
src <- merge(agg(own_share ~ m, last), agg(loser_share ~ m, last), by = "m",
             suffixes = c(".winner", ".loser"))
out$urban_sources <- src
show("Urban Fig. 3: share of each species in its own source at 2000 generations", src)

us <- get("urban_sweep")
us$first_wins <- ifelse(us$delay > 0, us$n1 > us$n2, us$n2 > us$n1)
us$share1 <- us$n1 / pmax(us$n1 + us$n2, 1)
us$mono <- pmax(us$n1, us$n2) / pmax(us$n1 + us$n2, 1) > 0.9
# species 1 alone evolving: its share when it enters 10 generations after species 2
s1 <- agg(share1 ~ genetic + vm + env_new, us[us$genetic %in% c("mutation_sp1", "standing_sp1") &
                                                us$delay == -10, ])
out$urban_sweep_sp1_late <- s1
show("Urban sweep: share of species 1 (alone evolving) when it enters 10 generations later",
     ftable(xtabs(mean ~ genetic + vm + env_new, s1)), digits = 2)
s1e <- agg(share1 ~ genetic + vm + env_new, us[us$genetic %in% c("mutation_sp1", "standing_sp1") &
                                                 us$delay == 10, ])
out$urban_sweep_sp1_early <- s1e
show("Urban sweep: share of species 1 (alone evolving) when it enters 10 generations earlier",
     ftable(xtabs(mean ~ genetic + vm + env_new, s1e)), digits = 2)
fw <- agg(first_wins ~ genetic + vm + env_new, us[us$genetic %in% c("mutation_both",
                                                                    "standing_both") &
                                                    abs(us$delay) == 10, ])
out$urban_sweep_first_wins <- fw
show("Urban sweep: first arrival wins, lead of 10 generations, both species evolving",
     ftable(xtabs(mean ~ genetic + vm + env_new, fw)), digits = 2)
mo <- agg(mono ~ genetic + vm + env_new, us[us$genetic %in% c("mutation_both", "standing_both"), ])
out$urban_sweep_mono <- mo
show("Urban sweep: one species above 90% of the new site, both evolving, all delays",
     ftable(xtabs(mean ~ genetic + vm + env_new, mo)), digits = 2)

# Hajian-Forooshani & Chase (2026) ---------------------------------------------------------
hr <- get("hajian_replicates")
hr$d_beta <- hr$d_gamma - hr$d_alpha
h700 <- hr[hr$step == 700, ]
keys <- c("niche", "dispersal", "n_disturbed")
ha_a <- agg(d_alpha ~ niche + dispersal + n_disturbed, h700)
ha_g <- agg(d_gamma ~ niche + dispersal + n_disturbed, h700)
ha_b <- agg(d_beta ~ niche + dispersal + n_disturbed, h700)
hb <- data.frame(ha_b[keys], d_alpha = ha_a$mean, d_gamma = ha_g$mean, d_beta = ha_b$mean,
                 lo_beta = ha_b$lo, hi_beta = ha_b$hi, n = ha_b$n)
out$hajian_step700 <- hb
show("Hajian: mean change in beta diversity at step 700, niche 0.25 (rows: extent)",
     xtabs(d_beta ~ n_disturbed + dispersal, hb[hb$niche == 0.25, ]))
ctrl <- hb[hb$n_disturbed == 0, c("niche", "dispersal", "d_beta", "lo_beta", "hi_beta")]
show("Hajian: change in beta diversity without disturbance (control), step 700", ctrl)
hb$excess <- hb$d_beta - ctrl$d_beta[match(paste(hb$niche, hb$dispersal),
                                           paste(ctrl$niche, ctrl$dispersal))]
show("Hajian: beta change minus the control, niche 0.25, step 700",
     xtabs(excess ~ n_disturbed + dispersal, hb[hb$niche == 0.25 & hb$n_disturbed > 0, ]))
# sign of the change at intermediate dispersal, with the 95% CI
mid <- hb[hb$niche == 0.25 & hb$dispersal == 1e-4, ]
mid$sign <- ifelse(mid$lo_beta > 0, "+", ifelse(mid$hi_beta < 0, "-", "0"))
show("Hajian: sign of the beta change at dispersal 1e-4, niche 0.25 (95% CI)",
     mid[, c("n_disturbed", "d_beta", "lo_beta", "hi_beta", "sign")])

hs <- get("hajian_sweep")
hs$d_beta <- hs$d_gamma - hs$d_alpha
ha <- agg(d_beta ~ when + duration + value + dispersal + n_disturbed, hs)
ha$sign <- ifelse(ha$lo > 0, "+", ifelse(ha$hi < 0, "-", "0"))
out$hajian_sweep <- ha
for (d in sort(unique(ha$dispersal))) {
  show(sprintf("Hajian sweep: beta change 100 steps after the pulse, dispersal %g, extent 40 (rows: duration)", d),
       round(xtabs(mean ~ duration + value, ha[ha$when == "after" & ha$dispersal == d &
                                                ha$n_disturbed == 40, ]), 2))
}
show("Hajian sweep: settings by sign of the beta change 100 steps after the pulse, by dispersal",
     table(ha$dispersal[ha$when == "after"], ha$sign[ha$when == "after"]))
show("Hajian sweep: same, pulses with value 0.1 to 0.3 only",
     with(ha[ha$when == "after" & ha$value <= 0.3, ], table(dispersal, sign)))

saveRDS(out, "data-raw/blockB/blockB_summary.rds")
