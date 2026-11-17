# Summaries of the simulations archived by Hajian-Forooshani & Chase (2026), for a quantitative
# comparison with data-raw/replication_hajian2026.rds. Their file
# 2025_03_18_metacom_disturb_sims_vA.RDS (106 MB, Figshare doi:10.6084/m9.figshare.30620603.v1,
# CC BY 4.0) is not kept in the package: download it and give its path in HFC_RDS.
# It has one row per replicate and step after the pulse (601 to 1000), 100 replicates for each
# of 7 niche breadths (0.1 to 0.4), 4 dispersal rates (1e-5 to 1e-2), 21 extents and 2 dispersal
# kernels (exp(-0.1 d), as in the article, and exp(-0.5 d)). The script of their Fig. 2 keeps
# the first kernel (0.1); that of their Fig. 4 overrides its kernel filter and pools both.
# Here both are summarised separately; the article uses kernel 0.1.
# Run from the root of the package; OUTPUT changes the file.
path <- Sys.getenv("HFC_RDS")
if (!nzchar(path) || !file.exists(path)) stop("Set HFC_RDS to the path of their RDS file.")
output <- Sys.getenv("OUTPUT", "data-raw/hajian2026_archived_summary.rds")

w <- readRDS(path)
w <- w[w$dispersal_prob %in% c(1e-5, 1e-4, 1e-3), ]
w$d_beta <- w$delta.gamma - w$delta.alpha
# a replicate is identified by its values before the pulse, as in the script of their Fig. 4
w$replicate <- interaction(w$exp_kernal, w$niche_breadth, w$dispersal_prob,
                           w$n.patches.pulsed, w$pre.regional.N.sum, w$pre.regional.N.sd,
                           drop = TRUE)
keys <- c("exp_kernal", "niche_breadth", "dispersal_prob", "n.patches.pulsed")

# changes from step 499 to each step after the pulse, as in data-raw/replication_hajian2026.R:
# replicates with gamma = 0 left out
# (only the steps used in the figures: 601 and every 50 steps from 650 to 1000)
ok <- w[w$post.gamma > 0 & w$post.data.ts %in% c(601, seq(650, 1000, 50)), ]
mean_se <- function(d, v, by) {
  m <- stats::aggregate(d[v], d[by], mean)
  se <- stats::aggregate(d[v], d[by], function(x) stats::sd(x) / sqrt(length(x)))
  names(se)[names(se) %in% v] <- paste0(v, "_se")
  merge(m, se)
}
changes <- mean_se(ok, c("delta.alpha", "delta.gamma", "d_beta"), c(keys, "post.data.ts"))
names(changes) <- sub("delta\\.", "d_", names(changes))

# their Fig. 4: log change in gamma during the pulse (499 to 601) and after it (601 to 1000)
g601 <- w[w$post.data.ts == 601, c(keys, "replicate", "pre.gamma", "post.gamma")]
g1000 <- w[w$post.data.ts == 1000, c("replicate", "post.gamma")]
names(g601)[names(g601) == "post.gamma"] <- "gamma_601"
names(g1000)[2] <- "gamma_1000"
ext <- merge(g601, g1000, by = "replicate")
ext$during <- log(ext$gamma_601 / ext$pre.gamma)
ext$after <- log(ext$gamma_1000 / ext$gamma_601)
ext <- ext[is.finite(ext$after), ] # as in their script
fig4 <- mean_se(ext, c("during", "after"), keys)

saveRDS(list(changes = changes, fig4 = fig4,
             n_replicates = stats::aggregate(replicate ~ exp_kernal + niche_breadth +
                                               dispersal_prob + n.patches.pulsed,
                                             unique(w[c(keys, "replicate")]), length)),
        output)
cat("saved to", output, "\n")
