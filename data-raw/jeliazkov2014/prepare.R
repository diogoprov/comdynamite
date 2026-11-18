# Pond amphibians of Brie (Seine-et-Marne, France), 135 ponds and 11 species (Jeliazkov et al.
# 2014, Hydrobiologia 723: 7-23, doi:10.1007/s10750-013-1503-z), as distributed in the CESTES
# database (Jeliazkov et al. 2020, Scientific Data 7: 6, doi:10.1038/s41597-019-0344-7;
# Figshare doi:10.6084/m9.figshare.c.4459637, file CESTES.zip, licence CC BY 4.0 in the
# archive, CC0 on Figshare). Jeliazkov2014_AJ.xlsx is the unchanged CESTES file. This script
# writes jeliazkov2014.rds: counts (sum of five 5 x 5 m transects per pond), coordinates in km
# (Lambert 93) and the environmental gradient used in the application, the first axis of a PCA
# of the eight quantitative environmental variables (standardised; pond area on the log scale),
# rescaled to 0-100. Note: in the sheet splist the families of Alytes obstetricans and Bufo bufo
# are swapped; the species names are right.
f <- "data-raw/jeliazkov2014/Jeliazkov2014_AJ.xlsx"
comm <- as.data.frame(readxl::read_excel(f, "comm"))
envir <- as.data.frame(readxl::read_excel(f, "envir"))
coord <- as.data.frame(readxl::read_excel(f, "coord"))
splist <- as.data.frame(readxl::read_excel(f, "splist"))[1:11, 1:3]
stopifnot(identical(comm$Sites, envir$Sites), identical(comm$Sites, coord$Sites))
counts <- as.matrix(comm[, -1])
rownames(counts) <- paste0("pond", comm$Sites)
colnames(counts) <- splist$Taxon[match(colnames(counts), splist$TaxCode)]
vars <- data.frame(Fish = envir$Fish, WQI = envir$WQI, PropWood1000 = envir$PropWood1000,
                   PondDens_1000 = envir$PondDens_1000, logArea = log10(envir$Area),
                   AgricGrad = envir$AgricGrad, UrbanGrad = envir$UrbanGrad,
                   SemiNatGrad = envir$SemiNatGrad)
pca <- stats::prcomp(vars, scale. = TRUE)
pc1 <- pca$x[, 1]
gradient <- 100 * (pc1 - min(pc1)) / diff(range(pc1))
coords <- cbind(x = coord$X, y = coord$Y) / 1000
rownames(coords) <- rownames(counts)
out <- list(counts = counts, coords = coords, env = stats::setNames(gradient, rownames(counts)),
            variables = vars, loadings = pca$rotation[, 1],
            variance_pc1 = summary(pca)$importance[2, 1])
saveRDS(out, "data-raw/jeliazkov2014/jeliazkov2014.rds")
print(round(out$loadings, 2)); print(out$variance_pc1)
print(round(apply(coords, 2, range), 1))
