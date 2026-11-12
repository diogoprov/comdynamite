suppressPackageStartupMessages(library(sf))
sf_use_s2(FALSE)
load("/Users/diogoborgesprovete/GitHub/comdynamite/.claude/worktrees/pcps-sim-package-continuation-5c5aa1/data/anurans.rda")
d <- "/private/tmp/claude-501/ecopool/downloads/ecoregions"
cg <- st_read(file.path(d, "cg.geojson"), quiet = TRUE)
sj <- st_read(file.path(d, "sjc.geojson"), quiet = TRUE)
print(st_drop_geometry(cg)); print(st_drop_geometry(sj))
eco <- rbind(cg, sj) |> st_make_valid()
eco <- eco[!duplicated(eco$eco_id), ]
pts <- st_as_sf(data.frame(anurans$coords, city = anurans$city, pond = rownames(anurans$coords)),
                coords = c("lon", "lat"), crs = 4326)
j <- st_join(pts, eco[, c("eco_name", "eco_id")])
print(table(j$city, j$eco_name, useNA = "ifany"))
keep <- eco[eco$eco_id %in% unique(j$eco_id), ]
cat("vertices:", sapply(st_geometry(keep), function(g) nrow(st_coordinates(g))), "\n")
cat("area km2 (geom):", as.numeric(st_area(st_transform(keep, "ESRI:54034")))/1e6, "\n")
simp <- st_simplify(keep, dTolerance = 0.05, preserveTopology = TRUE) |> st_make_valid()
cat("vertices simp:", sapply(st_geometry(simp), function(g) nrow(st_coordinates(g))), "\n")
simp$city <- sapply(simp$eco_id, function(i) paste(unique(as.character(j$city[j$eco_id == i])), collapse = "; "))
st_write(keep, "out/ecoregions_full.gpkg", delete_dsn = TRUE, quiet = TRUE)
st_write(simp, "out/ecoregions_simplified.gpkg", delete_dsn = TRUE, quiet = TRUE)
print(st_bbox(simp[1,])); print(st_bbox(simp[2,]))
print(st_drop_geometry(simp))
