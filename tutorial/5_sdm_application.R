#############################################
# 1) SPECIES SDM FIT + PREDICT (dismo::bioclim)
#    - For each species:
#        a) get GBIF occurrences
#        b) train a simple BIOCLIM model on present climate
#        c) predict suitability for present and future climate
#############################################

## ---- BUFO BUFO -------------------------------------------------------

# 1.1 Extract occurrence points (sf) for Bufo bufo
occ_bufo <- occ_list[["Bufo bufo"]]

# 1.2 Force to WGS84 lon/lat coordinates (EPSG:4326)
occ_bufo <- st_transform(occ_bufo, crs = "EPSG:4326")

# 1.3 Extract lon/lat matrix for dismo
pres_xy_bufo <- st_coordinates(occ_bufo)

# 1.4 Fit a simple BIOCLIM SDM using current climate predictors
bc_bufo <- dismo::bioclim(env_present_rs, pres_xy_bufo)

# 1.5 Predict habitat suitability under present climate
suit_present_bufo <- raster::predict(env_present_rs, bc_bufo)

# 1.6 Predict habitat suitability under future climate
suit_future_bufo  <- raster::predict(env_future_rs,  bc_bufo)
dev.off()
# (optional quick plot)
par(mfrow = c(1,2))
plot(suit_present_bufo, main = "Bufo bufo - suitability (present)")
points(pres_xy_bufo, pch = 20, cex = 0.4)
plot(suit_future_bufo, main = "Bufo bufo - suitability (future)")
par(mfrow = c(1,1))


## ---- BUFOTES VIRIDIS -------------------------------------------------

occ_viridis <- occ_list[["Bufotes viridis"]]
occ_viridis <- st_transform(occ_viridis, crs = "EPSG:4326")
pres_xy_viridis <- st_coordinates(occ_viridis)

bc_viridis <- dismo::bioclim(env_present_rs, pres_xy_viridis)

suit_present_viridis <- raster::predict(env_present_rs, bc_viridis)

suit_future_viridis  <- raster::predict(env_future_rs,  bc_viridis)

par(mfrow = c(1,2))
plot(suit_present_viridis, main = "Bufotes viridis - suitability (present)")
points(pres_xy_viridis, pch = 20, cex = 0.4)
plot(suit_future_viridis, main = "Bufotes viridis - suitability (future)")
par(mfrow = c(1,1))


## ---- BOMBINA VARIEGATA ------------------------------------------

occ_bombina <- occ_list[["Bombina variegata"]]
occ_bombina <- st_transform(occ_bombina, crs = "EPSG:4326")
pres_xy_bombina <- st_coordinates(occ_bombina)

bc_bombina <- dismo::bioclim(env_present_rs, pres_xy_bombina)

suit_present_bombina <- raster::predict(env_present_rs, bc_bombina)
suit_future_bombina  <- raster::predict(env_future_rs,  bc_bombina)

par(mfrow = c(1,2))
plot(suit_present_bombina, main = "Bombina variegata - suitability (present)")
points(pres_xy_bombina, pch = 20, cex = 0.4)
plot(suit_future_bombina, main = "Bombina variegata - suitability (future)")
par(mfrow = c(1,1))

dev.off()
#############################################
# 2) BUILD SUITABILITY CUBE
#    Goal: create a unified stars data cube with dimensions:
#      - cell (spatial grid cell polygons)
#      - taxon (species)
#      - time ("present", "future")
#    and a single attribute: "suitability"
#
#    Strategy:
#      - aggregate each suitability raster to the analysis grid (grid_cells)
#      - stack species -> taxon dimension
#      - stack present/future -> time dimension
#############################################

# 2.1 Species order (must match everything else in the project)
species_vec <- params$species
stopifnot(all(species_vec == c("Bufo bufo", "Bufotes viridis", "Bombina variegata")))

# 2.2 Put all suitability rasters into named lists (one for present, one for future)
#     NOTE: these are RasterLayer objects right now
suitability_present_list <- list(
  "Bufo bufo"              = suit_present_bufo,
  "Bufotes viridis"        = suit_present_viridis,
  "Bombina variegata"  = suit_present_bombina
)

suitability_future_list <- list(
  "Bufo bufo"              = suit_future_bufo,
  "Bufotes viridis"        = suit_future_viridis,
  "Bombina variegata"  = suit_future_bombina
)

suitability_present_list
species_vec
# 2.3 Aggregate suitability to the analysis grid for each species and time.
#     The helper 'as_stars_on_grid()' already:
#       - converts to stars / aligns CRS
#       - aggregates values over each polygon in grid_cells
#       - returns a stars object with a "cell" dimension
#
#     We convert RasterLayer -> SpatRaster on the fly via terra::rast()

suit_present_grid_list <- lapply(species_vec, function(sp) {
  as_stars_on_grid(
    sr   = terra::rast(suitability_present_list[[sp]]),  # RasterLayer -> SpatRaster
    grid = grid_cells,                                   # sf grid of polygons (cell IDs)
    fun  = mean,                                         # average suitability per grid cell
    name = "suitability",
    na.rm = TRUE
  )
})

suit_future_grid_list <- lapply(species_vec, function(sp) {
  as_stars_on_grid(
    sr   = terra::rast(suitability_future_list[[sp]]),
    grid = grid_cells,
    fun  = mean,
    name = "suitability",
    na.rm = TRUE
  )
})

# 2.4 Stack species along a new "taxon" dimension for PRESENT
suit_present_cube <- do.call(c, suit_present_grid_list) |>
  stars::st_redimension() |>
  stars::st_set_dimensions(2, values = species_vec, names = "taxon")

# 2.5 Stack species along "taxon" for FUTURE
suit_future_cube <- do.call(c, suit_future_grid_list) |>
  stars::st_redimension() |>
  stars::st_set_dimensions(2, values = species_vec, names = "taxon")

# 2.6 Stack PRESENT and FUTURE along a new "time" dimension
suit_cube <- c(
  suit_present_cube,
  suit_future_cube,
  along = list(time = c("present", "future"))
)

# 2.7 Make sure dimension names and attribute name are clean/standard
suit_cube <- stars::st_set_dimensions(suit_cube, 1, names = "cell")
names(suit_cube) <- "suitability"

# Quick sanity check
suit_cube
stars::st_dimensions(suit_cube)

# Optional: visualize (by default stars will facet by taxon/time)
plot(suit_cube)
suit_cube

#############################################
# 3) MERGE SUITABILITY INTO THE GLOBAL data_cube
#    (Skip this if you just want the standalone suitability cube.
#     Keep this if you already have 'data_cube' with AOA / DI / HV
#     and you want to add suitability as a new attribute.)
#############################################

# 3.1 Align taxon/time labels to match the existing data_cube
suit_cube <- stars::st_set_dimensions(
  suit_cube,
  "taxon",
  values = stars::st_dimensions(data_cube)$taxon$values
)

suit_cube <- stars::st_set_dimensions(
  suit_cube,
  "time",
  values = stars::st_dimensions(data_cube)$time$values
)

# 3.2 Add "suitability" as a new attribute/band in data_cube
data_cube <- c(data_cube, suit_cube)

# Check that dimensions are still what we expect
print(stars::st_dimensions(data_cube))
names(data_cube)  # should now include "suitability"

data_cube

AOA_cube
#############################################
# 4) CLIP SUITABILITY BY AOA
#    Goal: produce an AoA-masked suitability, i.e.
#    suitability is kept only where the model is considered valid.
#
#    Inputs:
#      - suit_cube$suitability : [cell, taxon, time]
#      - AOA_cube$AOA          : [cell, taxon, time] binary (1 = inside AoA, 0 = outside)
#
#    Output:
#      - suit_cube_masked with attribute "suitability_masked"
#        where cells outside AoA are set to NA
#############################################

# 4.1 First, be sure suit_cube and AOA_cube share the same dimension labels
suit_cube <- stars::st_set_dimensions(
  suit_cube,
  "taxon",
  values = stars::st_dimensions(AOA_cube)$taxon$values
)

suit_cube <- stars::st_set_dimensions(
  suit_cube,
  "time",
  values = stars::st_dimensions(AOA_cube)$time$values
)

# Optional checks
stars::st_dimensions(suit_cube)
stars::st_dimensions(AOA_cube)

# 4.2 Extract raw arrays
suit_arr <- suit_cube$suitability   # numeric array [cell, taxon, time]
aoa_arr  <- AOA_cube$AOA            # 0/1 or TRUE/FALSE array [cell, taxon, time]

# 4.3 Mask: set suitability to NA wherever AoA == 0 (or AoA is NA)
suit_arr_masked <- suit_arr
suit_arr_masked[ aoa_arr == 0 | is.na(aoa_arr) ] <- NA

# 4.4 Rebuild a stars object with the masked suitability
suit_cube_masked <- stars::st_as_stars(
  list(suitability_masked = suit_arr_masked),
  dimensions = stars::st_dimensions(suit_cube)
)

# Visual check: this should show suitability only where AoA says "valid"
plot(suit_cube_masked)
suit_cube_masked
