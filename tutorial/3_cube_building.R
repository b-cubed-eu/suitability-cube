# =============================================================+===
# 3) STARS CUBE CREATION
#    - Build a cell grid over the country
#    - Aggregate AOA/DI to the grid (species × time)
#    - Build AOA/DI cubes (cell × taxon × time)
#    - Build HV cube (present only; NA future), aligned to AOA/DI
#    - Merge into a single multi-attribute stars object
# ================================================================

## 3.1 Inputs & grid params ----------------------------------------------------
species_vec  <- params$species
cellsize_deg <- params$grid_cellsize_deg
make_square  <- isTRUE(params$grid_square)

## 3.2 Country boundary (sf) ---------------------------------------------------
country_vec <- if (exists("prelim")) prelim$country_vec else geodata::gadm(params$country_name, 0, params$outdir)
country_sf  <- sf::st_as_sf(country_vec) |> sf::st_buffer(0)

## 3.3 Build the cell grid -----------------------------------------------------
grid_cells <- sf::st_make_grid(country_sf, cellsize = cellsize_deg,
                               what = "polygons", square = make_square) |>
  sf::st_as_sf() |>
  dplyr::mutate(cell = seq_len(dplyr::n()))
sf::st_crs(grid_cells) <- 4326


## 3.4 Aggregate helper (single-layer raster → grid) ---------------------------
as_stars_on_grid <- function(sr, grid, fun = mean, name = NULL, na.rm = TRUE) {
  # sr: SpatRaster (terra) or stars (single-layer)
  # grid: sf polygons (cell grid)
  
  # --- helper: mode for binary/categorical data ---
  mode_fun <- function(x, na.rm = TRUE) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  # if it's AOA (0/1), replace 'fun' with mode
  if (!is.null(name) && grepl("AOA", name, ignore.case = TRUE)) {
    fun <- mode_fun
  }
  
  # 1) convert to stars if needed
  sr_st <- if (inherits(sr, "SpatRaster")) stars::st_as_stars(sr) else sr
  if (length(names(sr_st)) != 1L) sr_st <- sr_st[1]  # keep one layer
  
  # 2) align CRS
  crs_sr   <- sf::st_crs(sr_st)
  crs_grid <- sf::st_crs(grid)
  if (!is.na(crs_sr) && !is.na(crs_grid) && crs_sr != crs_grid) {
    sr_st <- stars::st_warp(sr_st, crs = crs_grid)
  }
  
  # 3) aggregate over polygons
  agg <- suppressWarnings(stats::aggregate(sr_st, by = sf::st_geometry(grid), FUN = fun, na.rm = na.rm))
  
  # 4) rename geometry dimension to "cell"
  if ("geometry" %in% names(stars::st_dimensions(agg))) {
    agg <- stars::st_set_dimensions(agg, "geometry", names = "cell")
  }
  agg <- stars::st_set_dimensions(agg, "cell", values = sf::st_geometry(grid))
  
  # 5) optional rename
  if (!is.null(name)) names(agg) <- name
  
  agg
}

## 3.5 Build AOA/DI cubes (cell × taxon × time) --------------------------------
# keep only species that actually have AOA/DI
species_vec <- species_vec[species_vec %in% names(aoa_di_by_species) & !vapply(aoa_di_by_species, is.null, TRUE)]
stopifnot(length(species_vec) > 0)

build_metric_cube <- function(metric = c("AOA","DI")) {
  metric <- match.arg(metric)
  pres_list <- lapply(species_vec, \(sp) as_stars_on_grid(aoa_di_by_species[[sp]]$present[[metric]], grid_cells, mean, metric))
  fut_list  <- lapply(species_vec, \(sp) as_stars_on_grid(aoa_di_by_species[[sp]]$future [[metric]], grid_cells, mean, metric))
  pres <- do.call(c, pres_list) |> stars::st_redimension() |> stars::st_set_dimensions(2, values = species_vec, names = "taxon")
  fut  <- do.call(c, fut_list)  |> stars::st_redimension() |> stars::st_set_dimensions(2, values = species_vec, names = "taxon")
  cube <- c(pres, fut, along = list(time = c("present","future"))) |> stars::st_set_dimensions(1, names = "cell")
  names(cube) <- metric
  cube
}

AOA_cube <- build_metric_cube("AOA")
DI_cube  <- build_metric_cube("DI")

## 3.6 Build HV cube (present only; NA in future), aligned to AOA/DI -----------
dims  <- stars::st_dimensions(AOA_cube)   # reuse geometry + labels
shape <- dim(AOA_cube$AOA)                # c(n_cell, n_taxa, n_time)

# gather HV scalars in species order (missing species → NA)
hv_vals <- vapply(species_vec, function(sp) as.numeric(hv_by_species[[sp]]), numeric(1))
hv_arr  <- array(NA_real_, shape,
                 dimnames = list(NULL, dims$taxon$values, dims$time$values))
i_present <- match("present", dims$time$values)
for (j in seq_along(species_vec)) hv_arr[, j, i_present] <- hv_vals[j]

HV_cube <- stars::st_as_stars(list(HV = hv_arr), dimensions = dims)

## 3.7 Merge into final multi-attribute cube -----------------------------------
data_cube <- c(AOA_cube, DI_cube, HV_cube)
data_cube <- stars::st_set_dimensions(data_cube, "taxon", values = species_vec)
data_cube <- stars::st_set_dimensions(data_cube, "time",  values = c("present","future"))

# sanity check
print(stars::st_dimensions(data_cube))

## plots
library(ggplot2)
library(sf)

# ensure sf types
country_sf <- if (!inherits(country_sf, "sf")) st_as_sf(country_sf) else country_sf
grid_sf    <- if (!inherits(grid_cells, "sf")) st_as_sf(grid_cells) else grid_cells

# same CRS
if (st_crs(grid_sf) != st_crs(country_sf)) {
  grid_sf <- st_transform(grid_sf, st_crs(country_sf))
}

# quick map
ggplot() +
  geom_sf(data = country_sf, fill = "plum2", color = "red", linewidth = 0.3) +
  geom_sf(data = grid_sf, fill = NA, color = "black", linewidth = 0.2) +
  labs(title = "Analysis grid", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

# stars summary (concise)
print(data_cube)


# explicit dimensions table
stars::st_dimensions(data_cube)

# list of attributes and basic stats
names(data_cube)                # AOA, DI, HV
summary(data_cube[c("AOA","DI","HV")])

# optional: a quick sanity check plot for one species & time
# plot(data_cube["DI"][, "Bufo bufo", "present"], main = "DI — Bufo bufo (present)")
plot(data_cube)
