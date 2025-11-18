############################################################
# ANALYSIS – Quick recipes over the stars data cube
# dims(data_cube): cell × taxon × time ; attrs: AOA, DI, HV
############################################################

# ============================================================
# 1) LOCATE A CELL FROM A LON/LAT AND PLOT IT (with raster bg)
# ============================================================
# Define a point (lon, lat) in EPSG:4326
pt <- st_sf(geometry = st_sfc(st_point(c(12.5, 42.5)), crs = 4326))

# Which grid cell contains the point?
which_cell <- suppressWarnings(st_join(pt, grid_cells, join = st_intersects))

if (is.na(which_cell$cell)) {
  message("❌ This point is NOT within the study area.")
} else {
  cell_id <- which_cell$cell
  message(sprintf("✅ Point falls inside cell #%d", cell_id))
  
  # Selected cell polygon
  selected_cell <- grid_cells[grid_cells$cell == cell_id, ]
  
  # Use any raster as background (here: first present bioclim layer)
  base_raster <- st_as_stars(bio_present_sel[[1]])
  
  # Optional Italy background if available
  italy_layer <- if (exists("italy")) st_as_sf(italy) else NULL
  
  ggplot() +
    geom_stars(data = base_raster) +
    scale_fill_viridis_c(option = "C") +
    { if (!is.null(italy_layer)) geom_sf(data = italy_layer, fill = NA, color = "black", linewidth = 0.3) } +
    geom_sf(data = grid_cells, fill = NA, color = "grey60", linewidth = 0.1) +
    geom_sf(data = selected_cell, fill = "orange", color = "red", linewidth = 0.4) +
    geom_sf(data = pt, color = "blue", size = 2) +
    labs(
      title = paste("Cell", cell_id, "containing the selected point"),
      subtitle = "Orange = selected cell | Blue = queried point"
    ) +
    theme_minimal()
}

# ============================================================
# 2) BASIC INTROSPECTION & SLICING
# ============================================================
# Spatial extent of the cube (bbox of all cells)
st_bbox(data_cube)

# Slice example (confirm dims order): [cell, taxon, time]
# AOA & DI arrays’ shape
dim(data_cube[c("AOA","DI")])

# Cell 1361, both species, PRESENT (time = 1)
data_cube[,1361, , 1]

# Cell 1361, both species, FUTURE (time = 2)
data_cube[,1361, , 2]

# Pull just DI array (all species, both times) for cell 1361
pull(data_cube["DI",1361, , ])

# ============================================================
# 3) CROP THE CUBE AROUND A SMALL WINDOW (centroid + buffer)
# ============================================================
b <- st_as_sfc(st_bbox(data_cube)) |>
  st_centroid() |>
  st_buffer(units::set_units(500, m))  # ~500 m buffer (adjust as needed)

plot(b)            # polygon used for cropping
data_cube[b]       # cropped cube

# ============================================================
# 4) PERMUTE DIMS OR SPLIT BY ATTRIBUTE (advanced)
# ============================================================
# Permute: e.g., put time first, then cell, then taxon
aperm(data_cube, c(3, 1, 2))

# Split into single-attribute cubes (list of stars objects)
rs <- split(data_cube)   # rs$AOA, rs$DI, rs$HV

# ============================================================
# 5) BUILD A PAIRWISE DI-DIFFERENCE CUBE (cell × comparison × time)
#    DI_diff = (species_i − species_j) for each pair
# ============================================================

build_DI_diff_cube <- function(data_cube, pairs = NULL) {
  dims <- st_dimensions(data_cube)
  tax  <- trimws(as.character(dims$taxon$values))
  di   <- data_cube["DI"]
  arr  <- di[[1]]  # array [cell, taxon, time]
  
  stopifnot(length(tax) >= 2)
  ncell <- dim(arr)[1]; ntax <- dim(arr)[2]; ntime <- dim(arr)[3]
  
  # Resolve pairs (default: all i<j)
  if (is.null(pairs)) {
    pairs_idx <- t(combn(ntax, 2))
  } else if (is.character(pairs)) {
    pairs_idx <- cbind(match(trimws(pairs[,1]), tax),
                       match(trimws(pairs[,2]), tax))
  } else {
    pairs_idx <- as.matrix(pairs)  # numeric indices
  }
  if (any(is.na(pairs_idx))) stop("Pairs include unknown species labels.")
  if (ncol(pairs_idx) != 2)  stop("`pairs` must have two columns (sp1, sp2).")
  
  K <- nrow(pairs_idx)
  arr_diff <- array(NA_real_, dim = c(ncell, K, ntime))
  comp_labels <- character(K)
  
  for (k in seq_len(K)) {
    i <- pairs_idx[k, 1]; j <- pairs_idx[k, 2]
    arr_diff[, k, ] <- arr[, i, ] - arr[, j, ]
    comp_labels[k]  <- paste0(tax[i], " - ", tax[j])
  }
  
  # Reuse dims but rename 2nd dim to 'comparison'
  dims_new <- dims
  names(dims_new)[2] <- "comparison"
  dims_new$comparison$values <- comp_labels
  
  st_as_stars(list(DI_diff = arr_diff), dimensions = dims_new)
}

# Example: all pairwise differences (default)
DI_diff_cube <- build_DI_diff_cube(data_cube)
DI_diff_cube

# Example: specific order/direction (indices are 1-based)
# pairs_idx <- rbind(c(1,2), c(2,1), c(2,3))
# DI_diff_cube <- build_DI_diff_cube(data_cube, pairs = pairs_idx)

# ============================================================
# 6) PLOT DI DIFFERENCES FOR ONE CELL (bars by comparison & time)
# ============================================================
cell_id <- 1361
slice_diff <- DI_diff_cube[,cell_id , , drop = FALSE]
df_diff <- as.data.frame(slice_diff) |> select(comparison, time, DI_diff)

ggplot(df_diff, aes(x = comparison, y = DI_diff, fill = time)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(
    title = paste("DI difference at cell", cell_id),
    subtitle = "Positive = higher DI for the first species in the pair",
    x = "Comparison (sp1 − sp2)", y = "DI difference"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# ============================================================
# 7) PLOT DI PER SPECIES FOR ONE CELL (present vs future)
# ============================================================
cell_id <- 1361
df_DI <- as.data.frame(data_cube["DI",cell_id , , drop = FALSE])

ggplot(df_DI, aes(x = taxon, y = DI, fill = time)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  labs(title = paste("DI at cell", cell_id), x = "Species", y = "DI") +
  theme_minimal()

# ============================================================
# 8) COUNT AOA CELLS (0/1) BY SPECIES AND TIME (whole area)
# ============================================================
aoa_df <- as.data.frame(data_cube["AOA"]) |>
  select(taxon, time, AOA) |>
  mutate(AOA = as.integer(round(AOA)))  # ensure 0/1

aoa_counts <- aoa_df |>
  group_by(taxon, time, AOA) |>
  summarise(n_cells = n(), .groups = "drop")

print(aoa_counts)

ggplot(aoa_counts, aes(x = taxon, y = n_cells, fill = factor(AOA))) +
  geom_col(width = 0.7) +
  facet_wrap(~ time) +
  scale_fill_manual(values = c("#bbbbbb", "#2c7fb8"), name = "AOA",
                    labels = c("0 (outside AOA)", "1 (inside AOA)")) +
  labs(title = "AOA cell counts by species", x = "Species", y = "Number of cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# ============================================================
# 9) HYPERVOLUME BARS (present only; constant across space)
# ============================================================
hv_df <- as.data.frame(data_cube["HV", , ,1, drop = FALSE]) |>
  select(taxon, HV) |>
  group_by(taxon) |>
  summarise(HV = suppressWarnings(first(na.omit(HV))), .groups = "drop")

print(hv_df)

ggplot(hv_df, aes(x = taxon, y = HV)) +
  geom_col(width = 0.7) +
  labs(title = "Environmental hypervolume (present)", x = "Species", y = "Hypervolume") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
