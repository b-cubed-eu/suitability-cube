# ================================================================
# 2) CONTENT
#    - Hypervolume on PRESENT (per species)
#    - AOA & DI on PRESENT and FUTURE (per species)
# ================================================================

stopifnot(exists("prelim"))
attach(prelim)  # brings: params, bio_present_sel, bio_future_sel, occ_list, ...

## 2.1 Helpers ---------------------------------------------------------------

# Extract predictors at occurrences (keep only selected vars; drop NA rows; drop zero-variance cols)
extract_predictors_at_points <- function(rst, occ_sf, pred_vars) {
  if (!is.na(sf::st_crs(occ_sf)) && !is.na(terra::crs(rst)) &&
      sf::st_crs(occ_sf)$wkt != terra::crs(rst)) {
    occ_sf <- sf::st_transform(occ_sf, crs = terra::crs(rst))
  }
  vals <- terra::extract(rst, occ_sf)
  df   <- dplyr::select(as.data.frame(vals), dplyr::all_of(pred_vars)) |> stats::na.omit()
  keep <- vapply(df, function(x) is.numeric(x) && sd(x, na.rm = TRUE) > 0, logical(1))
  df[, keep, drop = FALSE]
}

# Z-standardization with mean imputation for NAs
z_transform <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      nas <- is.na(col); m <- mean(col, na.rm = TRUE); col[nas] <- m
      s <- sd(col, na.rm = TRUE); if (is.na(s) || s == 0) return(rep(0, length(col)))
      (col - m) / s
    } else col
  }); df
}

# Hypervolume helpers
compute_global_bw <- function(env_df) {
  bw <- try(hypervolume::estimate_bandwidth(as.matrix(env_df)), silent = TRUE)
  if (inherits(bw, "try-error") || any(!is.finite(bw))) 1 else bw
}
hyp_calc <- function(env_df, bw, spp = 50L) {
  X <- as.matrix(env_df)
  if (nrow(unique(X)) < (ncol(X) + 1L)) return(NA_real_)
  hv <- try(hypervolume::hypervolume_gaussian(X, kde.bandwidth=bw, samples.per.point=spp, verbose=FALSE), silent=TRUE)
  if (inherits(hv, "try-error")) return(NA_real_)
  hv@Volume
}

# Compute AOA/DI for present & future given a training DF
compute_aoa_pair <- function(env_train_df, new_present, new_future) {
  keep <- vapply(env_train_df, function(x) is.numeric(x) && sd(x, na.rm = TRUE) > 0, logical(1))
  trn  <- env_train_df[, keep, drop = FALSE]; vars <- colnames(trn)
  list(
    present = CAST::aoa(newdata = new_present, train = trn, variables = vars),
    future  = CAST::aoa(newdata = new_future,  train = trn, variables = vars)
  )
}

## 2.2 Hypervolume (PRESENT) ---------------------------------------------------
pred_vars_present <- names(bio_present_sel)
hv_by_species     <- setNames(vector("list", length(params$species)), params$species)

for (sp in params$species) {
  occ_sf <- occ_list[[sp]]
  if (is.null(occ_sf) || nrow(occ_sf) == 0) { hv_by_species[[sp]] <- NA_real_; next }
  train_df  <- extract_predictors_at_points(bio_present_sel, occ_sf, pred_vars_present)
  if (nrow(train_df) < (ncol(train_df)+1L)) { hv_by_species[[sp]] <- NA_real_; next }
  z_df <- z_transform(train_df)
  hv_by_species[[sp]] <- hyp_calc(z_df, compute_global_bw(z_df))
}

hv_by_species
## 2.3 AOA & DI (PRESENT & FUTURE) --------------------------------------------
aoa_di_by_species <- setNames(vector("list", length(params$species)), params$species)

for (sp in params$species) {
  occ_sf <- occ_list[[sp]]
  if (is.null(occ_sf) || nrow(occ_sf) == 0) { aoa_di_by_species[[sp]] <- NULL; next }
  train_df <- extract_predictors_at_points(bio_present_sel, occ_sf, pred_vars_present)
  if (nrow(train_df) < 5 || ncol(train_df) < 1) { aoa_di_by_species[[sp]] <- NULL; next }
  aoa_di_by_species[[sp]] <- compute_aoa_pair(train_df, new_present = bio_present_sel, new_future = bio_future_sel)
}

detach(prelim)

plot(aoa_di_by_species$`Bufo bufo`$present$DI)
plot(aoa_di_by_species$`Bufo bufo`$present$AOA)
                 
