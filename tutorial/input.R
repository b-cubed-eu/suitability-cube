# ================================================================
# 1) PRELIMINARIES
#    - User inputs
#    - Helpers
#    - Country boundary & climate rasters (present/future)
#    - GBIF occurrences (multi-species)
#    - Correlation-based variable selection (on PRESENT)
# ================================================================

## 1.1 Packages ---------------------------------------------------------------
library(geodata)
library(terra)
library(sf)
library(tidyverse)
library(rgbif)
library(corrplot)
library(viridis)
library(stars)
library(CAST)
library(hypervolume)

## 1.2 User inputs (EDIT here) ------------------------------------------------
params <- list(
  species      = c("Bufo bufo", "Bufotes viridis", "Pelophylax esculentus"),
  country_name = "Italy",
  country_iso  = "IT",
  res_arcmin   = 2.5,
  ssp_code     = "245",
  gcm_model    = "BCC-CSM2-MR",
  period       = "2041-2060",
  outdir       = tempdir(),      # change to a persistent path if desired
  gbif_years   = c(2010, 2020),
  gbif_limit   = 20000,
  cor_thr      = 0.7,
  cor_frac     = 0.10,
  grid_cellsize_deg = 0.25,      # ~25 km
  grid_square      = FALSE       # FALSE = hex, TRUE = square
)

## 1.3 Helpers ---------------------------------------------------------------

# Align a raster stack to a target grid (bilinear)
align_to <- function(src, target) terra::resample(src, target, method = "bilinear")

# Drop highly correlated layers using a simple greedy scheme
drop_high_corr <- function(rst, thr = 0.7, frac = 0.1, seed = 42) {
  set.seed(seed)
  sz  <- max(1000, round(terra::ncell(rst) * frac))
  smp <- terra::spatSample(rst, size = sz, method = "random", na.rm = TRUE, as.points = FALSE)
  mat <- as.matrix(smp)
  keep_cols <- which(colSums(!is.na(mat)) > 0)
  mat <- mat[, keep_cols, drop = FALSE]
  cm  <- suppressWarnings(cor(mat, use = "pairwise.complete.obs", method = "pearson"))
  to_drop <- character(0); vars <- colnames(cm)
  avg_abs <- sort(colMeans(abs(cm), na.rm = TRUE), decreasing = TRUE)
  for (v in names(avg_abs)) {
    if (v %in% to_drop) next
    high <- setdiff(names(which(abs(cm[v, ]) > thr)), v)
    to_drop <- union(to_drop, high)
  }
  list(cor = cm, selected = setdiff(vars, to_drop), dropped = to_drop)
}

# GBIF search → sf for a single species
gbif_occ_sf <- function(scientific_name, country_iso = "IT", years = c(2000, 2020), limit = 20000) {
  key <- rgbif::name_backbone(name = scientific_name)$usageKey
  dat <- rgbif::occ_search(
    taxonKey = key,
    country  = country_iso,
    year     = paste0(years[1], ",", years[2]),
    basisOfRecord = "HUMAN_OBSERVATION",
    occurrenceStatus = "PRESENT",
    hasCoordinate   = TRUE,
    limit = limit
  )$data
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  dat <- dat[!is.na(dat$decimalLongitude) & !is.na(dat$decimalLatitude), ]
  dat <- dat[!duplicated(dat[, c("decimalLongitude","decimalLatitude")]), ]
  sf::st_as_sf(dat, coords = c("decimalLongitude","decimalLatitude"), crs = 4326, remove = FALSE)
}

# Multi-species wrapper (named list)
gbif_occ_list <- function(species_vec, country_iso, years, limit) {
  setNames(lapply(species_vec, \(sp) gbif_occ_sf(sp, country_iso, years, limit)), species_vec)
}

## 1.4 Country boundary & climate rasters -------------------------------------
## you can replace it with your own vector
country_vec <- geodata::gadm(params$country_name, level = 0, path = params$outdir)

bio_present <- geodata::worldclim_country(
  country = params$country_name, var = "bio", res = params$res_arcmin, path = params$outdir
) |> terra::crop(country_vec) |> terra::mask(country_vec)

bio_future <- geodata::cmip6_world(
  model = params$gcm_model, ssp = params$ssp_code, time = params$period,
  var = "bio", res = params$res_arcmin, path = params$outdir
) |> terra::crop(country_vec) |> terra::mask(country_vec)

bio_present_aligned <- align_to(bio_present, bio_future)
names(bio_future) <- names(bio_present_aligned)

message("Present (aligned) res: ", paste(terra::res(bio_present_aligned), collapse = ", "))
message("Future res:            ", paste(terra::res(bio_future),         collapse = ", "))
print(terra::compareGeom(bio_present_aligned, bio_future, stopOnError = FALSE))

## 1.5 GBIF occurrences --------------------------------------------------------
occ_list <- gbif_occ_list(params$species, params$country_iso, params$gbif_years, params$gbif_limit)

## 1.6 Variable selection (on PRESENT) ----------------------------------------
cor_res   <- drop_high_corr(bio_present_aligned, thr = params$cor_thr, frac = params$cor_frac)
cmat      <- cor_res$cor
vars_keep <- cor_res$selected
vars_drop <- cor_res$dropped

bio_present_sel <- bio_present_aligned[[vars_keep]]
bio_future_sel  <- bio_future[[vars_keep]]

# (optional) visualize correlation matrix
corrplot::corrplot(cmat, method="number", col=viridis::magma(30), type="lower", tl.pos="ld", number.cex=0.6)

# Hand-off bundle for next sections
prelim <- list(
  params=params, country_vec=country_vec,
  bio_present_sel=bio_present_sel, bio_future_sel=bio_future_sel,
  occ_list=occ_list, vars_keep=vars_keep, vars_drop=vars_drop
)
