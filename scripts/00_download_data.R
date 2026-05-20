# =============================================================================
# 00_download_data.R
# -----------------------------------------------------------------------------
# Day-0 data preparation for the ENM-Thailand workshop.
#
# Outputs (saved under data/processed/):
#   * thailand_boundary.gpkg    — Thailand national polygon (GADM v4.1)
#   * bioclim_thailand.tif      — WorldClim v2.1 bioclim, 10 arc-min, cropped
#   * sambar_gbif_raw.csv       — raw GBIF download for Rusa unicolor
#   * sambar_clean_thailand.csv — cleaned, Thailand-only occurrences
#
# Coordinate reference system: WGS84 / EPSG:4326
# Bioclim resolution:          10 arc-minutes (~18 km at the equator)
#
# Run once. Day 2 and Day 3 Rmds load directly from data/processed/.
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(geodata)
  library(rgbif)
  library(CoordinateCleaner)
  library(dplyr)
  library(readr)
})

# ---- Project paths ----------------------------------------------------------
raw_dir       <- "data/raw"
processed_dir <- "data/processed"
dir.create(raw_dir,       recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Thailand boundary (GADM) --------------------------------------------
message("[1/3] Downloading Thailand boundary from GADM ...")
tha <- geodata::gadm(country = "THA", level = 0, path = raw_dir)
tha_sf <- sf::st_as_sf(tha)
sf::st_write(tha_sf,
             file.path(processed_dir, "thailand_boundary.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ---- 2. WorldClim bioclim (10 arc-min, global → cropped) -------------------
message("[2/3] Downloading WorldClim v2.1 bioclim @ 10 arc-min ...")
bio_global <- geodata::worldclim_global(var = "bio", res = 10, path = raw_dir)

message("       Cropping + masking to Thailand ...")
bio_tha <- terra::crop(bio_global, terra::ext(tha))
bio_tha <- terra::mask(bio_tha, tha)

# Tidy names: wc2.1_10m_bio_1 -> bio1
names(bio_tha) <- paste0("bio", seq_len(terra::nlyr(bio_tha)))

terra::writeRaster(bio_tha,
                   filename  = file.path(processed_dir, "bioclim_thailand.tif"),
                   overwrite = TRUE)

# ---- 3. GBIF occurrences for Sambar deer (Rusa unicolor) -------------------
message("[3/3] Downloading GBIF occurrences for Rusa unicolor ...")

# Resolve the taxon key once, then pull records limited to Thailand.
key <- rgbif::name_backbone(name = "Rusa unicolor")$usageKey

raw <- rgbif::occ_search(
  taxonKey            = key,
  country             = "TH",
  hasCoordinate       = TRUE,
  hasGeospatialIssue  = FALSE,
  limit               = 5000
)$data

readr::write_csv(raw, file.path(processed_dir, "sambar_gbif_raw.csv"))
message("       Raw rows: ", nrow(raw))

# ---- Clean coordinates ------------------------------------------------------
sambar <- raw %>%
  dplyr::select(species,
                decimalLongitude, decimalLatitude,
                year, month,
                basisOfRecord, gbifID) %>%
  dplyr::filter(!is.na(decimalLongitude),
                !is.na(decimalLatitude),
                !is.na(year)) %>%
  dplyr::distinct(decimalLongitude, decimalLatitude, year, .keep_all = TRUE)

# CoordinateCleaner flags zero coords, centroids, capitals, GBIF HQ, etc.
flags <- CoordinateCleaner::clean_coordinates(
  x       = sambar,
  lon     = "decimalLongitude",
  lat     = "decimalLatitude",
  species = "species",
  tests   = c("zeros", "capitals", "centroids", "equal", "gbif", "institutions"),
  verbose = FALSE
)

sambar_clean <- sambar[flags$.summary, , drop = FALSE]

# Keep only points that actually fall inside Thailand (just in case)
sambar_sf <- sf::st_as_sf(sambar_clean,
                          coords = c("decimalLongitude", "decimalLatitude"),
                          crs    = 4326, remove = FALSE)
inside       <- as.logical(sf::st_intersects(sambar_sf, tha_sf, sparse = FALSE))
sambar_clean <- sambar_clean[inside, , drop = FALSE]

readr::write_csv(sambar_clean,
                 file.path(processed_dir, "sambar_clean_thailand.csv"))
message("       Clean rows in Thailand: ", nrow(sambar_clean))

message("\nDone. Files written to ", processed_dir)
