# =============================================================================
# download_temporal_rasters.R
# -----------------------------------------------------------------------------
# Fetch the annual LST + Precip + NDVI rasters (2010–2025) that Day 3 needs.
#
# The rasters are too large for git, so they live in an external archive.
# Choose ONE host and uncomment the matching URL. Then run:
#     source("scripts/download_temporal_rasters.R")
#
# The script downloads a zip, unpacks it into `temporal_rasters/`, and
# verifies the expected files are present.
# =============================================================================

# ---- 1. Source URL ---------------------------------------------------------

zenodo_url <- "https://zenodo.org/records/20344510/files/processed.zip?download=1"

# Pick whichever you filled in
url <- if (nzchar(zenodo_url)) zenodo_url else gh_release_url
if (!nzchar(url)) {
  stop(
    "No download URL configured.\n",
    "Edit scripts/download_temporal_rasters.R and set either ",
    "`zenodo_url` or `gh_release_url`."
  )
}

# ---- 2. Download + unpack --------------------------------------------------
target_dir <- "data/"
dir.create(target_dir, showWarnings = FALSE, recursive = TRUE)

zip_path <- tempfile(fileext = ".zip")
message("Downloading from: ", url)
utils::download.file(url, destfile = zip_path, mode = "wb")

message("Unpacking into: ", normalizePath(target_dir))
utils::unzip(zip_path, exdir = target_dir, overwrite = TRUE)
unlink(zip_path)
