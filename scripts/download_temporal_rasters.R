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
# OPTION A — Zenodo (recommended for academic / DOI-citable hosting).
#   After uploading temporal_rasters_2010-2025.zip to Zenodo, paste the
#   direct-download URL here. It looks like:
#     https://zenodo.org/records/XXXXXXX/files/temporal_rasters_2010-2025.zip
zenodo_url <- ""

# OPTION B — GitHub release attached to this repo.
#   Tag a release in the GitHub UI (e.g. "data-v1"), upload the zip as an
#   asset. The URL pattern is:
#     https://github.com/paanwaris/ENM-Thailand/releases/download/<tag>/<file>
gh_release_url <- ""

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
target_dir <- "temporal_rasters"
dir.create(target_dir, showWarnings = FALSE, recursive = TRUE)

zip_path <- tempfile(fileext = ".zip")
message("Downloading from: ", url)
utils::download.file(url, destfile = zip_path, mode = "wb")

message("Unpacking into: ", normalizePath(target_dir))
utils::unzip(zip_path, exdir = target_dir, junkpaths = TRUE)
unlink(zip_path)

# ---- 3. Verify -------------------------------------------------------------
years    <- 2010:2025
expected <- file.path(target_dir,
                      sprintf("Thailand_LST_Precip_1km_%d.tif", years))
missing  <- expected[!file.exists(expected)]

if (length(missing)) {
  warning(
    "The following expected files are still missing after download:\n  ",
    paste(missing, collapse = "\n  "),
    "\nCheck the zip contents on the host (filenames, folder structure)."
  )
} else {
  message("All ", length(expected), " annual rasters present.")
}
