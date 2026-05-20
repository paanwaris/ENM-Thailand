# =============================================================================
# 00_install_packages.R
# -----------------------------------------------------------------------------
# Install every R package required by the ENM-Thailand workshop.
# Run once before starting Day 1. Re-run if you upgrade R or change machines.
# =============================================================================

# ---- CRAN dependencies ------------------------------------------------------
cran_pkgs <- c(
  # core spatial stack
  "terra", "sf", "geodata",
  # data wrangling and plotting
  "dplyr", "tidyr", "readr", "ggplot2", "viridis",
  # occurrence data + downloads
  "rgbif", "CoordinateCleaner",
  # modelling helpers (used by TemporalModelR + extras)
  "mgcv", "randomForest", "hypervolume", "fastcpd",
  # Rmd / build tools
  "rmarkdown", "knitr", "devtools", "remotes"
)

missing <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(missing) > 0) {
  message("Installing CRAN packages: ", paste(missing, collapse = ", "))
  install.packages(missing, dependencies = TRUE)
} else {
  message("All CRAN packages already installed.")
}

# ---- GitHub packages (workshop core) ----------------------------------------
gh_pkgs <- list(
  nicheR         = "castanedaM/nicheR",
  bean           = "paanwaris/bean",
  TemporalModelR = "CJHughes926/TemporalModelR"
)

for (pkg in names(gh_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, " from ", gh_pkgs[[pkg]])
    remotes::install_github(gh_pkgs[[pkg]], upgrade = "never")
  } else {
    message(pkg, " already installed.")
  }
}

message("\nAll dependencies ready. Next: source('scripts/00_download_data.R').")
