# =============================================================================
# build_pkgdown_site.R
# -----------------------------------------------------------------------------
# Build the public pkgdown website for the ENM-Thailand workshop.
#
# The repository is "package-shaped" (it has a DESCRIPTION so pkgdown gets
# clean metadata + ORCID rendering), but there is no `R/` folder with code.
# So we call the individual pkgdown build steps explicitly and skip
# `build_reference()` (which would try to `library(ENMThailand)`).
#
# Prerequisites:
#   * The shared inputs (data/processed/, temporal_rasters/) already exist —
#     knit Day2a_Data_Download.Rmd first so every vignette can run end-to-end.
#   * `pkgdown` is installed.
# =============================================================================

needed  <- c("pkgdown", "knitr", "rmarkdown")
missing <- needed[!needed %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing, dependencies = TRUE)

notebook_rmds <- file.path(
  "notebooks",
  c("Day1_nicheR.Rmd",
    "Day2a_Data_Download.Rmd",
    "Day2b_Bean_Processing.Rmd",
    "Day3_TemporalModelR.Rmd")
)
missing_rmds <- notebook_rmds[!file.exists(notebook_rmds)]
if (length(missing_rmds)) {
  stop("Missing canonical Rmds in notebooks/: ",
       paste(missing_rmds, collapse = ", "))
}

# Tell Day 2a to render its code without hitting GADM / WorldClim / GBIF.
# Day 1, 2b, 3 use their own file-existence guards.
Sys.setenv(PKGDOWN_BUILD = "1")
on.exit(Sys.unsetenv("PKGDOWN_BUILD"), add = TRUE)

# ---- Build the site, step by step -------------------------------------------
# init_site()  copies assets (CSS, JS, favicons) into docs/.
# build_home() renders README.md + LICENSE.md as index.html / LICENSE.html.
# build_articles() renders every vignette under vignettes/ (the four workshop
#                  days, via the `child = "../DayX_*.Rmd"` wrappers).
# NOTE: We deliberately skip build_reference() because the repository has no
#       R/ functions to document.
# build_news()   is a no-op unless NEWS.md exists.
# build_search() builds the in-site search index.

pkgdown::init_site()
pkgdown::build_home(preview = FALSE)
pkgdown::build_articles(preview = FALSE)
pkgdown::build_news(preview = FALSE)
pkgdown::build_search()

message("\npkgdown site built at: ", normalizePath("docs"))
