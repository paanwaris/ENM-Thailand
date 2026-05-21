# =============================================================================
# build_pkgdown_site.R
# -----------------------------------------------------------------------------
# Build the public pkgdown website for the ENM-Thailand workshop.
#
# Prerequisites:
#   * The shared inputs (data/processed/, temporal_rasters/) already exist —
#     knit Day2a_Data_Download.Rmd first so every vignette can run end-to-end.
#   * `pkgdown` and the three workshop packages are installed.
#
# What this script does:
#   1. Installs missing build-time dependencies (pkgdown, knitr, rmarkdown).
#   2. Calls `pkgdown::build_site()` against the wrappers in `vignettes/`,
#      which use `child = "../DayX_*.Rmd"` chunks to render the canonical
#      top-level Rmds — there is no duplicated source.
#   3. Writes the rendered HTML into `docs/`.
#
# Deploy via GitHub Pages: push the contents of `docs/` to the `gh-pages`
# branch with `pkgdown::deploy_to_branch()`, or wire up the canonical
# `usethis::use_pkgdown_github_pages()` workflow once.
# =============================================================================

needed  <- c("pkgdown", "knitr", "rmarkdown")
missing <- needed[!needed %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing, dependencies = TRUE)

# Sanity-check that the canonical Rmds exist at the repo root
root_rmds <- c(
  "Day1_nicheR.Rmd",
  "Day2a_Data_Download.Rmd",
  "Day2b_Bean_Processing.Rmd",
  "Day3_TemporalModelR.Rmd"
)

missing_rmds <- root_rmds[!file.exists(root_rmds)]
if (length(missing_rmds)) {
  stop("Missing canonical Rmds at the repository root: ",
       paste(missing_rmds, collapse = ", "))
}

# Build the site. `new_process = FALSE` keeps the inherited environment so
# packages like geodata don't have to be re-resolved per vignette.
pkgdown::build_site(
  override    = list(destination = "docs"),
  new_process = FALSE,
  install     = FALSE
)

message("\npkgdown site built at: ", normalizePath("docs"))
