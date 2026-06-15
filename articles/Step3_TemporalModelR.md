# Step 3 — Temporally-explicit SDMs with TemporalModelR

## Overview

> **Workshop navigation**
>
> - [**Pre-workshop** — Downloading the workshop
>   data.](https://paanwaris.github.io/ENM-Thailand/articles/Pre-workshop.html)
>   Get the Thailand boundary, WorldClim climate layers, and GBIF Sambar
>   records before the workshop starts.
> - [**Step 1: nicheR** — Virtual species with
>   nicheR.](https://paanwaris.github.io/ENM-Thailand/articles/Step1_nicheR.html)
>   Build a virtual species niche and map it.
> - [**Step 2: bean** — Environmental thinning with
>   bean.](https://paanwaris.github.io/ENM-Thailand/articles/Step2_bean.html)
>   Reduce environmental sampling bias before fitting a real niche.
> - **Step 3: TemporalModelR — Temporally-explicit models with
>   TemporalModelR.** Link species occurrences with environmental data
>   through time.

`TemporalModelR` builds **temporally-explicit** species distribution
models. Each occurrence is paired with the environmental conditions *at
the time it was recorded*, rather than against a long-term average.
Models can then project predictions as time-series surfaces and detect
pixel-level temporal trajectories (always-suitable, never-suitable,
increasing, decreasing, or fluctuating).

> **Data.** Step 3 uses the pre-prepared annual rasters that live under
> `data/processed/temporal_split/` — one precipitation
> (`prec_YYYY.tif`), one temperature (`temp_YYYY.tif`, derived from
> land-surface temperature), and one NDVI (`NDVI_YYYY.tif`) file per
> year. **No new climate data is downloaded here** — we simply read the
> local files and stack them. The `Pre-workshop.Rmd` notebook produced
> the Sambar deer GBIF table and the Thailand polygon we use alongside
> them.

In this Step we will:

1.  Read the per-year precipitation, temperature, and NDVI rasters and
    align them to a common grid.
2.  Build a temporal-mean E-space, thin the GBIF Sambar deer occurrences
    **once** with `bean` against it, then rarefy with
    [`spatiotemporal_rarefaction()`](https://cjhughes926.github.io/TemporalModelR/reference/spatiotemporal_rarefaction.html),
    and pair each surviving record with the environment it actually
    experienced using
    [`temporally_explicit_extraction()`](https://cjhughes926.github.io/TemporalModelR/reference/temporally_explicit_extraction.html).
3.  Z-score the rasters, summarise their temporal mean, build
    spatiotemporal CV folds, and create fold-stratified pseudoabsences.
4.  Fit all four `TemporalModelR` algorithms
    ([`build_temporal_glm()`](https://cjhughes926.github.io/TemporalModelR/reference/build_temporal_glm.html),
    [`build_temporal_gam()`](https://cjhughes926.github.io/TemporalModelR/reference/build_temporal_gam.html),
    [`build_temporal_rf()`](https://cjhughes926.github.io/TemporalModelR/reference/build_temporal_rf.html),
    [`build_temporal_hv()`](https://cjhughes926.github.io/TemporalModelR/reference/build_temporal_hv.html))
    against the same partition and project predictions across years with
    [`generate_spatiotemporal_predictions()`](https://cjhughes926.github.io/TemporalModelR/reference/generate_spatiotemporal_predictions.html).
5.  Postprocess **every model** with
    [`summarize_raster_outputs()`](https://cjhughes926.github.io/TemporalModelR/reference/summarize_raster_outputs.html),
    [`analyze_temporal_patterns()`](https://cjhughes926.github.io/TemporalModelR/reference/analyze_temporal_patterns.html),
    and
    [`analyze_trends_by_spatial_unit()`](https://cjhughes926.github.io/TemporalModelR/reference/analyze_trends_by_spatial_unit.html)
    so consensus surfaces, trajectory classifications, and provincial
    summaries are directly comparable across algorithms.

All function names and arguments below follow the official package
documentation: <https://cjhughes926.github.io/TemporalModelR/>

## Install and load packages

``` r

# ---- CRAN dependencies -----------------------------------------------------
cran_pkgs <- c(
  "terra",
  "sf",
  "dplyr",
  "readr",
  "viridis",
  "mgcv",          # backend for build_temporal_gam()
  "randomForest",  # backend for build_temporal_rf()
  "hypervolume",   # backend for build_temporal_hv() (presence-only)
  "fastcpd",
  "rmarkdown",
  "knitr",
  "remotes",
  "geodata",
  "bean"
)
to_install <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

# ---- TemporalModelR from GitHub --------------------------------------------
if (!requireNamespace("TemporalModelR", quietly = TRUE)) {
  remotes::install_github("CJHughes926/TemporalModelR")
}
```

``` r

library(TemporalModelR)  # temporally-explicit SDM pipeline
library(bean)            # per-year environmental thinning
library(terra)           # modern raster + vector spatial backend
library(sf)              # tidy vector spatial data ("simple features")
library(dplyr)           # tidy data manipulation
library(readr)           # fast, friendly CSV reading and writing

# Make every random draw in the rest of the file reproducible.
set.seed(2026)
```

## Inputs

``` r

# Thailand polygon (from the Pre-workshop). `quiet = TRUE` suppresses the
# usual GDAL chatter so the rendered page stays clean.
thailand <- st_read("data/processed/boundary/thailand_boundary.gpkg",
                    quiet = TRUE)

# Clean GBIF Sambar deer table from the Pre-workshop. We:
#   * restrict to 2015-2025 (the range covered by every per-year raster);
#   * rename the GBIF coordinate columns to the conventional `longitude` /
#     `latitude` so they match the names bean and TemporalModelR expect;
#   * add a `presence = 1L` column so the table is ready to be combined
#     with pseudoabsences (presence = 0) later on.
sambar <- read_csv("data/processed/gbif/sambar_thailand.csv",
                   show_col_types = FALSE) %>%
  filter(year >= 2015, year <= 2025) %>% # keep modelled years only
  rename(longitude = decimalLongitude, # GBIF -> conventional name
         latitude = decimalLatitude) %>% # GBIF -> conventional name
  mutate(presence = 1L) # every GBIF record is a presence

# How many records per year survive that filter?
table(sambar$year)
```

    ## 
    ## 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 
    ##   31   30   31   40   55   19   17   68  106  143  226

## Preprocessing

### Align rasters

[`raster_align()`](https://cjhughes926.github.io/TemporalModelR/reference/raster_align.html)
reprojects and resamples every raster found in `input_dir` against a
reference, then writes the aligned copies into `output_dir`. Because the
three variables come from different sources, they do not share an
extent, resolution, or CRS out of the box — we use one of the per-year
NDVI files as the reference grid for everything.

``` r

# Quick visual sanity-check of the three "current" layers.
prec_ex <- rast("data/processed/temporal_split/prec_2025.tif")
temp_ex <- rast("data/processed/temporal_split/temp_2025.tif")
NDVI_ex <- rast("data/processed/temporal_split/NDVI_2025.tif")

par(mfrow = c(1, 3))
plot(prec_ex, main = "Precipitation 2025")
plot(temp_ex, main = "Temperature 2025")
plot(NDVI_ex, main = "NDVI 2025")
```

![](Step3_TemporalModelR_files/figure-html/raster-align-1.png)

``` r

par(mfrow = c(1, 1))     # reset the plot region back to a single panel

# Source directory of the per-variable per-year split rasters, and target
# directory where the aligned copies will be written.
split_dir <- "data/processed/temporal_split"
aligned_dir <- "data/processed/rasters_aligned"
dir.create(aligned_dir, showWarnings = FALSE, recursive = TRUE)

# Use NDVI as the reference grid — every per-year raster will be reprojected
# / resampled to match its extent, resolution, and CRS.
ref_raster <- NDVI_ex

raster_align(
  input_dir = split_dir, # where the per-year rasters currently live
  output_dir = aligned_dir, # where the aligned copies will be written
  reference_raster = ref_raster, # target grid (extent + resolution + CRS)
  overwrite = TRUE # replace any aligned files left from a prior run
)
```

### Build a temporal-mean E-space for thinning

Before any thinning, we build a single **temporal-mean** reference
E-space that summarises the environment Sambar deer experienced across
the entire 2010–2025 window. We compute the per-pixel mean across years
for each of the three aligned variables, then stack the three mean
surfaces into one 3-layer `SpatRaster`.

This temporal-mean stack will be the **reference environmental space**
that `bean` uses to thin every year. Using a single, period-wide E-space
(rather than the year-specific raster of each occurrence) keeps the
thinning grid *consistent across years*: a 2014 record and a 2020 record
that fall in the same temporal-mean pod are treated as environmentally
redundant, even though their individual years’ rasters may differ
slightly.

``` r

# List every aligned per-year raster, then split by variable name.
aligned_files <- list.files(aligned_dir, pattern = "\\.tif$", full.names = TRUE)

# Group the file list by variable name (prec_ / temp_ / NDVI_) so each
# group becomes a multi-layer SpatRaster spanning every year.
prec_files <- file.path(aligned_dir,
                        grep("prec_", basename(aligned_files), value = TRUE))
temp_files <- file.path(aligned_dir,
                        grep("temp_", basename(aligned_files), value = TRUE))
NDVI_files <- file.path(aligned_dir,
                        grep("NDVI_", basename(aligned_files), value = TRUE))

# Per-pixel mean across years for each variable — terra::mean() over a
# multi-layer SpatRaster collapses the time dimension into a single layer.
mean_prec <- terra::mean(rast(prec_files), na.rm = TRUE)
mean_temp <- terra::mean(rast(temp_files), na.rm = TRUE)
mean_NDVI <- terra::mean(rast(NDVI_files), na.rm = TRUE)

# Stack the three mean surfaces into one 3-layer reference E-space, and
# rename the layers to the variable names bean will look up.
env_ref <- c(mean_prec, mean_temp, mean_NDVI)
names(env_ref) <- c("mean_prec", "mean_temp", "mean_NDVI")

par(mfrow = c(1, 3))
plot(mean_prec, main = "Mean precipitation (2010–2025)")
plot(mean_temp, main = "Mean temperature (2010–2025)")
plot(mean_NDVI, main = "Mean NDVI (2010–2025)")
```

![](Step3_TemporalModelR_files/figure-html/temporal-mean-env-1.png)

``` r

par(mfrow = c(1, 1)) # reset the plot region back to a single panel
```

Reading the panels: each map is the per-pixel average of the
corresponding variable across the 16 modelled years, in the variable’s
natural units (mm for precipitation, °C-equivalent units for
land-surface-temperature-derived `temp`, and a unitless 0–1 index for
NDVI). Bright patches stand for places that were *consistently* wet,
warm, or vegetated; dark patches stand for places that were
*consistently* dry, cool, or sparse. This is the environmental landscape
`bean` will use to define the thinning grid.

### Environmental thinning with `bean`

We reduce **environmental sampling bias** in the Sambar deer occurrences
using `bean`, with the temporal-mean stack `env_ref` from the previous
section as the reference E-space. Because the reference E-space is a
single, period-wide summary (not year-specific), we run `bean` **once on
the full occurrence table** rather than year by year — every record is
judged against the same environmental landscape, so two records that
fall in the same pod are environmentally redundant regardless of which
year they were collected in.

We use the **`"sheather-jones"`** bandwidth selector inside
[`find_env_resolution()`](https://paanwaris.github.io/bean/reference/find_env_resolution.html).
Sheather–Jones is a data-adaptive, plug-in selector (Sheather & Jones,
1991) that is robust to the non-Gaussian shapes typical of real
occurrence data — it is the modern recommended default and is what the
workshop uses everywhere bean appears. We pair it with the stochastic
[`thin_env_nd()`](https://paanwaris.github.io/bean/reference/thin_env_nd.html)
so each pod retains one of the original GBIF coordinates rather than a
synthetic centroid.

Every column of the original `sambar` table is preserved on the output.
After thinning bean returns only the longitude/latitude pair, so we
recover the full original attributes (including the `year` column the
temporal modelling needs) by matching on the coordinate pair.

``` r

# Pass only the coordinate pair into bean; the full attribute table is
# recovered below by joining the thinned coordinates back onto `sambar`.
sambar_xy <- sambar[, c("longitude", "latitude")]

# Extract the three temporal-mean env values at every occurrence, and
# z-score them so all three axes share the same numeric scale.
prepped <- prepare_bean(
  data = sambar_xy,
  env_rasters = env_ref, # temporal-mean 3-layer stack
  longitude = "longitude",
  latitude = "latitude",
  transform = "scale"
)

# Sheather-Jones plug-in bandwidth as the cell width (one width per
# variable, in z-score units because prepare_bean() scaled the data).
res <- find_env_resolution(
  data = prepped,
  env_vars = c("mean_prec", "mean_temp", "mean_NDVI"),
  method = "sheather-jones"
)

# Build the environmental grid using `res$suggested_resolution` and
# randomly retain one of the original occurrences in every occupied pod.
thin_out <- thin_env_nd(
  data = prepped,
  env_vars = c("mean_prec", "mean_temp", "mean_NDVI"),
  grid_resolution = res$suggested_resolution
)

# Recover every original column by matching on the coordinate pair.
keep_key <- paste(thin_out$thinned_data$longitude,
                  thin_out$thinned_data$latitude)
sambar_thinned <- sambar[
  paste(sambar$longitude, sambar$latitude) %in% keep_key, , drop = FALSE
]

# Before / after counts per year — verifies that the year column survived.
data.frame(
  year   = sort(unique(sambar$year)),
  before = as.integer(table(sambar$year)),
  after  = as.integer(table(factor(sambar_thinned$year,
                                   levels = sort(unique(sambar$year)))))
)
```

    ##    year before after
    ## 1  2015     31     7
    ## 2  2016     30     6
    ## 3  2017     31     7
    ## 4  2018     40     5
    ## 5  2019     55     3
    ## 6  2020     19     5
    ## 7  2021     17     3
    ## 8  2022     68     8
    ## 9  2023    106    22
    ## 10 2024    143    16
    ## 11 2025    226    29

``` r

# Use the thinned table for everything that follows.
sambar <- sambar_thinned
```

Every pod of the temporal-mean E-space now contributes at most one
Sambar deer occurrence, so the modelling table no longer over-represents
the heavily-surveyed parts of the environment.

### Spatio-temporal rarefaction

Once the per-year thinning has reduced *environmental* duplication, we
still want to reduce *spatial* duplication: keep at most one record per
raster pixel **per year**, so a repeated visit to the same place in a
later year is preserved as a temporally-independent observation.

``` r

# Working directory for everything point-related downstream.
points_dir <- "data/processed/points"
dir.create(points_dir, showWarnings = FALSE, recursive = TRUE)

# Persist the thinned table to disk so spatiotemporal_rarefaction() can
# pick it up by file path.
points_csv <- file.path(points_dir, "sambar_temporal.csv")
write_csv(sambar, points_csv)

# Use the CRS of the reference raster so points and rasters live in the
# same coordinate system throughout the pipeline.
study_crs <- st_crs(rast(ref_raster))

rare_out <- spatiotemporal_rarefaction(
  points_sp = points_csv, # input CSV of presences
  output_dir = file.path(points_dir, "rarefied"), # where rarefied tables are written
  reference_raster = ref_raster, # defines the pixel grid for rarefaction
  time_cols = "year", # column carrying the time stamp
  xcol = "longitude", # column with x coordinate
  ycol = "latitude", # column with y coordinate
  points_crs = study_crs, # CRS of the input points
  output_prefix = "Sambar_annual", # prefix for the written files
  verbose = FALSE # silence per-record progress messages
)

# The three counts below show how aggressive rarefaction was: how many
# records went in, how many survive pure spatial rarefaction, and how
# many survive spatio-temporal rarefaction (one record per pixel per
# year).
rare_out$input_points # n records that entered rarefaction
```

    ## [1] 111

``` r

rare_out$spatial_points # n surviving pure spatial rarefaction (one per pixel)
```

    ## [1] 110

``` r

rare_out$spatiotemporal_points # n surviving spatio-temporal (one per pixel per year)
```

    ## [1] 110

### Match each observation to the environment it experienced

[`temporally_explicit_extraction()`](https://cjhughes926.github.io/TemporalModelR/reference/temporally_explicit_extraction.html)
walks the aligned raster directory using `variable_patterns`: clean
variable names on the left, filename patterns with time placeholders
(`YEAR`) on the right. For each surviving occurrence, it pulls the
precipitation, temperature, and NDVI values from the raster of the
matching year — so a 2014 record is paired with the 2014 climate, not
the long-term average.

``` r

# Output directory for the extracted tables and the scaling parameters.
ext_dir <- file.path(points_dir, "extracted")

ext_out <- temporally_explicit_extraction(
  points_sp = rare_out$files_created$spatiotemporal, # rarefied points to extract at
  raster_dir = aligned_dir, # raster source folder
  variable_patterns = c( # variable -> filename pattern map
    "prec" = "prec_YEAR",
    "temp" = "temp_YEAR",
    "NDVI" = "NDVI_YEAR"
  ),
  time_cols = "year", # time-stamp column on the points
  xcol = "X", # rarefied table uses upper-case X / Y
  ycol = "Y",
  points_crs = study_crs, # CRS shared by points and rasters
  output_dir = ext_dir,
  output_prefix = "sambar_extracted",
  save_raw = TRUE, # write the raw (un-scaled) extracted table
  save_scaled = TRUE, # also write a z-scored copy
  save_scaling_params = TRUE, # write the per-variable means + SDs we'll reuse
  verbose = FALSE # silence per-record extraction progress
)

# Inspect the first few rows of the raw (un-scaled) extracted table.
head(ext_out$raw_values)
```

    ##   year pixel_id     prec     temp      NDVI         x        y
    ## 1 2025     1698 1636.465 22.73027 0.7164522  99.38985 18.82608
    ## 2 2022     1785 1563.137 28.29210 0.7183695  98.82302 18.68756
    ## 3 2021     3210 1561.478 25.62316 0.7308000 101.57942 17.36530
    ## 4 2025     3574 1649.073 27.18566 0.6480218 100.87285 16.99498
    ## 5 2025     3676 1558.965 25.64610 0.5752261 101.62258 16.97002
    ## 6 2025     3677 1590.410 24.14503 0.6517609 101.69314 16.96520

``` r

# The means and SDs used to z-score the points. We will re-use these
# below so the rasters and points share the same transformation.
ext_out$scaling_params
```

    ##   variable         mean          sd
    ## 1     prec 1563.6098544 452.4398546
    ## 2     temp   28.0566521   2.6989743
    ## 3     NDVI    0.6732006   0.1062606

### Scale rasters

[`scale_rasters()`](https://cjhughes926.github.io/TemporalModelR/reference/scale_rasters.html)
z-scores every raster in `input_dir` using the per-variable means and
SDs from extraction, so that the rasters and the occurrence points share
the *same* transformation. This is critical: a GLM or RF coefficient
learned on z-scored points only makes sense when applied to z-scored
rasters.

``` r

# Destination directory for the scaled (z-scored) per-year rasters.
scaled_dir <- "data/processed/rasters_scaled"
dir.create(scaled_dir, showWarnings = FALSE, recursive = TRUE)

scale_rasters(
  input_dir = aligned_dir, # raw aligned rasters
  output_dir = scaled_dir, # write the z-scored copies here
  scaling_params_file = ext_out$files_created$scaling_params, # reuse the means + SDs from extraction
  variable_patterns = c(
    "prec" = "prec_YEAR",
    "temp" = "temp_YEAR",
    "NDVI" = "NDVI_YEAR"
  ),
  time_cols = "year",
  overwrite = TRUE, # replace any scaled files left from a previous run
  verbose = FALSE # silence per-raster scaling progress
)
```

### Spatiotemporal partitioning

[`spatiotemporal_partition()`](https://cjhughes926.github.io/TemporalModelR/reference/spatiotemporal_partition.html)
splits the occurrence + raster history into cross-validation folds in
**both** space and time, so a fold left out for testing is novel along
both dimensions — a much harder, more honest evaluation than random
k-fold CV.

``` r

par(mfrow = c(1, 3))

partition <- spatiotemporal_partition(
  reference_shapefile_path = thailand, # study-area polygon (Thailand)
  points_file_path = ext_out$files_created$scaled, # presence table with extracted env
  xcol = "x", # scaled table uses lower-case x / y
  ycol = "y",
  points_crs = study_crs,
  time_cols = "year",
  n_spatial_folds = 2, # how many blocks to split the country into
  n_temporal_folds = 2, # how many blocks to split the years into
  max_imbalance = 0.15, # max allowed presence imbalance between folds
  create_plot = TRUE, # produce the fold-map diagnostic
  verbose = FALSE # silence per-fold construction messages
)
```

![](Step3_TemporalModelR_files/figure-html/partition-1.png)

``` r

# Summary table: which fold has which years and which spatial block.
partition$summary
```

    ##                        parameter          value
    ## 1                    total_folds              4
    ## 2                n_spatial_folds              2
    ## 3               n_temporal_folds              2
    ## 4               n_balanced_folds              0
    ## 5                 n_random_folds              0
    ## 6                 partition_mode spatiotemporal
    ## 7                   total_points            110
    ## 8                 points_removed              0
    ## 9               pct_rows_removed              0
    ## 10           final_imbalance_pct          16.36
    ## 11 temporal_partitioning_enabled           TRUE

``` r

par(mfrow = c(1, 1)) # reset the plot region back to a single panel
```

### Generate pseudoabsences

A binary GLM or RF needs *absences* as well as presences. We generate
**fold-stratified pseudoabsences** — random non-presence locations drawn
at a 2:1 absence:presence ratio from anywhere outside a 50 km buffer
around the presences in the same fold.

``` r

absences <- generate_absences(
  partition_result = partition, # fold structure built above
  reference_shapefile_path = thailand, # absences are drawn inside this polygon
  raster_dir = scaled_dir, # env values are extracted at each absence
  variable_patterns = c(
    "prec" = "prec_YEAR",
    "temp" = "temp_YEAR",
    "NDVI" = "NDVI_YEAR"
  ),
  method = "buffer", # draw absences outside a buffer around presences
  buffer_distance = 50000, # 50 km buffer around each presence
  ratio = 2, # 2 absences per presence
  time_cols = "year",
  create_plot = TRUE, # produce a diagnostic absence map
  plot_by_fold = TRUE, # one panel per fold
  verbose = FALSE # silence per-fold progress messages
)
```

![](Step3_TemporalModelR_files/figure-html/absences-1.png)![](Step3_TemporalModelR_files/figure-html/absences-2.png)![](Step3_TemporalModelR_files/figure-html/absences-3.png)![](Step3_TemporalModelR_files/figure-html/absences-4.png)![](Step3_TemporalModelR_files/figure-html/absences-5.png)

``` r

# Per-fold counts of presences vs absences.
absences$summary
```

    ##   fold n_presences n_pseudoabsences ratio_achieved
    ## 1    1          27               54              2
    ## 2    2          27               54              2
    ## 3    3          32               64              2
    ## 4    4          24               48              2

## Modelling

`TemporalModelR` ships four model families, all with the same
`partition_result` + (optionally) `pseudoabsence_result` interface so
swapping algorithms is a one-chunk change:

- **GLM** — generalised linear model (`build_temporal_glm`). Linear /
  polynomial responses, easiest to interpret.
- **GAM** — generalised additive model (`build_temporal_gam`, via
  `mgcv`). Smooth, data-adaptive nonlinear responses while still
  per-variable interpretable.
- **Random forest** — ensemble of decision trees (`build_temporal_rf`,
  via `randomForest`). Captures non-linear interactions automatically;
  less interpretable than GLM/GAM.
- **Hypervolume** — *n*-dimensional kernel-density envelope
  (`build_temporal_hv`, via `hypervolume`). **Presence-only** — no
  pseudoabsences needed and no probability threshold.

We will fit all four against the same partition so the per-fold metrics
are directly comparable. Three of them (GLM, GAM, RF) consume the same
presence/absence table; only the hypervolume omits
`pseudoabsence_result`.

### GLM (linear + a quadratic on temperature)

The temperature niche is expected to peak at some intermediate value and
fall away on either side, so we add a quadratic term on `temp` to
capture that hump. The `tss` threshold method picks the cutoff that
maximises the True Skill Statistic on each fold’s test set.

``` r

glm_res <- build_temporal_glm(
  partition_result = partition, # CV folds in space + time
  pseudoabsence_result = absences, # presences + absences per fold
  model_formula = ~ temp + I(temp^2) + prec + NDVI, # quadratic on temp, linear on prec / NDVI
  link = "logit", # binary GLM with logit link
  threshold_method = "tss", # pick threshold maximising True Skill Statistic
  output_dir = "outputs/GLM_Models",
  create_plot = FALSE,
  overwrite = TRUE,
  time_cols = "year",
  verbose = FALSE
)

# Out-of-fold E-space test metrics — one row per fold.
glm_res$fold_test_metrics
```

    ##   Fold Threshold Testing_TP Testing_FN Testing_TN Testing_FP Sensitivity
    ## 1    1      0.32         17         10         26         28      0.6296
    ## 2    2      0.42         17         10         43         11      0.6296
    ## 3    3      0.41         20         12         41         20      0.6250
    ## 4    4      0.37         14         10         23         24      0.5833
    ##   Specificity    TSS  Kappa    AUC
    ## 1      0.4815 0.1111 0.0952 0.5734
    ## 2      0.7963 0.4259 0.4220 0.7277
    ## 3      0.6721 0.2971 0.2805 0.7172
    ## 4      0.4894 0.0727 0.0636 0.5949

### GAM (smooth nonlinear responses)

A GAM lets each predictor have a *smooth, data-adaptive* response curve
without hand-coding polynomials. We wrap every predictor in `s()`
(thin-plate regression spline) and let `mgcv` choose the basis
dimension. The link function, threshold rule, and fold structure are the
same as the GLM, so the two are directly comparable.

``` r

gam_res <- build_temporal_gam(
  partition_result = partition, # CV folds in space + time
  pseudoabsence_result = absences, # presences + absences per fold
  model_formula = ~ s(temp) + s(prec) + s(NDVI), # `s()` = univariate smooth (mgcv thin-plate spline)
  link = "logit", # binary GAM with logit link
  gam_params = list(method = "REML"), # smoothing-parameter selection rule (mgcv default REML)
  threshold_method = "tss", # pick threshold maximising True Skill Statistic
  output_dir = "outputs/GAM_Models",
  create_plot = FALSE, # skip the smooth-response curves to keep output compact
  overwrite = TRUE,
  time_cols = "year",
  verbose = FALSE # silence per-fold smoothing progress
)

# Out-of-fold E-space test metrics — one row per fold.
gam_res$fold_test_metrics
```

    ##   Fold Threshold Testing_TP Testing_FN Testing_TN Testing_FP Sensitivity
    ## 1    1      0.35         16         11         28         26      0.5926
    ## 2    2      0.45         10         17         48          6      0.3704
    ## 3    3      0.29         23          9         31         30      0.7188
    ## 4    4      0.42         11         13         35         12      0.4583
    ##   Specificity    TSS  Kappa    AUC
    ## 1      0.5185 0.1111 0.0976 0.5748
    ## 2      0.8889 0.2593 0.2887 0.7181
    ## 3      0.5082 0.2269 0.1963 0.7059
    ## 4      0.7447 0.2030 0.2051 0.5417

### Random forest (flexible ensemble)

A 500-tree random forest. RFs capture non-linear interactions
automatically and usually outperform the GLM/GAM at the cost of a less
interpretable model.

``` r

rf_res <- build_temporal_rf(
  partition_result = partition, # CV folds in space + time
  pseudoabsence_result = absences, # presences + absences per fold
  model_vars = c("temp", "prec", "NDVI"), # RF needs predictor names, not a formula
  rf_params = list(ntree = 500), # number of decision trees grown per fold
  threshold_method = "tss", # pick threshold maximising True Skill Statistic
  output_dir = "outputs/RF_Models",
  create_plot = FALSE, # skip variable-importance plot
  overwrite = TRUE,
  time_cols = "year",
  verbose = FALSE # silence per-fold tree-growth progress
)

# Out-of-fold E-space test metrics — one row per fold.
rf_res$fold_test_metrics
```

    ##   Fold Threshold Testing_TP Testing_FN Testing_TN Testing_FP Sensitivity
    ## 1    1      0.34         15         12         32         22      0.5556
    ## 2    2      0.34         13         14         35         19      0.4815
    ## 3    3      0.33         23          9         39         22      0.7188
    ## 4    4      0.36         10         14         33         14      0.4167
    ##   Specificity    TSS  Kappa    AUC
    ## 1      0.5926 0.1481 0.1356 0.6063
    ## 2      0.6481 0.1296 0.1239 0.6070
    ## 3      0.6393 0.3581 0.3266 0.7316
    ## 4      0.7021 0.1188 0.1188 0.6055

### Hypervolume (presence-only kernel density)

The hypervolume builds an *n*-dimensional kernel-density envelope around
the presence points alone — it does not need pseudoabsences and does not
produce a probability scale, so there is no `pseudoabsence_result` or
`threshold_method` argument. Only sensitivity-based metrics are reported
(specificity / AUC / TSS require negative training data).

> **Heads-up.** Hypervolume fitting can be slow (several minutes per
> fold on real data) because it samples millions of random points per
> envelope. If you are running the workshop on a laptop, you can comment
> this chunk out and skip the HV columns in the plots below.

``` r

hv_res <- build_temporal_hv(
  partition_result = partition, # CV folds in space + time (presences only)
  model_vars = c("temp", "prec", "NDVI"), # predictors defining the envelope axes
  method = "gaussian", # Gaussian KDE envelope (alternative: "svm" = one-class SVM)
  hypervolume_params = list(), # accept hypervolume::hypervolume_gaussian() defaults
  output_dir = "outputs/HV_Models",
  create_plot = FALSE, # skip the per-fold pairplots
  verbose = FALSE # silence per-fold KDE-construction progress
)

# Presence-only E-space metrics: sensitivity, envelope volume, per-fold overlap.
hv_res$fold_test_metrics
```

    ##   Fold N_Test E_Volume Testing_TP Testing_FN Sensitivity
    ## 1    1     30  87.5480         30          0      1.0000
    ## 2    2     29  85.9348         29          0      1.0000
    ## 3    3     34  59.0094         30          4      0.8824
    ## 4    4     26  62.5078         21          5      0.8077

``` r

hv_res$volumes
```

    ##    fold1    fold2    fold3    fold4 
    ## 87.54798 85.93478 59.00938 62.50780

## Project predictions through time

[`generate_spatiotemporal_predictions()`](https://cjhughes926.github.io/TemporalModelR/reference/generate_spatiotemporal_predictions.html)
walks every requested time step, loads the matching rasters, applies the
fitted models, and writes **one fold-vote consensus raster per time
step**. Every cell of the output carries a *count* (`0`, `1`, …, *N*) of
how many of the *N* folds voted that the pixel was suitable in that time
step — these are the package’s native fold-vote outputs, which we will
visualise with the package’s default
[`terra::plot()`](https://rspatial.github.io/terra/reference/plot.html)
styling.

The call signature is identical for all four model types; `model_type`
inside `model_result` tells the function which projection logic to
apply. The hypervolume projection simply omits `pseudoabsence_result`
because the model was never fit against absences.

Key arguments shared by every call below:

- **`raster_dir`** — the folder of *scaled* per-year predictor rasters
  written by
  [`scale_rasters()`](https://cjhughes926.github.io/TemporalModelR/reference/scale_rasters.html).
  Predictions must use the same transformation as the training points.
- **`variable_patterns`** — a named vector mapping each predictor name
  (left) to its filename pattern (right) with the `YEAR` placeholder, so
  `temp_YEAR` resolves to `temp_2015.tif` for year 2015.
- **`time_steps`** — a `data.frame` giving every (year, …) combination
  to project. We use 2015–2025 to match the `sambar` filter and the 3 ×
  4 plot grid below.
- **`output_dir`** — where the per-fold per-time-step prediction rasters
  get written. Each model writes to its own folder so the
  post-processing later can address them separately.

``` r

time_steps <- data.frame(year = 2015:2025) # 11 years -> fits a 3 x 4 plot grid

# A single named vector used by every projection call.
var_patterns <- c("prec" = "prec_YEAR",
                  "temp" = "temp_YEAR",
                  "NDVI" = "NDVI_YEAR")

# GLM projection: walk through every (year) time step and apply each fold's
# fitted GLM to the matching scaled rasters, writing one fold-vote raster
# per fold per year into output_dir.
glm_preds <- generate_spatiotemporal_predictions(
  partition_result = partition, # CV fold geometry built earlier
  model_result = glm_res, # the fitted GLM (one model per fold)
  pseudoabsence_result = absences, # tells the function which area each fold "owns"
  raster_dir = scaled_dir, # folder of z-scored per-year predictor rasters
  variable_patterns = var_patterns, # maps predictor name -> filename pattern (with YEAR)
  time_cols = "year", # name of the time column on the points
  time_steps = time_steps, # the data.frame of years (2015:2025) to project
  output_dir = "outputs/glm_predictions", # per-fold per-year rasters land here
  overwrite = TRUE, # replace any rasters left from a prior run
  verbose = FALSE # silence the per-fold progress messages
)

# GAM projection: same inputs as the GLM call above, only `model_result`
# and `output_dir` change. The function reads `model_type` inside `gam_res`
# and automatically uses the correct projection logic for a GAM.
gam_preds <- generate_spatiotemporal_predictions(
  partition_result = partition, # same fold geometry
  model_result = gam_res, # the fitted GAM (one per fold)
  pseudoabsence_result = absences, # same fold-stratified absences
  raster_dir = scaled_dir, # same scaled per-year predictor rasters
  variable_patterns = var_patterns, # same name -> filename map
  time_cols = "year", # same time-column name
  time_steps = time_steps, # same 11-year grid
  output_dir = "outputs/gam_predictions", # write into a *separate* folder per model
  overwrite = TRUE, # replace previous-run rasters
  verbose = FALSE # silence per-fold progress
)

# Random forest projection: identical interface; the function dispatches to
# the RF projection logic via `model_type` stored inside `rf_res`.
rf_preds <- generate_spatiotemporal_predictions(
  partition_result = partition, # same fold geometry
  model_result = rf_res, # the fitted RF (500 trees per fold)
  pseudoabsence_result = absences, # same fold-stratified absences
  raster_dir = scaled_dir, # same scaled per-year predictor rasters
  variable_patterns = var_patterns, # same name -> filename map
  time_cols = "year", # same time-column name
  time_steps = time_steps, # same 11-year grid
  output_dir = "outputs/rf_predictions", # separate output folder for RF
  overwrite = TRUE, # replace previous-run rasters
  verbose = FALSE # silence per-fold progress
)

# Hypervolume projection: NOTE the missing `pseudoabsence_result` argument.
# Hypervolume is presence-only, so the function never needs absences when
# projecting an HV model — passing them in is unnecessary.
hv_preds <- generate_spatiotemporal_predictions(
  partition_result = partition, # same fold geometry
  model_result = hv_res, # the fitted hypervolume (Gaussian KDE per fold)
  raster_dir = scaled_dir, # same scaled per-year predictor rasters
  variable_patterns = var_patterns, # same name -> filename map
  time_cols = "year", # same time-column name
  time_steps = time_steps, # same 11-year grid
  output_dir = "outputs/hv_predictions", # separate output folder for HV
  overwrite = TRUE, # replace previous-run rasters
  verbose = FALSE # silence per-fold progress
)
```

### Visualise the consensus maps

For each model we stack the per-time-step prediction rasters straight
from `preds$prediction_files` and pass them to
[`terra::plot()`](https://rspatial.github.io/terra/reference/plot.html)
— the package’s default visualisation. Every raster cell is a
**fold-vote count** (0 to *N*) of how many folds agreed the pixel was
suitable in that year. With 11 modelled years (2015–2025), each stack
lays out as a 3 × 4 grid.

``` r

# Load all per-year GLM prediction rasters into a single multi-layer
# SpatRaster. `terra::rast()` accepts a character vector of file paths
# and stacks them in the order given.
glm_pred_stack <- terra::rast(glm_preds$prediction_files)

# Set human-readable layer names (e.g. "Prediction_..._2015") by stripping
# the `.tif` extension from each filename. These names become the panel
# titles in the plot below.
names(glm_pred_stack) <- sub("\\.tif$", "", basename(glm_preds$prediction_files))

# Draw the 11 yearly fold-vote maps as a 3 x 4 grid.
# Arguments:
#   - First positional argument: the multi-layer SpatRaster to plot.
#   - nr = 3, nc = 4: 3 rows by 4 columns of panels.
#   - mar = c(1.0, 1.0, 1.5, 3.0): inner margins for each panel
#     (bottom, left, top, right) in lines of text; the wider right margin
#     reserves room for the per-panel legend strip.
#   - legend = FALSE: suppress per-panel legends so the grid stays tidy
#     (we explain the colour scale in the prose below).
terra::plot(glm_pred_stack, nr = 3, nc = 4,
            mar    = c(1.0, 1.0, 1.5, 3.0),
            legend = FALSE)
```

![](Step3_TemporalModelR_files/figure-html/plot-glm-1.png)

``` r

# Same three-step pattern as the GLM chunk above: build a multi-layer
# SpatRaster from the GAM prediction files, label each layer, and draw a
# 3 x 4 grid of yearly fold-vote maps.
gam_pred_stack <- terra::rast(gam_preds$prediction_files)
names(gam_pred_stack) <- sub("\\.tif$", "", basename(gam_preds$prediction_files))
terra::plot(gam_pred_stack, nr = 3, nc = 4,
            mar    = c(1.0, 1.0, 1.5, 3.0),
            legend = FALSE)
```

![](Step3_TemporalModelR_files/figure-html/plot-gam-1.png)

``` r

# Same pattern, but for the random forest prediction files.
rf_pred_stack <- terra::rast(rf_preds$prediction_files)
names(rf_pred_stack) <- sub("\\.tif$", "", basename(rf_preds$prediction_files))
terra::plot(rf_pred_stack, nr = 3, nc = 4,
            mar    = c(1.0, 1.0, 1.5, 3.0),
            legend = FALSE)
```

![](Step3_TemporalModelR_files/figure-html/plot-rf-1.png)

``` r

# Same pattern, but for the hypervolume prediction files.
hv_pred_stack <- terra::rast(hv_preds$prediction_files)
names(hv_pred_stack) <- sub("\\.tif$", "", basename(hv_preds$prediction_files))
terra::plot(hv_pred_stack, nr = 3, nc = 4,
            mar    = c(1.0, 1.0, 1.5, 3.0),
            legend = FALSE)
```

![](Step3_TemporalModelR_files/figure-html/plot-hv-1.png)

Reading these grids: each panel is one year (2015 in the top-left, 2025
in the bottom-right of each 3 × 4 grid). Brighter pixels mean more folds
agreed the pixel was suitable in that year; darker pixels mean fewer or
no folds did. Any *shift* in the bright patches across panels is the
model’s view of how Sambar deer habitat changed during 2015–2025.
Comparing the four model grids side-by-side highlights where the
algorithms agree on the niche and where they disagree.

### Per-model performance diagnostics

[`plot_model_assessment()`](https://cjhughes926.github.io/TemporalModelR/reference/plot_model_assessment.html)
is the package’s default diagnostic plot. It produces per-fold
time-series of percent suitable, sensitivity (and, for presence/absence
models, specificity), and the cumulative-binomial-probability (CBP)
test, overlaid with the overall E-space metrics from
`$fold_test_metrics` when `model_result` is supplied. We run it once per
model below so the four algorithms are directly comparable.

Key arguments:

- **`predictions`** — the object returned by
  [`generate_spatiotemporal_predictions()`](https://cjhughes926.github.io/TemporalModelR/reference/generate_spatiotemporal_predictions.html)
  for that model.
- **`time_column`** — the name(s) of the time column(s) in the
  partition. Single annual axis here, so `"year"`.
- **`model_result`** — the fitted-model object; supplying it overlays
  the E-space `$fold_test_metrics` as reference lines on the
  per-time-step series.

``` r

# Draw the GLM's per-fold diagnostic time series.
# - predictions:  the object returned by generate_spatiotemporal_predictions()
#                 for this model; contains the per-time-step G-space metrics.
# - time_column:  name of the time axis used for the x-axis of the panels.
# - model_result: the fitted-model object; supplying it overlays the
#                 overall E-space fold metrics as horizontal reference lines.
plot_model_assessment(
  predictions = glm_preds,
  time_column = "year",
  model_result = glm_res
)
```

![](Step3_TemporalModelR_files/figure-html/assess-glm-1.png)![](Step3_TemporalModelR_files/figure-html/assess-glm-2.png)![](Step3_TemporalModelR_files/figure-html/assess-glm-3.png)![](Step3_TemporalModelR_files/figure-html/assess-glm-4.png)![](Step3_TemporalModelR_files/figure-html/assess-glm-5.png)![](Step3_TemporalModelR_files/figure-html/assess-glm-6.png)

``` r

# Same diagnostic for the GAM. Arguments behave identically to the GLM
# call above, only the `predictions` and `model_result` objects differ.
plot_model_assessment(
  predictions = gam_preds,
  time_column = "year",
  model_result = gam_res
)
```

![](Step3_TemporalModelR_files/figure-html/assess-gam-1.png)![](Step3_TemporalModelR_files/figure-html/assess-gam-2.png)![](Step3_TemporalModelR_files/figure-html/assess-gam-3.png)![](Step3_TemporalModelR_files/figure-html/assess-gam-4.png)![](Step3_TemporalModelR_files/figure-html/assess-gam-5.png)![](Step3_TemporalModelR_files/figure-html/assess-gam-6.png)

``` r

# Same diagnostic for the random forest.
plot_model_assessment(
  predictions = rf_preds,
  time_column = "year",
  model_result = rf_res
)
```

![](Step3_TemporalModelR_files/figure-html/assess-rf-1.png)![](Step3_TemporalModelR_files/figure-html/assess-rf-2.png)![](Step3_TemporalModelR_files/figure-html/assess-rf-3.png)![](Step3_TemporalModelR_files/figure-html/assess-rf-4.png)![](Step3_TemporalModelR_files/figure-html/assess-rf-5.png)![](Step3_TemporalModelR_files/figure-html/assess-rf-6.png)

``` r

# Same diagnostic for the hypervolume. The function detects that the
# model is presence-only and automatically drops the specificity panel
# (no negative training data => no specificity).
plot_model_assessment(
  predictions = hv_preds,
  time_column = "year",
  model_result = hv_res
)
```

![](Step3_TemporalModelR_files/figure-html/assess-hv-1.png)![](Step3_TemporalModelR_files/figure-html/assess-hv-2.png)![](Step3_TemporalModelR_files/figure-html/assess-hv-3.png)![](Step3_TemporalModelR_files/figure-html/assess-hv-4.png)

The hypervolume diagnostics show only the sensitivity, percent-suitable,
and CBP panels — specificity is undefined for a presence-only model.

## Postprocessing

Every model produced its own folder of per-fold per-year prediction
rasters during projection (`outputs/glm_predictions/`,
`outputs/gam_predictions/`, `outputs/rf_predictions/`,
`outputs/hv_predictions/`). The three post-processing functions below —
[`summarize_raster_outputs()`](https://cjhughes926.github.io/TemporalModelR/reference/summarize_raster_outputs.html),
[`analyze_temporal_patterns()`](https://cjhughes926.github.io/TemporalModelR/reference/analyze_temporal_patterns.html),
and
[`analyze_trends_by_spatial_unit()`](https://cjhughes926.github.io/TemporalModelR/reference/analyze_trends_by_spatial_unit.html)
— are model-agnostic: they read a prediction folder and operate on
whatever the prediction rasters contain. We therefore run **one call per
model** in each step so the consensus surfaces, trajectory
classifications, and provincial summaries are directly comparable across
GLM / GAM / RF / Hypervolume.

### Consensus surfaces across folds

[`summarize_raster_outputs()`](https://cjhughes926.github.io/TemporalModelR/reference/summarize_raster_outputs.html)
reads the per-time-step prediction rasters from `predictions_dir` and
applies a **consensus rule** across folds to flag pixels that *at least
`consensus` folds* agreed were suitable. Setting `consensus = 2` is a
majority rule for a 4-fold partition. It returns a binary
`consensus_stack` (one layer per year) plus a `frequency_raster` (the
proportion of years each pixel was suitable).

Key arguments:

- **`predictions_dir`** — folder of fold/year prediction rasters from
  [`generate_spatiotemporal_predictions()`](https://cjhughes926.github.io/TemporalModelR/reference/generate_spatiotemporal_predictions.html).
- **`output_dir`** — where to write the binary stack and frequency
  raster.
- **`consensus`** — minimum number of folds that must agree before a
  pixel is called suitable; here, `2` of `4`.

``` r

# Collapse the GLM's per-fold per-year prediction rasters into two
# products: a binary `consensus_stack` (one layer per year) and a single
# 0-1 `frequency_raster` (proportion of years suitable).
glm_consensus <- summarize_raster_outputs(
  predictions_dir = "outputs/glm_predictions", # source folder of per-fold rasters
  output_dir = "outputs/glm_consensus", # where to write the consensus/frequency outputs
  consensus = 2, # majority rule: >=2 of the 4 folds must agree
  overwrite = TRUE, # replace any outputs from a previous run
  verbose = FALSE # silence progress messages
)
```

``` r

# Same consensus rule for the GAM. Only the source / target folders change.
gam_consensus <- summarize_raster_outputs(
  predictions_dir = "outputs/gam_predictions", # source folder of per-fold rasters
  output_dir = "outputs/gam_consensus", # GAM-specific output folder
  consensus = 2, # same majority rule
  overwrite = TRUE,
  verbose = FALSE
)
```

``` r

# Same consensus rule for the random forest.
rf_consensus <- summarize_raster_outputs(
  predictions_dir = "outputs/rf_predictions", # source folder of per-fold rasters
  output_dir = "outputs/rf_consensus", # RF-specific output folder
  consensus = 2, # same majority rule
  overwrite = TRUE,
  verbose = FALSE
)
```

``` r

# Same consensus rule for the hypervolume.
hv_consensus <- summarize_raster_outputs(
  predictions_dir = "outputs/hv_predictions", # source folder of per-fold rasters
  output_dir = "outputs/hv_consensus", # HV-specific output folder
  consensus = 2, # same majority rule
  overwrite = TRUE,
  verbose = FALSE
)
```

#### Per-year binary suitability — the `consensus_stack`

Each `consensus_stack` is a multi-layer `SpatRaster` carrying **one
binary layer per year** (yellow = suitable, purple = unsuitable). Plot
it the same way we plotted the prediction stacks earlier — a 3 × 4 grid
of the 11 modelled years for the model in question.

``` r

# Draw the GLM's per-year binary suitability surfaces as a 3 x 4 grid.
# Arguments mirror the prediction-stack plot above:
#   - First positional argument: the binary multi-layer SpatRaster.
#   - nr / nc: 3 rows by 4 columns of yearly panels.
#   - mar:     inner margins per panel (bottom, left, top, right).
#   - legend = FALSE: drop per-panel legends because the values are binary
#     (0 = unsuitable, 1 = suitable) and the prose below covers the
#     interpretation.
terra::plot(glm_consensus$consensus_stack, nr = 3, nc = 4,
            mar = c(1.0, 1.0, 1.5, 3.0),
            legend = FALSE)
```

![](Step3_TemporalModelR_files/figure-html/consensus-stack-glm-1.png)

``` r

# Same plotting recipe applied to the GAM's binary consensus stack.
terra::plot(gam_consensus$consensus_stack, nr = 3, nc = 4,
            mar = c(1.0, 1.0, 1.5, 3.0),
            legend = FALSE)
```

![](Step3_TemporalModelR_files/figure-html/consensus-stack-gam-1.png)

``` r

# Same plotting recipe applied to the random forest's binary consensus stack.
terra::plot(rf_consensus$consensus_stack, nr = 3, nc = 4,
            mar = c(1.0, 1.0, 1.5, 3.0),
            legend = FALSE)
```

![](Step3_TemporalModelR_files/figure-html/consensus-stack-rf-1.png)

``` r

# Same plotting recipe applied to the hypervolume's binary consensus stack.
terra::plot(hv_consensus$consensus_stack, nr = 3, nc = 4,
            mar = c(1.0, 1.0, 1.5, 3.0),
            legend = FALSE)
```

![](Step3_TemporalModelR_files/figure-html/consensus-stack-hv-1.png)

Reading these grids: yellow pixels were called suitable by **at least 2
of the 4 folds** (the `consensus = 2` majority rule), purple pixels were
not. Compared to the raw fold-vote maps shown earlier under *Visualise
the consensus maps*, these are the **binary** version after the
consensus threshold has been applied.

#### Frequency across years

The frequency raster collapses each consensus stack across years into a
single 0–1 surface — the proportion of years a pixel was called
suitable. Drawn per model as a 2 × 2 panel for easy side-by-side
comparison.

``` r

# Set up a 2 x 2 plot grid so each model's frequency raster occupies one
# panel, making side-by-side comparison straightforward.
par(mfrow = c(2, 2))

# Each $frequency_raster is a single 0-1 SpatRaster (proportion of years
# suitable). plot() draws it with the default viridis palette and `main`
# sets the panel title.
plot(glm_consensus$frequency_raster,
     main = "GLM — proportion of years suitable")

plot(gam_consensus$frequency_raster,
     main = "GAM — proportion of years suitable")

plot(rf_consensus$frequency_raster,
     main = "RF — proportion of years suitable")

plot(hv_consensus$frequency_raster,
     main = "Hypervolume — proportion of years suitable")
```

![](Step3_TemporalModelR_files/figure-html/consensus-plot-1.png)

``` r

# Reset the plot region back to a single panel for anything that follows.
par(mfrow = c(1, 1))
```

Reading the panels: pixels close to 1 (bright) were predicted suitable
in every modelled year — the *persistent core* of the niche for that
model. Pixels close to 0 (dark) were never suitable. Intermediate values
mark areas that gained or lost suitability across 2015–2025. Differences
between the four panels show where the algorithms disagree on long-term
suitability.

### Pixel-level temporal trajectories (changepoint analysis)

[`analyze_temporal_patterns()`](https://cjhughes926.github.io/TemporalModelR/reference/analyze_temporal_patterns.html)
runs `fastcpd` per pixel and classifies each pixel’s trajectory across
years into one of six categories: never-suitable, always-suitable,
no-pattern (statistically flat), increasing, decreasing, or fluctuating.
The result is a categorical map of **how** habitat is changing rather
than just *whether* it is suitable.

Key arguments:

- **`binary_stack`** — the per-time-step binary `consensus_stack`
  produced above.
- **`summary_raster`** — the matching 0–1 frequency raster; used to
  short-circuit pixels that are always or never suitable.
- **`time_steps`** — the same `data.frame` of time-step labels used
  during projection.
- **`spatial_autocorrelation` / `estimate_time`** — both `FALSE` here to
  keep the workshop fast; enable on real analyses if you need them.

``` r

# Run changepoint detection on the GLM's binary consensus stack to
# classify each pixel's trajectory across years (never-, always-,
# increasing, decreasing, fluctuating, no-pattern).
glm_patterns <- analyze_temporal_patterns(
  binary_stack = glm_consensus$consensus_stack, # binary per-year stack from summarize_raster_outputs()
  summary_raster = glm_consensus$frequency_raster, # 0-1 frequency raster, short-circuits always/never pixels
  time_steps = time_steps, # same data.frame of years used everywhere else
  output_dir = "outputs/glm_patterns", # where the pattern raster and CSV go
  spatial_autocorrelation = FALSE, # do NOT add a neighbours covariate (faster)
  estimate_time = FALSE, # skip the runtime-estimation pre-pass (faster)
  overwrite = TRUE, # replace any patterns left from a prior run
  verbose = FALSE # silence per-pixel progress messages
)
```

![](Step3_TemporalModelR_files/figure-html/changepoints-glm-1.png)

``` r

# Same call, applied to the GAM's consensus stack. Only the input rasters
# and the output folder change.
gam_patterns <- analyze_temporal_patterns(
  binary_stack = gam_consensus$consensus_stack,
  summary_raster = gam_consensus$frequency_raster,
  time_steps = time_steps,
  output_dir = "outputs/gam_patterns",
  spatial_autocorrelation = FALSE,
  estimate_time = FALSE,
  overwrite = TRUE,
  verbose = FALSE
)
```

![](Step3_TemporalModelR_files/figure-html/changepoints-gam-1.png)

``` r

# Same call, applied to the random forest's consensus stack.
rf_patterns <- analyze_temporal_patterns(
  binary_stack = rf_consensus$consensus_stack,
  summary_raster = rf_consensus$frequency_raster,
  time_steps = time_steps,
  output_dir = "outputs/rf_patterns",
  spatial_autocorrelation = FALSE,
  estimate_time = FALSE,
  overwrite = TRUE,
  verbose = FALSE
)
```

![](Step3_TemporalModelR_files/figure-html/changepoints-rf-1.png)

``` r

# Same call, applied to the hypervolume's consensus stack.
hv_patterns <- analyze_temporal_patterns(
  binary_stack = hv_consensus$consensus_stack,
  summary_raster = hv_consensus$frequency_raster,
  time_steps = time_steps,
  output_dir = "outputs/hv_patterns",
  spatial_autocorrelation = FALSE,
  estimate_time = FALSE,
  overwrite = TRUE,
  verbose = FALSE
)
```

![](Step3_TemporalModelR_files/figure-html/changepoints-hv-1.png)

The four trajectory-class maps, side-by-side:

``` r

# Set up a 2 x 2 plot region so the four models' trajectory maps appear
# side by side.
par(mfrow = c(2, 2))

# Each $pattern is a categorical SpatRaster whose pixel values encode the
# trajectory class (1 = never-suitable, 2 = always, 3 = no pattern,
# 4 = increasing, 5 = decreasing, 6 = fluctuating).
# `main` titles each panel; `legend = FALSE` keeps the grid tidy and we
# describe the categories in the prose below.
plot(glm_patterns$pattern,
     main = "GLM — temporal trajectory class",
     legend = FALSE)

plot(gam_patterns$pattern,
     main = "GAM — temporal trajectory class",
     legend = FALSE)

plot(rf_patterns$pattern,
     main = "RF — temporal trajectory class",
     legend = FALSE)

plot(hv_patterns$pattern,
     main = "Hypervolume — temporal trajectory class",
     legend = FALSE)
```

![](Step3_TemporalModelR_files/figure-html/changepoints-plot-1.png)

``` r

# Reset the plot region back to a single panel for anything that follows.
par(mfrow = c(1, 1))
```

Reading the panels: each colour is a categorical trajectory class.
Comparing the four maps shows how much each algorithm agrees on the
*shape* of change (increasing / decreasing / fluctuating), not just its
location.

### Aggregate by spatial unit (Thailand provinces)

[`analyze_trends_by_spatial_unit()`](https://cjhughes926.github.io/TemporalModelR/reference/analyze_trends_by_spatial_unit.html)
aggregates the pixel trajectories to a province polygon set (GADM level
1 — Thailand’s provinces) so trends can be discussed region-by-region
rather than pixel-by-pixel.

Key arguments:

- **`shapefile_path`** — the polygon set to aggregate against.
- **`name_field`** — column carrying the human-readable unit name; we
  use `NAME_1` (province name in GADM).
- **`binary_stack` / `pattern_raster` / `time_decrease_raster` /
  `time_increase_raster`** — outputs from the consensus and pattern
  steps for the model in question.

``` r

# Download Thailand provinces (GADM level 1) and convert SpatVector -> sf
# so analyze_trends_by_spatial_unit() can read it.
provinces <- geodata::gadm(country = "THA", level = 1, path = "data/raw") |>
  st_as_sf()
```

``` r

# Aggregate the GLM's pixel-level trajectories to Thailand's provinces.
# The function overlays the polygons on the rasters, counts pixels in
# each trajectory class per province, and writes summary CSVs.
glm_regional <- analyze_trends_by_spatial_unit(
  shapefile_path = provinces, # sf object: Thailand provinces
  name_field = "NAME_1", # column carrying the province name
  binary_stack = glm_consensus$consensus_stack, # per-year binary suitability stack
  pattern_raster = glm_patterns$pattern, # trajectory class per pixel
  time_decrease_raster = glm_patterns$time_decrease, # year of first detected decrease (per pixel)
  time_increase_raster = glm_patterns$time_increase, # year of first detected increase (per pixel)
  time_steps = time_steps, # same data.frame of modelled years
  output_dir = "outputs/glm_regional", # GLM-specific output folder
  overwrite = TRUE, # replace previous-run outputs
  create_plot = TRUE, # produce the diagnostic plots returned in `$plots`
  verbose = FALSE # silence per-province progress
)
```

![](Step3_TemporalModelR_files/figure-html/regional-glm-1.png)![](Step3_TemporalModelR_files/figure-html/regional-glm-2.png)

``` r

# Show the diagnostic plot(s) attached to the result object: e.g.
# scatterpie map of pattern composition, per-province time series, etc.
glm_regional$plots
```

    ## $pattern_map

![](Step3_TemporalModelR_files/figure-html/regional-glm-data-1.png)

    ## 
    ## $time_series

``` r

# First three rows of the per-province summary table (one row per province,
# columns for pixel counts in each trajectory class plus percentages).
head(glm_regional$overall_summary, 3)
```

    ##         Spatial_Unit Always_Absent Always_Present No_Pattern Increasing
    ## 1      Amnat Charoen            52              0          0          0
    ## 2          Ang Thong            17              0          0          0
    ## 3 Bangkok Metropolis            32              0          0          0
    ##   Decreasing Fluctuating Failed Total_Pixels Pct_Always_Absent
    ## 1          0           0      0           52               100
    ## 2          0           0      0           17               100
    ## 3          0           0      0           32               100
    ##   Pct_Always_Present Pct_No_Pattern Pct_Increasing Pct_Decreasing
    ## 1                  0              0              0              0
    ## 2                  0              0              0              0
    ## 3                  0              0              0              0
    ##   Pct_Fluctuating Prop_Increasing Prop_Stable_Suitable Prop_Decreasing
    ## 1               0              NA                   NA               0
    ## 2               0              NA                   NA               0
    ## 3               0              NA                   NA               0
    ##   Prop_Stable_Unsuitable
    ## 1                    100
    ## 2                    100
    ## 3                    100

``` r

# Same aggregation applied to the GAM outputs. Only the input rasters and
# the output folder differ from the GLM call above.
gam_regional <- analyze_trends_by_spatial_unit(
  shapefile_path = provinces,
  name_field = "NAME_1",
  binary_stack = gam_consensus$consensus_stack,
  pattern_raster = gam_patterns$pattern,
  time_decrease_raster = gam_patterns$time_decrease,
  time_increase_raster = gam_patterns$time_increase,
  time_steps = time_steps,
  output_dir = "outputs/gam_regional",
  overwrite = TRUE,
  create_plot = TRUE,
  verbose = FALSE
)
```

![](Step3_TemporalModelR_files/figure-html/regional-gam-1.png)![](Step3_TemporalModelR_files/figure-html/regional-gam-2.png)

``` r

# Diagnostic plots and the first three rows of the per-province summary
# for the GAM (same layout as the GLM results above).
gam_regional$plots
```

    ## $pattern_map

![](Step3_TemporalModelR_files/figure-html/regional-gam-data-1.png)

    ## 
    ## $time_series

``` r

head(gam_regional$overall_summary, 3)
```

    ##         Spatial_Unit Always_Absent Always_Present No_Pattern Increasing
    ## 1      Amnat Charoen            52              0          0          0
    ## 2          Ang Thong            17              0          0          0
    ## 3 Bangkok Metropolis            30              1          0          0
    ##   Decreasing Fluctuating Failed Total_Pixels Pct_Always_Absent
    ## 1          0           0      0           52            100.00
    ## 2          0           0      0           17            100.00
    ## 3          0           0      1           32             93.75
    ##   Pct_Always_Present Pct_No_Pattern Pct_Increasing Pct_Decreasing
    ## 1               0.00              0              0              0
    ## 2               0.00              0              0              0
    ## 3               3.12              0              0              0
    ##   Pct_Fluctuating Prop_Increasing Prop_Stable_Suitable Prop_Decreasing
    ## 1               0              NA                   NA               0
    ## 2               0              NA                   NA               0
    ## 3               0               0                   50               0
    ##   Prop_Stable_Unsuitable
    ## 1                 100.00
    ## 2                 100.00
    ## 3                  96.77

``` r

# Same aggregation applied to the random forest outputs.
rf_regional <- analyze_trends_by_spatial_unit(
  shapefile_path = provinces,
  name_field = "NAME_1",
  binary_stack = rf_consensus$consensus_stack,
  pattern_raster = rf_patterns$pattern,
  time_decrease_raster = rf_patterns$time_decrease,
  time_increase_raster = rf_patterns$time_increase,
  time_steps = time_steps,
  output_dir = "outputs/rf_regional",
  overwrite = TRUE,
  create_plot = TRUE,
  verbose = FALSE
)
```

![](Step3_TemporalModelR_files/figure-html/regional-rf-1.png)![](Step3_TemporalModelR_files/figure-html/regional-rf-2.png)

``` r

# Diagnostic plots and the first three rows of the per-province summary
# for the random forest.
rf_regional$plots
```

    ## $pattern_map

![](Step3_TemporalModelR_files/figure-html/regional-rf-data-1.png)

    ## 
    ## $time_series

``` r

head(rf_regional$overall_summary, 3)
```

    ##         Spatial_Unit Always_Absent Always_Present No_Pattern Increasing
    ## 1      Amnat Charoen            20              0          0          0
    ## 2          Ang Thong            14              0          0          0
    ## 3 Bangkok Metropolis             4              0          0          0
    ##   Decreasing Fluctuating Failed Total_Pixels Pct_Always_Absent
    ## 1          0           0     32           52             38.46
    ## 2          0           0      3           17             82.35
    ## 3          0           0     28           32             12.50
    ##   Pct_Always_Present Pct_No_Pattern Pct_Increasing Pct_Decreasing
    ## 1                  0              0              0              0
    ## 2                  0              0              0              0
    ## 3                  0              0              0              0
    ##   Pct_Fluctuating Prop_Increasing Prop_Stable_Suitable Prop_Decreasing
    ## 1               0               0                    0               0
    ## 2               0               0                    0               0
    ## 3               0               0                    0               0
    ##   Prop_Stable_Unsuitable
    ## 1                  38.46
    ## 2                  82.35
    ## 3                  12.50

``` r

# Same aggregation applied to the hypervolume outputs.
hv_regional <- analyze_trends_by_spatial_unit(
  shapefile_path = provinces,
  name_field = "NAME_1",
  binary_stack = hv_consensus$consensus_stack,
  pattern_raster = hv_patterns$pattern,
  time_decrease_raster = hv_patterns$time_decrease,
  time_increase_raster = hv_patterns$time_increase,
  time_steps = time_steps,
  output_dir = "outputs/hv_regional",
  overwrite = TRUE,
  create_plot = TRUE,
  verbose = FALSE
)
```

![](Step3_TemporalModelR_files/figure-html/regional-hv-1.png)![](Step3_TemporalModelR_files/figure-html/regional-hv-2.png)

``` r

# First three rows of the per-province summary for the hypervolume.
head(hv_regional$overall_summary, 3)
```

    ##         Spatial_Unit Always_Absent Always_Present No_Pattern Increasing
    ## 1      Amnat Charoen             0              8          0          0
    ## 2          Ang Thong             0              0          0          0
    ## 3 Bangkok Metropolis             6              0          0          0
    ##   Decreasing Fluctuating Failed Total_Pixels Pct_Always_Absent
    ## 1          0           0     44           52              0.00
    ## 2          0           0     17           17              0.00
    ## 3          0           0     27           33             18.18
    ##   Pct_Always_Present Pct_No_Pattern Pct_Increasing Pct_Decreasing
    ## 1              15.38              0              0              0
    ## 2               0.00              0              0              0
    ## 3               0.00              0              0              0
    ##   Pct_Fluctuating Prop_Increasing Prop_Stable_Suitable Prop_Decreasing
    ## 1               0               0                15.38               0
    ## 2               0               0                 0.00               0
    ## 3               0               0                 0.00               0
    ##   Prop_Stable_Unsuitable
    ## 1                   0.00
    ## 2                   0.00
    ## 3                  18.18

The per-province `overall_summary` tables let you rank provinces by the
magnitude or sign of change and quickly spot regions where the four
models tell different stories.

## Wrap-up

What we covered:

- Read the pre-prepared annual precipitation, temperature, and NDVI
  rasters directly from `data/processed/temporal_split/`.
- Summarised the per-year aligned rasters as a single **temporal-mean**
  3-layer stack (`mean_prec`, `mean_temp`, `mean_NDVI`) and used it as
  the reference E-space for thinning.
- Applied **`bean`** once against that temporal-mean E-space
  (Sheather–Jones bandwidth, stochastic per-pod thinning) to reduce
  environmental sampling bias across the full Sambar deer occurrence
  table before any temporal modelling.
- Built a fully temporally-explicit pipeline: time-stamped occurrences ↔︎
  time-stamped rasters via
  [`temporally_explicit_extraction()`](https://cjhughes926.github.io/TemporalModelR/reference/temporally_explicit_extraction.html).
- Used spatio-temporal cross-validation folds and fold-stratified
  pseudoabsences for an honest evaluation.
- Fitted **all four `TemporalModelR` algorithms** — GLM, GAM, random
  forest, and hypervolume — against the same fold structure and compared
  them on the same held-out test sets.
- Projected every model across 2015–2025 with
  [`generate_spatiotemporal_predictions()`](https://cjhughes926.github.io/TemporalModelR/reference/generate_spatiotemporal_predictions.html)
  and visualised the resulting fold-vote consensus rasters using the
  package’s default
  [`terra::plot()`](https://rspatial.github.io/terra/reference/plot.html)
  styling (3 × 4 grid per model).
- Inspected per-fold E-space and per-time-step G-space performance with
  the package’s default
  [`plot_model_assessment()`](https://cjhughes926.github.io/TemporalModelR/reference/plot_model_assessment.html)
  diagnostic for each of the four models.
- Ran the full post-processing pipeline —
  [`summarize_raster_outputs()`](https://cjhughes926.github.io/TemporalModelR/reference/summarize_raster_outputs.html),
  [`analyze_temporal_patterns()`](https://cjhughes926.github.io/TemporalModelR/reference/analyze_temporal_patterns.html),
  [`analyze_trends_by_spatial_unit()`](https://cjhughes926.github.io/TemporalModelR/reference/analyze_trends_by_spatial_unit.html)
  — **once per model** so consensus frequency rasters, pixel-level
  trajectory classes, and per-province summary tables are directly
  comparable across GLM / GAM / RF / Hypervolume.

## References

- Castaneda-Guzman M, Hughes C, Paansri P, Cobos M (2026). *nicheR:
  Ellipsoid-Based Virtual Niches and Visualization*. R package version
  0.1.0, <https://github.com/castanedaM/nicheR>.
- Chamberlain S, Barve V, Mcglinn D, Oldoni D, Desmet P, Geffert L, Ram
  K (2026). *rgbif: Interface to the Global Biodiversity Information
  Facility API*. R package version 3.8.5.2,
  <https://CRAN.R-project.org/package=rgbif>.
- Hijmans R (2026). *geodata: Access Geographic Data*. R package version
  0.6-10, <https://rspatial.github.io/geodata/>.
- Hughes C, Castaneda-Guzman M, E. Escobar L (2026). *TemporalModelR:
  Temporally Explicit Species Distribution Modelling in R*. R package
  version 0.2.0, <https://github.com/CJHughes926/TemporalModelR>.
- Paansri P, Escobar L (2026). *bean: Data Thinning of Species
  Occurrences in Environmental Space*. R package version 0.2.1,
  <https://github.com/paanwaris/bean>.
