# ENM-Thailand

**A three-day R Markdown workshop on ecological niche modeling (ENM) and species distribution modeling (SDM) in Thailand.**

This repository contains the materials for a hands-on workshop that walks through three modern R packages for ENM/SDM, using the Sambar deer (*Rusa unicolor*) as a case study species and Thailand as the study region.

| Day | Package | Repository | Focus |
|-----|---------|-----------|-------|
| 1 | **nicheR** | <https://github.com/castanedaM/nicheR> | Ellipsoid-based virtual niches (virtual species — no occurrence data needed) |
| 2 | **bean**   | <https://github.com/paanwaris/bean> | Reducing environmental sampling bias by thinning in E-space |
| 3 | **TemporalModelR** | <https://github.com/CJHughes926/TemporalModelR> | Temporally explicit SDMs |

---

## Repository structure

```
ENM-Thailand/
├── README.md                   ← this file
├── .gitignore
├── ENM-Thailand.Rproj          ← create on first open in RStudio
├── scripts/
│   ├── 00_install_packages.R   ← installs every R dependency
│   └── 00_download_data.R      ← downloads GBIF + WorldClim, crops to Thailand
├── Day1_nicheR.Rmd             ← Day 1 lesson
├── Day2_bean.Rmd               ← Day 2 lesson
├── Day3_TemporalModelR.Rmd     ← Day 3 lesson
├── data/
│   ├── raw/                    ← raw downloads (kept out of git)
│   └── processed/              ← cleaned GBIF + cropped bioclim (kept out of git)
├── outputs/                    ← model objects, prediction rasters
└── figures/                    ← rendered plots
```

The `data/`, `outputs/`, and `figures/` folders are populated by the Day-0 download script and by knitting the three lessons.

---

## How to run the workshop

### 1. Clone the repo

```bash
git clone https://github.com/<your-username>/ENM-Thailand.git
cd ENM-Thailand
```

Or use **GitHub Desktop** → *File → Clone repository*.

### 2. Open the R project

Open `ENM-Thailand.Rproj` in RStudio. If the `.Rproj` file does not exist yet, RStudio will create one for you (*File → New Project → Existing Directory*).

### 3. Install dependencies (one time)

```r
source("scripts/00_install_packages.R")
```

This installs CRAN dependencies plus the three workshop packages from GitHub.

### 4. Download the data (one time, ~15 min)

```r
source("scripts/00_download_data.R")
```

This downloads:

- GBIF occurrence records for **Sambar deer** (*Rusa unicolor*) cropped to Thailand.
- **WorldClim v2.1** 19 bioclimatic variables at **10 arc-minute** resolution, cropped to Thailand (WGS84 / EPSG:4326).

Results land in `data/processed/`. The Day 2 and Day 3 Rmd files load directly from this folder, so you only need to do this once.

### 5. Knit each day

In RStudio, open each `.Rmd` and click **Knit**, or:

```r
rmarkdown::render("Day1_nicheR.Rmd")
rmarkdown::render("Day2_bean.Rmd")
rmarkdown::render("Day3_TemporalModelR.Rmd")
```

---

## Workshop content

### Day 1 — `Day1_nicheR.Rmd`

Build a **virtual species** using `nicheR`. We define an ellipsoid in environmental space, project it onto Thailand bioclim rasters, sample virtual occurrences, and simulate virtual communities. No empirical occurrence data are used on Day 1 — the goal is to understand ENM mechanics by working with a niche we control.

### Day 2 — `Day2_bean.Rmd`

Use `bean` to reduce environmental sampling bias in **real** Sambar deer occurrence data. We prepare GBIF records, find an objective grid resolution in E-space, thin clusters, fit an ellipsoid niche, and project suitability back to Thailand.

### Day 3 — `Day3_TemporalModelR.Rmd`

Use `TemporalModelR` to build a **temporally explicit SDM** for Sambar deer. We align rasters, perform spatiotemporal rarefaction, generate fold-stratified pseudoabsences, fit a GLM and a random forest, project predictions through time, and classify pixel-level temporal trajectories.

---

## Data sources & citations

- **GBIF**: GBIF.org occurrence download for *Rusa unicolor* (downloaded fresh each run via the `rgbif` package). Please cite the DOI returned by your download.
- **WorldClim v2.1**: Fick, S.E. and R.J. Hijmans, 2017. WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. *International Journal of Climatology* 37 (12): 4302–4315. <https://www.worldclim.org/>
- **Thailand boundary**: GADM v4.1 (free for academic and non-commercial use).
- **Packages**:
  - `nicheR` — Castaneda-Guzman, Hughes, Paansri, Cobos
  - `bean` — Paansri & Escobar
  - `TemporalModelR` — Hughes, Castaneda-Guzman, Escobar (2026)

---

## License

MIT for the workshop materials. Each underlying R package retains its own license (all MIT at the time of writing).

## Contributing

Issues and pull requests welcome. If a code chunk does not run on your system, please open an Issue including your `sessionInfo()` output.
