<p align="center">
  <img src="man/figures/logo.png" alt="ENM-Thailand workshop 2026 hex sticker" width="240"/>
</p>

<h1 align="center">ENM-Thailand workshop 2026</h1>

<p align="center">
  <a href="https://paanwaris.github.io/ENM-Thailand/">
    <img src="https://img.shields.io/badge/pkgdown-site-1B6E3F?logo=r&logoColor=white" alt="pkgdown site"/>
  </a>
  <a href="https://github.com/paanwaris/ENM-Thailand/blob/main/LICENSE.md">
    <img src="https://img.shields.io/badge/license-MIT-F2B705" alt="MIT license"/>
  </a>
</p>

A three-day hands-on workshop on ecological niche modelling (ENM) and species distribution modelling (SDM), using Sambar deer (*Rusa unicolor*) in Thailand as the case study.

| Day | Package | Repository | Focus |
|-----|---------|-----------|-------|
| 1 | **nicheR** | <https://github.com/castanedaM/nicheR> | Ellipsoid-based virtual niches (virtual species — no occurrence data needed) |
| 2a | — (data prep) | — | Download GBIF + WorldClim + GADM into `data/processed/` |
| 2b | **bean**   | <https://github.com/paanwaris/bean> | Reducing environmental sampling bias by thinning in E-space |
| 3 | **TemporalModelR** | <https://github.com/CJHughes926/TemporalModelR> | Temporally explicit SDMs using local annual rasters (2010–2025) |

The full rendered website lives at **<https://paanwaris.github.io/ENM-Thailand/>**.

---

## Authors

| | Role | ORCID |
|---|---|---|
| **Paanwaris Paansri** | Main author / Maintainer · `paanwaris@vt.edu` | <https://orcid.org/0000-0001-9992-098X> |
| **Luis E. Escobar** | Co-author | <https://orcid.org/0000-0001-5735-2750> |

---

## Repository structure

```
ENM-Thailand/
├── README.md
├── DESCRIPTION                         ← package metadata + authorship
├── LICENSE / LICENSE.md                ← MIT
├── _pkgdown.yml                        ← pkgdown site config
├── .gitignore
├── .Renviron.example                   ← template for GBIF credentials
├── ENM-Thailand.Rproj                  ← create on first open in RStudio
├── Day1_nicheR.Rmd                     ← canonical Day 1 lesson
├── Day2a_Data_Download.Rmd             ← canonical Day 2a — downloads
├── Day2b_Bean_Processing.Rmd           ← canonical Day 2b lesson
├── Day3_TemporalModelR.Rmd             ← canonical Day 3 lesson
├── vignettes/                          ← thin wrappers consumed by pkgdown
│   ├── Day1_nicheR.Rmd
│   ├── Day2a_Data_Download.Rmd
│   ├── Day2b_Bean_Processing.Rmd
│   └── Day3_TemporalModelR.Rmd
├── scripts/
│   ├── build_pkgdown_site.R            ← builds the public website
│   └── make_hex_logo.R                 ← regenerates man/figures/logo.png
├── man/figures/
│   └── logo.png                        ← workshop hex logo (referenced in README)
├── data/
│   ├── raw/                            ← raw downloads (gitignored)
│   └── processed/                      ← cleaned outputs (gitignored)
├── temporal_rasters/                   ← annual LST + Precip rasters 2010–2025 (provided)
├── outputs/                            ← model objects, prediction rasters
└── figures/                            ← rendered plots
```

Package installation is embedded in each Rmd — there is no central install script. The first chunk of every Rmd installs whichever workshop package that day needs.

The `vignettes/` folder contains thin wrappers (`child = "../DayX_*.Rmd"`) so pkgdown renders the **same source** as the canonical top-level files — no duplication.

---

## Quick start

### 1. Clone the repo

```bash
git clone https://github.com/paanwaris/ENM-Thailand.git
cd ENM-Thailand
```

### 2. Set your GBIF credentials (one time)

Day 2a reads GBIF credentials from environment variables so nothing confidential ever lands in git. Copy the template and fill in your details:

```bash
cp .Renviron.example .Renviron
```

```r
usethis::edit_r_environ()   # opens .Renviron in RStudio
```

```
GBIF_USER  = "paanwaris"
GBIF_EMAIL = "paanwaris@vt.edu"
GBIF_PWD   = "********"
```

Restart R after editing so the variables are loaded. `.Renviron` is in `.gitignore`.

### 3. Knit Day 2a first (≈15 min)

Day 2a populates `data/processed/` with everything the other days read. Run it once:

```r
rmarkdown::render("Day2a_Data_Download.Rmd")
```

### 4. Knit each lesson day

Each Rmd is self-contained and can be knit independently after Day 2a has been run once:

```r
rmarkdown::render("Day1_nicheR.Rmd")
rmarkdown::render("Day2b_Bean_Processing.Rmd")
rmarkdown::render("Day3_TemporalModelR.Rmd")
```

Equivalent from the shell:

```bash
Rscript -e 'rmarkdown::render("Day1_nicheR.Rmd")'
Rscript -e 'rmarkdown::render("Day2b_Bean_Processing.Rmd")'
Rscript -e 'rmarkdown::render("Day3_TemporalModelR.Rmd")'
```

---

## Building the pkgdown website

The public site lives at `https://paanwaris.github.io/ENM-Thailand/` and is regenerated from the same Rmds via the wrappers in `vignettes/`.

### One-time setup

```r
install.packages(c("pkgdown", "rmarkdown", "knitr"))
usethis::use_pkgdown_github_pages()   # optional: wires up the gh-pages workflow
```

### Build the site locally

```r
source("scripts/build_pkgdown_site.R")
```

This calls `pkgdown::build_site()` and writes the rendered HTML into `docs/` (gitignored). The four workshop Rmds appear under the **Workshop days** menu in the navbar via the `articles:` section of `_pkgdown.yml`.

### Re-render a single article

```r
pkgdown::build_article("Day1_nicheR")
pkgdown::build_article("Day2a_Data_Download")
pkgdown::build_article("Day2b_Bean_Processing")
pkgdown::build_article("Day3_TemporalModelR")
```

---

## Regenerating the hex logo

The current `man/figures/logo.png` ships with the repo. If you want to
re-render it (e.g. swap the year on the next workshop), run:

```r
source("scripts/make_hex_logo.R")
```

This uses the [`hexSticker`](https://github.com/GuangchuangYu/hexSticker)
package to render a fresh `man/figures/logo.png`.

---

## Data sources & citations

- **GBIF**: GBIF.org occurrence download for *Rusa unicolor* (downloaded fresh each run via the `rgbif` package). Please cite the DOI returned by your download.
- **WorldClim v2.1**: Fick, S.E. and R.J. Hijmans, 2017. WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. *International Journal of Climatology* 37 (12): 4302–4315. <https://www.worldclim.org/>
- **Thailand boundary**: GADM v4.1 (free for academic and non-commercial use).
- **Workshop packages**:
  - `nicheR` — Castaneda-Guzman, Hughes, Paansri, Cobos
  - `bean` — Paansri & Escobar
  - `TemporalModelR` — Hughes, Castaneda-Guzman, Escobar (2026)

---

## License

MIT — see [`LICENSE.md`](LICENSE.md). Each underlying R package retains its own license.

## Contributing

Issues and pull requests welcome. If a code chunk does not run on your system, please open an Issue including your `sessionInfo()` output.
