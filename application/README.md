# Data Application Reproducibility Materials

This folder contains all code to reproduce the data application results reported
in the manuscript, including the confidence interval plot, DXA body fat
diagnostic figure, BMI distribution figure, sample summary table, and CI
estimates table.

---

## Software Requirements

- **R** (>= 4.2.0)
- **R packages**: `haven`, `readr`, `stringr`, `dplyr`, `survey`, `numDeriv`,
  `ggplot2`, `ggh4x`, `broom`, `writexl`, `tidyr`

Install missing packages with:

```r
install.packages(c("haven", "readr", "stringr", "dplyr", "survey", "numDeriv",
                   "ggplot2", "ggh4x", "broom", "writexl", "tidyr"))
```

---

## Data Sources

All raw data files must be placed in `data/raw/` before running the scripts.
The data are publicly available from the U.S. Centers for Disease Control and
Prevention (CDC) and the National Center for Health Statistics (NCHS).

### NHANES (National Health and Nutrition Examination Survey)

Public-use data files for cycles 1999–2000 (cycle A), 2001–2002 (cycle B),
2003–2004 (cycle C), and 2005–2006 (cycle D). Download XPT files from:

  <https://www.cdc.gov/nchs/nhanes/index.htm>

Place files in the following subfolders:

| Subfolder                | Contents                                       |
|--------------------------|------------------------------------------------|
| `data/raw/nhanes/99/`    | NHANES 1999-2000 main exam files               |
| `data/raw/nhanes/99x/`   | NHANES 1999-2000 DXA body composition files    |
| `data/raw/nhanes/99p/`   | NHANES 1999-2000 physical activity files       |
| `data/raw/nhanes/01/`    | NHANES 2001-2002 main exam files               |
| `data/raw/nhanes/01x/`   | NHANES 2001-2002 DXA body composition files    |
| `data/raw/nhanes/01p/`   | NHANES 2001-2002 physical activity files       |
| `data/raw/nhanes/03/`    | NHANES 2003-2004 main exam files               |
| `data/raw/nhanes/03x/`   | NHANES 2003-2004 DXA body composition files    |
| `data/raw/nhanes/03p/`   | NHANES 2003-2004 physical activity files       |
| `data/raw/nhanes/05/`    | NHANES 2005-2006 main exam files               |
| `data/raw/nhanes/05x/`   | NHANES 2005-2006 DXA body composition files    |
| `data/raw/nhanes/05p/`   | NHANES 2005-2006 physical activity files       |

> **Note on 1999–2000:** The 1999–2000 cycle is **excluded from the main
> calibration analysis** due to distributional differences in DXA-measured
> total body fat (see supplementary diagnostic figure). However, Script 01
> reads the 1999–2000 data to build the `nhanes_4cycles.rds` cache used for
> the EDA summary table and the body fat diagnostic figure. Place those XPT
> files in `data/raw/nhanes/99/`, `99x/`, `99p/` and place
> `NHANES_1999_2000_MORT_2019_PUBLIC.dat` in `data/raw/mortality/`.

#### Key variables used

| Domain              | Variables                                                    |
|---------------------|--------------------------------------------------------------|
| Survey design       | `SDMVPSU`, `SDMVSTRA`, `WTINT2YR`, `WTMEC2YR`              |
| Demographics        | `RIDAGEYR`, `RIAGENDR`, `RIDRETH1`, `DMDEDUC2`, `DMDEDUC3` |
| Body measures       | `BMXBMI`, `BMXHT`, `BMXWT`, `WHD010`, `WHD020`             |
| DXA body fat        | `DXDTOFAT`, `DXXTRFAT`, `DXDTRPF`, `DXDTOPF`               |
| Physical activity   | `PADTIMES`, `PADDURAT`, `PADMETS` (PAQIAF files in `*p/`)   |
| Smoking             | `SMQ020`, `SMD030`, `SMQ040`, `SMQ050Q`, `SMQ640`, `SMQ650` |
| Alcohol             | `ALQ101`/`ALD100`, `ALQ110`, `ALQ120Q/U`, `ALQ130`, etc.   |

### NHANES Mortality Linkage

Public-use linked mortality files (NCHS 2019 follow-up). Download from:

  <https://www.cdc.gov/nchs/data-linkage/mortality-public.htm>

Place in `data/raw/mortality/`:

- `NHANES_1999_2000_MORT_2019_PUBLIC.dat`
- `NHANES_2001_2002_MORT_2019_PUBLIC.dat`
- `NHANES_2003_2004_MORT_2019_PUBLIC.dat`
- `NHANES_2005_2006_MORT_2019_PUBLIC.dat`

Fixed-width columns used: `seqn` (1–6), `eligstat` (15), `mortstat` (16),
`permth_exm` (46–48).

### NHIS (National Health Interview Survey)

NHIS 2003–2004 public-use person files. Place a combined two-year CSV at:

  `data/raw/nhis/nhis_2003_2004.csv`

Required columns: `YEAR`, `HHX`, `FMX`, `PX`, `SAMPWEIGHT`, `PSU`, `STRATA`,
`AGE`, `SEX`, `HEIGHT`, `WEIGHT`, `HISPYN`, `RACENEW`, `EDUC`,
`SMOKEV`, `SMOKFREQNOW`, `QUITNO`, `CIGDAYMO`,
`ALC1YR`, `ALCLIFE`, `ALCDAYSYR`, `ALCAMT`,
`MOD10FNO`, `MOD10DNO`, `MOD10FTP`, `MOD10DTP`,
`VIG10FNO`, `VIG10DNO`, `VIG10FTP`, `VIG10DTP`,
`PREGNANTNOW`, `BMICALC`.

### NHIS Mortality Linkage

NCHS 2019 NHIS mortality follow-up files. Place in `data/raw/mortality/`:

- `NHIS_2003_MORT_2019_PUBLIC.dat`
- `NHIS_2004_MORT_2019_PUBLIC.dat`

Fixed-width columns used: `publicid` (1–14), `eligstat` (15), `mortstat` (16),
`dodyear` (23–26).

---

## Analysis Population

- Adults aged 40–69 at the time of NHANES/NHIS interview
- NHANES: complete DXA whole-body scan available (`DXDTOFAT` not missing)
- NHANES: self-reported BMI ≥ 18.5 kg/m²; height 120–220 cm
- NHIS: self-reported BMI ≥ 18.5 kg/m²; pregnant women excluded
- Mortality outcome: all-cause death within 10 years of interview (using
  `permth_exm` ≤ 120 months for NHANES; `dodyear` ≤ 2013/2014 for NHIS)

---

## Running the Code — Required Workflow

Set R's working directory to this folder (`data_application/`) and run scripts
in order:

```r
setwd("path/to/application")

source("code/01_prepare_data_application_data.R")   # ~5-10 min
source("code/02_run_data_application_models.R")      # ~10-20 min
source("code/03_make_data_application_figures_tables.R")
```

Each script must be run from the `data_application/` directory.  Scripts 02
and 03 will stop with an informative error if the required input files are
missing.

---

## Required Outputs (Manuscript and Supplementary Material)

### Primary figures — `figures/primary/`

| File                                                 | Description                                           |
|------------------------------------------------------|-------------------------------------------------------|
| `realdata_ci_plot_all.{pdf,png}`                     | Main CI plot (paper Figure)                           |
| `DXA_bodyfat_by_NHANES_cycle_boxplot.{pdf,png}`      | Supplement: DXA body fat by NHANES cycle and mortality|
| `BMI_distribution_NHANES_NHIS.{pdf,png}`             | Supplement: Weighted BMI distributions                |

> The DXA body fat boxplot (`DXA_bodyfat_by_NHANES_cycle_boxplot.*`)
> corresponds to `Fat.png` from the original analysis.  It shows DXA total
> body fat distributions across all four NHANES cycles (including 1999–2000),
> grouped by 10-year mortality status, motivating the exclusion of 1999–2000
> from the calibration.

### Tables — `results/`

| File                                               | Description                                  |
|----------------------------------------------------|----------------------------------------------|
| `results/model_results_table4.csv`                 | Full results in CSV format                   |
| `results/CI_table_tidy.csv`                        | Tidy CI table in CSV format                  |
| `results/CI_table_wide.csv`                        | CI table wide format                         |
| `results/CI_estimates_table.csv`                   | Formatted CI estimates table (paper Table)   |
| `results/Model_Bias_Final.csv`                     | Multi-model covariate deviation table        |
| `results/eda_diagnostics/sample_summary_table.csv` | Weighted sample summary (Supplementary)      |
| `results/eda_diagnostics/tab_eda_sample_summary.tex` | LaTeX version of the summary table         |

### Intermediate processed data — `data/processed/`

| File                              | Description                                          |
|-----------------------------------|------------------------------------------------------|
| `nhanes_analysis.rds`             | NHANES 3-cycle data (01+03+05); used for calibration |
| `nhis_analysis.rds`               | NHIS 2003-2004 data                                  |
| `nhanes_4cycles.rds`              | All 4 NHANES cycles; used only for EDA/diagnostics   |

These files are analysis-ready caches derived from the publicly available NHANES
public-use XPT files, the NHIS public-use person files, and the NCHS Public-Use
Linked Mortality Files (2019 follow-up). They are included to support a faster
reproducibility workflow (Scripts 02 and 03 can be run directly without
downloading raw data). Raw source files are not redistributed here; if you wish
to rebuild these files from scratch, follow the data download instructions in
the **Data Sources** section above and run `code/01_prepare_data_application_data.R`.

---

## Method Labels

| Code             | Description                                                    |
|------------------|----------------------------------------------------------------|
| Gold             | Benchmark: combined NHANES 2001–2006 (gold standard)          |
| NHANES05         | Internal-only: NHANES 2003–2004 (naive estimator)             |
| calib.NHANES     | NHANES pooled calibration (Chen et al. approach)               |
| M-Based.NHANES   | NHANES external calibration (model-based)                      |
| calib.NHIS       | NHIS pooled calibration                                        |
| M-Based.NHIS     | NHIS external calibration (model-based)                        |

---

## Session Info

To record the R environment used for the analysis:

```r
sessionInfo()
```
