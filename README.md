# Model-Assisted Calibration for Integrating Probability-Based Survey Samples

This repository contains reproducibility materials for a manuscript currently
under review. The manuscript proposes model-assisted calibration estimators for
integrating two independently drawn probability-based surveys, where the internal
sample has a variable not observed in the external sample. Analytical Taylor
linearization variance estimators are derived for all proposed methods.

---

## Overview

When a key variable (e.g., a biomarker from DXA body composition scans) is
observed only in a smaller internal survey, larger external surveys that measure
overlapping covariates can be used to improve estimation efficiency. We develop
**GREG (generalized regression) weight adjustment** estimators that reweight the
internal sample so its weighted estimating equations match external information,
without requiring a correctly specified outcome model.

Two settings are considered:

1. **Pooled-sample calibration** — when individual-level microdata from the
   external survey are available, allowing a combined analysis.
2. **External-summary calibration** — when only fitted coefficients and their
   covariance matrix are available from a published external analysis.

Taylor linearization variance estimators account for the complex stratified
multistage cluster designs of both surveys.

---

## Application

We apply the proposed methods to estimate the association between all-cause
mortality (10-year follow-up) and DXA-measured total body fat using:

- **Internal sample**: NHANES 2003–2004 (n₁ = 1,543) — includes DXA body fat
- **External samples**: NHANES 2001–2002 and 2005–2006 (n₂ = 3,107 combined) and
  NHIS 2003–2004 (n₂ = 24,086) — observe BMI, physical activity, smoking, and
  alcohol, but not DXA body fat

---

## Repository Structure

```
survey-calibration/
├── simulation/                   # Simulation reproducibility materials
│   ├── code/
│   │   ├── 01_make_table1_from_canonical_results.R   # Fast route: Table 1
│   │   ├── 02_make_table2_from_canonical_results.R   # Fast route: Table 2
│   │   ├── 03_make_simulation_figure_from_canonical_results.R  # Fast route: Figure 3
│   │   ├── 04_rerun_table1_simulation_optional.R     # Slow route (optional)
│   │   ├── 05_rerun_table2_simulation_optional.R     # Slow route (optional)
│   │   └── helper_functions/
│   │       ├── subfunctions20250403_V1_0.R            # Core functions for Table 1
│   │       ├── subfunctions20250417_V1_3.R            # Core functions for Table 2
│   │       └── seed.txt                               # Monte Carlo seeds (1000 rows)
│   ├── results/
│   │   └── canonical/            # Pre-saved summaries (table1/, table2/, figure/)
│   ├── figures/
│   │   └── Relative Variance.png # Canonical Figure 3
│   └── README.md                 # Detailed instructions for simulation materials
│
├── application/                  # Data application reproducibility materials
│   ├── code/
│   │   ├── 01_prepare_data_application_data.R   # Read raw XPT files → .rds cache
│   │   ├── 02_run_data_application_models.R     # Run calibration models
│   │   ├── 03_make_data_application_figures_tables.R  # Generate figures and tables
│   │   ├── beta_est_var_Biometrics_2.R          # Core variance estimation function
│   │   └── subfunctions20250417_V1_3.R          # Helper functions
│   ├── data/processed/           # Pre-processed analysis-ready data (.rds)
│   ├── figures/primary/          # Manuscript figures (pre-generated)
│   ├── results/                  # Pre-generated tables (CSV)
│   └── README.md                 # Detailed instructions for data application
│
├── README.md                     # This file
├── LICENSE
└── .gitignore
```

---

## Reproducing Results

### Simulation study (fast route — seconds)

Set R's working directory to `simulation/` and run:

```r
source("code/01_make_table1_from_canonical_results.R")   # Table 1
source("code/02_make_table2_from_canonical_results.R")   # Table 2
source("code/03_make_simulation_figure_from_canonical_results.R")  # Figure 3
```

These scripts read pre-saved canonical Monte Carlo summaries and complete in
seconds. See `simulation/README.md` for full details.

### Data application

Raw NHANES and NHIS public-use files must be downloaded from CDC/NCHS and
placed in `application/data/raw/` (see `application/README.md` for file layout).
Then set R's working directory to `application/` and run:

```r
source("code/01_prepare_data_application_data.R")   # ~5–10 min
source("code/02_run_data_application_models.R")     # ~10–20 min
source("code/03_make_data_application_figures_tables.R")
```

### Re-running simulations from scratch (optional; slow)

```r
# Table 1 rerun (several hours)
source("simulation/code/04_rerun_table1_simulation_optional.R")

# Table 2 rerun (many hours; fits glm on N = 1,000,000 per replicate)
source("simulation/code/05_rerun_table2_simulation_optional.R")
```

Rerun outputs go to `simulation/results/rerun_outputs/` and do not overwrite
the canonical results.

---

## Software

All code is written in **R** (≥ 4.1.0). No build system is required.

**Required packages:**

```r
install.packages(c(
  "survey",    # complex survey design analysis
  "MASS",      # matrix operations, ginv
  "numDeriv",  # numerical Jacobians for variance terms
  "ICC",       # intraclass correlation (simulation DGP)
  "haven",     # read NHANES/NHIS XPT files
  "dplyr", "tidyr",
  "ggplot2",
  "writexl", "readxl",
  "ggh4x", "broom", "readr", "stringr"
))
```

---

## Manuscript

This repository contains reproducibility materials for a manuscript currently
under review. Details will be updated upon publication.

---

## Citation

Citation information will be added upon publication.

---

## Authors

Yanhao Lu and Lingxiao Wang  
Department of Statistics, University of Virginia

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
