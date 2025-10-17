# Gun Violence and the Political Economy of Forced Labor Investigations in Brazil

**Replication Package for World Development (2025)**  
Gabriel Cepaluni (UNESP & University of Notre Dame)  
Jamil Civitarese (New York University)  
Version: October 2025

---

## Overview

This repository reproduces all analyses, figures, and tables for:

> Cepaluni, G., & Civitarese, J. (2025). *Gun Violence and the Political Economy of Forced Labor Investigations in Brazil.* World Development.

The project examines how gun availability affects the efficiency of Brazil's anti-slavery enforcement, combining a rational-choice model and a shift–share empirical design.

All results, tables, and figures in the paper and Appendix can be generated from the scripts provided here.

---

## Repository Structure
```
guns-forced-labor-replication/
│
├── 01_data/
│   ├── raw/                          # (Not included; instructions below)
│   ├── processed/                    # Cleaned datasets (.rds, .csv)
│   │   ├── slavery_guns_census_panel.rds
│   │   └── slavery_guns_census_panel_with_statecapacity_cleaned.rds
│   └── codebook.pdf                  # Variable descriptions and sources
│
├── 02_code/
│   ├── slaves_guns_shift_share_initial.R           # Lagged & initial exposure instruments
│   ├── slaves_guns_baseline.R                      # Baseline (2000) × t trends
│   ├── slaves_guns_main.R                          # Master analysis pipeline (~1000 lines)
│   ├── slaves_guns_main_replication_tables.R       # Comprehensive LaTeX tables
│   ├── slaves_guns_state_capacity.R                # Hit rates & capacity mechanisms
│   ├── slaves_guns_fiscal_adm_capacity_quartiles.R # Heterogeneity: Revenue/GDP & HDI
│   ├── slaves_guns_state_capacity_quartiles_code.R # Heterogeneity: RECORM quartiles
│   └── slaves_guns_trends.R                        # State×year FE & 5-year muni trends
│
├── 03_output/
│   ├── tables/                       # LaTeX (.tex) and CSV regression tables
│   ├── figures/                      # PDF/PNG coefficient plots and trends
│   └── logs/                         # Diagnostics (RI tests, Conley SEs)
│
├── LICENSE
└── README.md
```

---

## Replication Steps

### Step 1 — Data Preparation

The raw data include both **public** and **restricted** components:

| Source | Access | Notes |
|--------|---------|-------|
| Norwegian Initiative on Small Arms Transfers (NISAT) | Public | Annual small arms imports (USD) |
| Taurus S.A. Financial Reports | Public | Revenue from domestic firearm sales |
| Ministry of Labor and Employment (MTE, GEFM dataset) | Restricted | Forced-labor audits and rescued workers (accessed under Lei nº 12.527/2011) |
| IPEA / IBGE | Public | Socioeconomic indicators and agricultural controls |

**Processed datasets required:**
- `slavery_guns_census_panel.rds` — Main panel (operations, slave, shift-share instruments, controls)
- `slavery_guns_census_panel_with_statecapacity_cleaned.rds` — Extended panel with state capacity variables (HDI, RECORM, security spending)

Users with MTE access should place raw files under `01_data/raw/` and run preprocessing scripts (available upon request).

---

### Step 2 — Execution Order

To reproduce all main and appendix results, run scripts in the following sequence:

#### **Core Results**
1. **`slaves_guns_main.R`** → Master analysis pipeline
   - 12 main regressions (bivariate, controls, state×year FE)
   - Coefficient plots (`main_results.pdf`)
   - Parallel trends tests (`bartik_common_trends_plots.pdf`)
   - Randomization inference (500 reps, `ritest_plots.png`)
   - Rural population robustness (`rural_pop_coef_plot.pdf`)
   - Conley spatial SEs (`shift_share_conley_cowplot.pdf`)
   - Poisson models (`main_results_poisson.pdf`)
   - Homicide quantile heterogeneity (`coefficient_plot_shift_share.png`)
   - Regional analysis (5 Brazilian regions)
   - Budget execution figures (`execution_ratios_trend.png`, `national_trend_appointments.png`)

2. **`slaves_guns_main_replication_tables.R`** → Comprehensive tables for reviewer response
   - **Tables 1–4:** Main results with stepwise controls (LaTeX + CSV + template versions)
   - **Tables 5–8:** Leave-one-state-out robustness (26 states)
   - **Tables 9–12:** Comprehensive robustness (Conley SEs, rural trimming, Poisson)

#### **Robustness Checks**
3. **`slaves_guns_baseline.R`** → Baseline (2000) × t municipality-specific trends
   - Tests 4 baseline characteristics: rural share, log population, log GDP per capita, agricultural intensity
   - Produces 2 compact tables (operations & slave outcomes)

4. **`slaves_guns_trends.R`** → Progressive fixed effects specifications
   - Municipality + Year FE
   - Municipality + State×Year FE
   - Municipality×5-year-group trends (collinearity workaround)
   - **Table 1:** Coefficient stability across specifications
   - **Table 2:** Hit rate analysis with capacity interactions

#### **Mechanism & Heterogeneity**
5. **`slaves_guns_shift_share_initial.R`** → Alternative instruments
   - 2-year lagged exposure
   - 3-year lagged exposure
   - Initial (year 2000) exposure
   - Produces coefficient plot grid (`shift_share_regression_plots_short_titles.pdf`)

6. **`slaves_guns_state_capacity.R`** → Hit rates & false leads tests
   - Linear probability models (LPM) for hit rates (probability audit finds slaves | audit occurred)
   - State×Year FE robustness
   - Capacity interaction effects
   - **Outputs:** `main_results_summary.tex` and `.html`

7. **`slaves_guns_fiscal_adm_capacity_quartiles.R`** → Heterogeneity by municipal characteristics
   - **Part 1:** Revenue/GDP ratio quartiles (fiscal efficiency)
   - **Part 2:** HDI (IDHM) quartiles (human development)
   - Uses `ss_imports` only
   - **Outputs:** CSV, Excel, and LaTeX tables

8. **`slaves_guns_state_capacity_quartiles_code.R`** → RECORM quartiles
   - Municipal current revenue (RECORM) at baseline 2000
   - Tests heterogeneity by absolute fiscal capacity
   - **Outputs:** `recorm_quartile_imports_results.csv/.xlsx/.tex`

**Note:** `slaves_guns_state_capacity_quartiles.R` is an exact duplicate of script #7 and can be ignored.

---

## Output Summary

All outputs are automatically written to `03_output/`:

| Output Type | Files | Description |
|-------------|-------|-------------|
| **Main results** | `tables/table01-04_main_*.tex/csv` | Stepwise controls (audits & slavery × imports/Taurus) |
| **Template tables** | `tables/table01-04_main_*_template.tex` | Same as above, showing all control coefficients |
| **Leave-one-out** | `tables/table05-08_loo_*.tex/csv` | State-by-state robustness (26 states) |
| **Robustness** | `tables/table09-12_robustness_*.tex/csv` | Conley SEs, rural trimming, Poisson |
| **Baseline trends** | `tables/revresp_*_baseline_trends_compact.tex/csv` | Municipality-specific baseline × t |
| **Capacity quartiles** | `tables/*_quartile_*_results.tex/csv/xlsx` | Revenue/GDP, HDI, RECORM heterogeneity |
| **Coefficient plots** | `figures/main_results.pdf` | Grid of 12 main models |
| | `figures/rural_pop_coef_plot.pdf` | Rural population robustness |
| | `figures/shift_share_conley_cowplot.pdf` | Conley spatial SEs |
| | `figures/main_results_poisson.pdf` | Poisson model results |
| **Trends & diagnostics** | `figures/bartik_common_trends_plots.pdf` | Parallel trends by exposure |
| | `figures/coefficient_plot_shift_share*.png` | Homicide quantile heterogeneity |
| | `figures/execution_ratios_trend.png` | Budget execution 2019–2025 |
| | `figures/national_trend_appointments.png` | Labor inspector appointments 2000–2019 |
| **Regional analysis** | `figures/Slavesbartik_*.pdf` | Regional coefficients (slaves outcome) |
| | `figures/Auditbartik_*.pdf` | Regional coefficients (audits outcome) |
| | `figures/sum_of_*_plot.pdf` | Bar charts by region |
| **Randomization inference** | `figures/ritest_plots.png` | RI distributions for 12 models |

---

## Environment and Dependencies

**Replication tested on:**
```
R version 4.4.1 (2024-06-14)
Platform: x86_64-w64-mingw32 / Windows 10

Required packages:
- fixest (v0.12.3)      # Fast fixed-effects models
- dplyr (v1.1.4)        # Data manipulation
- ggplot2 (v3.5.1)      # Visualization
- cowplot (v1.1.3)      # Plot grids
- broom (v1.0.5)        # Model tidying
- data.table (v1.14.8)  # Fast data operations
- xtable (v1.8-4)       # LaTeX tables
- modelsummary (v2.0.0) # Regression tables
- conleyreg (v0.3.2)    # Conley spatial SEs
- ritest (v0.2.0)       # Randomization inference
- writexl (v1.4.2)      # Excel export
- tidylog (v1.0.2)      # Verbose dplyr operations

Optional:
- etwfe (v0.5.0)        # Event-study designs
- haven (v2.5.3)        # Stata file import
- patchwork (v1.1.3)    # Advanced plot composition
```

**Installation:**
```r
install.packages(c("fixest", "dplyr", "ggplot2", "cowplot", "broom", 
                   "data.table", "xtable", "modelsummary", "conleyreg", 
                   "writexl", "tidylog"))
remotes::install_github("grantmcdermott/ritest")
```

---

## Notes on Restricted Data

Due to legal restrictions under **Brazil's Freedom of Information Law (Lei nº 12.527/2011)**, raw data from the Ministry of Labor and Employment (MTE/GEFM) containing municipality-level forced labor audits and worker rescues cannot be redistributed.

**All code is fully functional** with user-supplied datasets containing the same variable structure:
- `ibge_cod` (municipality identifier)
- `year` (2000–2019)
- `operations` (number of labor audits)
- `slave` (number of workers rescued)
- `ss_imports`, `ss_taurus` (shift-share instruments)
- Socioeconomic controls (harvested area, GDP per capita, cattle, municipal revenue)

**To request access:**
> Secretaria de Inspeção do Trabalho (SIT)  
> Ministério do Trabalho e Emprego  
> Portal da Transparência: [https://www.gov.br/pt-br/acesso-a-informacao](https://www.gov.br/pt-br/acesso-a-informacao)

---

## Data and Code Availability Statement

**All replication code and documentation are available at:**  
[https://github.com/gcepaluni/guns-forced-labor-replication](https://github.com/gcepaluni/guns-forced-labor-replication)

**Data sources:**
- **Public datasets** (NISAT, Taurus financials, IBGE/IPEA indicators) are openly available through their respective portals
- **Restricted MTE data** can be obtained upon request from the Ministry of Labor under Lei nº 12.527/2011

---

## License

This replication package is released under the **MIT License**, permitting reuse with attribution.

© 2025 Gabriel Cepaluni & Jamil Civitarese

---

## Citation

If you use these materials, please cite:

> Cepaluni, G., & Civitarese, J. (2025). *Gun Violence and the Political Economy of Forced Labor Investigations in Brazil.* World Development.

---

## Contact

- **Gabriel Cepaluni**  
  Department of International Politics, São Paulo State University (UNESP)  
  Kellogg Institute for International Studies, University of Notre Dame  
  Email: gabriel.cepaluni@unesp.br

- **Jamil Civitarese**  
  The Wilf Family Department of Politics, New York University  
  Email: jc9663@nyu.edu

---

## Troubleshooting

**Common issues:**

1. **"File not found" errors:** Verify paths in each script match your directory structure. Update `in_path` and `out_dir` variables as needed.

2. **Package conflicts:** If `fixest` clustering fails, update to version ≥0.12.0: `install.packages("fixest")`

3. **Memory limitations:** `slaves_guns_main.R` processes ~100,000 municipality-year observations. Recommended RAM: 8GB+

4. **Conley SE warnings:** Spatial autocorrelation calculations require longitude/latitude. Verify these columns exist in your data.

5. **Missing DT_clean/DT_hit objects:** `slaves_guns_trends.R` requires pre-loaded data objects. Run data preparation section from `slaves_guns_main.R` first, or modify script to load data directly.

---
