# Gun Violence and the Political Economy of Forced Labor Investigations in Brazil

**Replication Package for World Development (2025)**  
Gabriel Cepaluni (UNESP)  
Jamil Civitarese (New York University)  
Version: October 2025

---

## Overview

This repository reproduces all analyses, figures, and tables for:

> Cepaluni, G., & Civitarese, J. (2025). *Gun Violence and the Political Economy of Forced Labor Investigations in Brazil.* R&R at World Development.

The project examines how gun availability affects the efficiency of Brazil's anti-slavery enforcement, combining a rational-choice model and a shift–share empirical design.

All results, tables, and figures in the paper and Appendix can be generated from the scripts provided here.

---

## Repository Structure
```
guns-forced-labor-replication/
│
├── 01_data/
│   ├── raw/                          # (Not included; instructions below)
│   └── processed/                    # Cleaned datasets (.rds, .csv)
│       ├── slavery_guns_census_panel.rds
│       ├── slavery_guns_census_panel_with_statecapacity.rds
│       └── slavery_guns_census_panel_with_statecapacity_cleaned.rds
│
├── 02_code/
│   ├── 01_slaves_guns_main.R                          # Master analysis pipeline (~1000 lines)
│   ├── 02_slaves_guns_main_replication_tables.R       # Comprehensive LaTeX tables
│   ├── 03_slaves_guns_baseline.R                      # Baseline (2000) × t trends
│   ├── 04_slaves_guns_trends.R                        # State×year FE & 5-year muni trends
│   ├── 05_slaves_guns_shift_share_initial.R           # Lagged & initial exposure instruments
│   ├── 06_slaves_guns_state_capacity.R                # Hit rates & capacity mechanisms
│   └── 07_slaves_guns_fiscal_adm_capacity_quartiles.R # Heterogeneity: Revenue/GDP, HDI & RECORM
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

#### Available Processed Datasets

The following cleaned datasets are included in `01_data/processed/`:

| File | Size | Last Updated | Description |
|------|------|--------------|-------------|
| `slavery_guns_census_panel.rds` | 28.6 KB | 2025-08-18 | **Main panel dataset** with municipality-year observations (2000-2019): operations, slave, shift-share instruments (ss_imports, ss_taurus), socioeconomic controls, geographic coordinates |
| `slavery_guns_census_panel_with_statecapacity.rds` | 29.8 KB | 2025-09-29 | **Extended panel** adding state capacity variables: HDI (p_idhm), security spending, baseline characteristics for trend analysis |
| `slavery_guns_census_panel_with_statecapacity_cleaned.rds` | 30.2 KB | 2025-09-29 | **Final cleaned version** used in all analyses; includes RECORM (municipal current revenue), complete variable set for heterogeneity analyses |

**Core Variables:**

**Outcomes:**
- `operations` — Number of forced labor audits per municipality-year
- `slave` — Number of workers rescued from slavery-like conditions

**Shift-Share Instruments:**
- `ss_imports` — Gun imports × homicide exposure
- `ss_taurus` — Taurus revenue × homicide exposure
- `exposure` — Municipality's time-averaged share of national homicides

**Controls:**
- `mean_harvested_area_interp` — Agricultural land (hectares)
- `gdp_per_capita` — Municipal GDP per capita (BRL, constant 2010)
- `mean_number_cattle_interp` — Cattle herd size
- `mean_revenue_interp` — Total municipal revenue (BRL)

**State Capacity (extended versions):**
- `p_idhm` — Municipal Human Development Index (0-1)
- `p_desp_funcao_seguranca_mun` — Security spending (% of budget)
- `RECORM` — Municipal current revenue (baseline 2000)

**Geography:**
- `ibge_cod` — 7-digit municipal identifier
- `year` — Calendar year (2000-2019)
- `state_code` — State identifier (UF code)
- `latitude`, `longitude` — Municipal centroids

Users with MTE access should place raw files under `01_data/raw/` and run preprocessing scripts (available upon request).

---

### Step 2 — Execution Order

To reproduce all main and appendix results, **run scripts 01–07 in numerical order:**

#### **Core Results**

**1. `01_slaves_guns_main.R`** → Master analysis pipeline
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

**2. `02_slaves_guns_main_replication_tables.R`** → Comprehensive tables for reviewer response
   - **Tables 1–4:** Main results with stepwise controls (LaTeX + CSV + template versions)
   - **Tables 5–8:** Leave-one-state-out robustness (26 states)
   - **Tables 9–12:** Comprehensive robustness (Conley SEs, rural trimming, Poisson)

#### **Robustness Checks**

**3. `03_slaves_guns_baseline.R`** → Baseline (2000) × t municipality-specific trends
   - Tests 4 baseline characteristics: rural share, log population, log GDP per capita, agricultural intensity
   - Produces 2 compact tables (operations & slave outcomes)

**4. `04_slaves_guns_trends.R`** → Progressive fixed effects specifications
   - Municipality + Year FE
   - Municipality + State×Year FE
   - Municipality×5-year-group trends (collinearity workaround)
   - **Table 1:** Coefficient stability across specifications
   - **Table 2:** Hit rate analysis with capacity interactions

#### **Mechanism & Heterogeneity**

**5. `05_slaves_guns_shift_share_initial.R`** → Alternative instruments
   - 2-year lagged exposure
   - 3-year lagged exposure
   - Initial (year 2000) exposure
   - Produces coefficient plot grid (`shift_share_regression_plots_short_titles.pdf`)

**6. `06_slaves_guns_state_capacity.R`** → Hit rates & false leads tests
   - Linear probability models (LPM) for hit rates (probability audit finds slaves | audit occurred)
   - State×Year FE robustness
   - Capacity interaction effects
   - **Outputs:** `main_results_summary.tex` and `.html`

**7. `07_slaves_guns_fiscal_adm_capacity_quartiles.R`** → Heterogeneity by municipal characteristics
   - **Part 1:** Revenue/GDP ratio quartiles (fiscal efficiency)
   - **Part 2:** HDI (IDHM) quartiles (human development)
   - **Part 3:** RECORM quartiles (absolute fiscal capacity)
   - Uses `ss_imports` only
   - **Outputs:** CSV, Excel, and LaTeX tables for all three heterogeneity analyses

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

**All code is fully functional** with user-supplied datasets containing the same variable structure.

**Minimum required variables:**
- `ibge_cod` (municipality identifier)
- `year` (2000–2019)
- `operations` (number of labor audits)
- `slave` (number of workers rescued)
- `ss_imports`, `ss_taurus` (shift-share instruments)
- Socioeconomic controls (harvested area, GDP per capita, cattle, municipal revenue)

**To request access to MTE data:**

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
- **Processed datasets** included in this repository contain municipality-year panels (2000-2019) with forced labor outcomes, shift-share instruments, and socioeconomic controls

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
  Email: gabriel.cepaluni@unesp.br

- **Jamil Civitarese**  
  The Wilf Family Department of Politics, New York University  
  Email: jc9663@nyu.edu

---
```

