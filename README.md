# Gun Violence and the Political Economy of Forced Labor Investigations in Brazil

**Replication Package (R&R at *World Development*)**  
Gabriel Cepaluni (UNESP) | Jamil Civitarese (New York University)  
Version: October 2025

---

## Overview

This repository reproduces the analyses, figures, and tables for:

> Cepaluni, G., & Civitarese, J. *Gun Violence and the Political Economy of Forced Labor Investigations in Brazil* (R&R at *World Development*).

The project studies how gun availability affects Brazil's anti-slavery enforcement via an informational mechanism, combining a rational-choice framework with a shift–share empirical design.

---

## Repository Structure
```
guns-forced-labor-replication/
│
├── Data/
│   ├── slavery_guns_census_panel.rds
│   ├── slavery_guns_census_panel_with_statecapacity.rds
│   ├── slavery_guns_census_panel_with_statecapacity_cleaned.rds
│   └── slavery_guns_budget_2000_2025.txt
│
├── Code/
│   ├── 01_slaves_guns_main.R
│   ├── 02_slaves_guns_main_replication_tables.R
│   ├── 03_slaves_guns_baseline.R
│   ├── 04_slaves_guns_trends.R
│   ├── 05_slaves_guns_shift_share_initial.R
│   ├── 06_slaves_guns_budget.R
│   ├── 07_slaves_guns_state_capacity.R
│   └── 08_slaves_guns_fiscal_adm_capacity_quartiles.R
│
└── Results/
    ├── tables/
    └── figures/
```

---

## Data Description (core variables)

**Outcomes**
- `operations`: number of anti-forced-labor audits (municipality–year)
- `slave`: number of workers rescued

**Shift–share variables**
- `ss_imports`: small-arms imports × homicide exposure
- `ss_taurus`: Taurus revenue × homicide exposure
- `exposure`: municipality's time-averaged share of national firearm homicides

**Controls**
- `mean_harvested_area_interp`, `gdp_per_capita`, `mean_number_cattle_interp`, `mean_revenue_interp`

**State capacity (when available)**
- `p_idhm` (HDI), `RECORM` (baseline municipal current revenue), and related proxies

**Budget file**
- `slavery_guns_budget_2000_2025.txt`: federal budget series used for execution-ratio plots

> Raw MTE/GEFM microdata are not redistributed. See **Data Access** below.

---

## Execution Order

Run scripts **in numerical order (01–08)**.

### Core analysis
- **01_slaves_guns_main.R**: main regressions (bivariate, +controls, +state×year FE), diagnostic plots (e.g., coefficient plots, parallel-trends, RI), robustness (Conley SEs, Poisson), and heterogeneity by homicide quantiles and regions.
- **02_slaves_guns_main_replication_tables.R**: LaTeX/CSV tables for main and robustness results (incl. leave-one-state-out).

### Additional robustness
- **03_slaves_guns_baseline.R**: baseline-year (2000) characteristics × linear trends robustness; compact tables for `operations` and `slave`.
- **04_slaves_guns_trends.R**: progression of fixed effects (base FE → state×year FE → municipality 5-year grouped trends) and hit-rate checks.

### Instruments & mechanisms
- **05_slaves_guns_shift_share_initial.R**: Bartik instruments using **2–3-period lags** and **baseline-2000 exposure shares**; coefficient plot grid.
- **06_slaves_guns_budget.R**: budget execution ratios (Committed/Liquidated/Paid) and national AFT appointments trend; reads `slavery_guns_budget_2000_2025.txt`.
- **07_slaves_guns_state_capacity.R**: hit-rate LPMs (conditional on audit), state×year FE robustness, and capacity interactions; single summary table.
- **08_slaves_guns_fiscal_adm_capacity_quartiles.R**: heterogeneity by Revenue/GDP and HDI quartiles; RECORM quartiles; exports CSV/LaTeX/Excel tables.

> Outputs are written to `Results/figures/` and `Results/tables/`.

---

## Environment

Tested with **R ≥ 4.4** on Windows. Required packages include:

- `fixest`, `dplyr`, `ggplot2`, `data.table`, `broom`, `tidyr`, `stringi`, `readr`, `scales`
- Tables/exports: `modelsummary`, `xtable`, `writexl`
- Robustness/utilities: `conleyreg`, `ritest`, `cowplot`, `patchwork`, `tidylog`, `haven`

Install (example):
```r
install.packages(c(
  "fixest","dplyr","ggplot2","data.table","broom","tidyr","stringi","readr","scales",
  "modelsummary","xtable","writexl","conleyreg","cowplot","patchwork","tidylog","haven"
))
remotes::install_github("grantmcdermott/ritest")
```

---

## Data Access

Due to **Lei nº 12.527/2011**, raw MTE/GEFM microdata used to construct `operations` and `slave` cannot be redistributed here.

**To request access:**

> Secretaria de Inspeção do Trabalho (SIT) – Ministério do Trabalho e Emprego  
> [https://www.gov.br/pt-br/acesso-a-informacao](https://www.gov.br/pt-br/acesso-a-informacao)

Public sources (e.g., NISAT small-arms transfers, Taurus financial reports, IBGE/IPEA indicators) should be obtained from their official portals.

---

## License

**MIT License** (reuse permitted with attribution).  
© 2025 Gabriel Cepaluni & Jamil Civitarese

---

## Contact

* **Gabriel Cepaluni** — UNESP — [gabriel.cepaluni@unesp.br](mailto:gabriel.cepaluni@unesp.br)
* **Jamil Civitarese** — NYU — [jc9663@nyu.edu](mailto:jc9663@nyu.edu)
