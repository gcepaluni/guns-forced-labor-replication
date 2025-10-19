# Gun Violence and the Political Economy of Forced Labor Investigations in Brazil
**Replication Package (R&R at *World Development*)**

---

## Overview

This repository reproduces all analyses, figures, and tables for:

> Cepaluni, G., & Civitarese, J. *Gun Violence and the Political Economy of Forced Labor Investigations in Brazil* (R&R at *World Development*).

The project examines how gun availability affects Brazil's anti-slavery enforcement efficiency through an informational mechanism, combining a rational-choice framework with a shift–share empirical design.

---

## Repository Structure

```
guns-forced-labor-replication/
│
├── Data/
│   ├── dfGunsPaperNovember24.rds
│   ├── slavery_guns_census_panel.rds
│   ├── slavery_guns_census_panel_with_statecapacity.rds
│   ├── slavery_guns_census_panel_with_statecapacity_cleaned.rds
│   └── slavery_guns_budget_2000_2025.txt
│
├── Code/
│   ├── 01_slaves_guns_main.R
│   ├── 02_slaves_guns_main_replication_tables.R
│   ├── 03_slaves_guns_baseline.R
│   ├── 04_slaves_guns_shift_share_initial.R
│   ├── 05_slaves_guns_budget.R
│   ├── 06_slaves_guns_state_capacity.R
│   └── 07_slaves_guns_fiscal_adm_capacity_quartiles.R
│
└── Results/
    ├── tables/
    └── figures/
```

---

## Code Files

| Script | Description |
|--------|-------------|
| `01_slaves_guns_main.R` | Main analysis: descriptive statistics, fixed-effects and PPML regressions, coefficient plots, randomization inference, Conley SEs, and regional heterogeneity |
| `02_slaves_guns_main_replication_tables.R` | Robustness tables: stepwise controls, leave-one-state-out, and comprehensive robustness checks |
| `03_slaves_guns_baseline.R` | Municipality-specific baseline (2000)×t trends with year and UF×year fixed effects |
| `04_slaves_guns_shift_share_initial.R` | Lagged (two- and three-period) and initial-exposure shift–share instruments |
| `05_slaves_guns_budget.R` | Budget execution (2019–2025) and AFT appointment (2000–2019) figures |
| `06_slaves_guns_state_capacity.R` | Municipal capacity interactions and conditional hit-rate analysis |
| `07_slaves_guns_fiscal_adm_capacity_quartiles.R` | Heterogeneity by fiscal (revenue/GDP) and HDI quartiles |

---

## License

**MIT License**  
© 2025 Gabriel Cepaluni & Jamil Civitarese

---

## Contact

**Gabriel Cepaluni** — UNESP — [gabriel.cepaluni@unesp.br](mailto:gabriel.cepaluni@unesp.br)  
**Jamil Civitarese** — NYU — [jc9663@nyu.edu](mailto:jc9663@nyu.edu)
