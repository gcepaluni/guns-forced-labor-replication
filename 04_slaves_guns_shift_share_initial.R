################################################################################
# GUN VIOLENCE AND THE POLITICAL ECONOMY OF FORCED LABOR INVESTIGATIONS IN BRAZIL
# TWO- AND THREE-PERIOD LAGGED SHIFT–SHARE ANALYSIS WITH INITIAL SHARES
# Author: Gabriel Cepaluni
# Date:   2025-05-16
################################################################################
# DESCRIPTION:
#   Implements lagged and initial exposure analyses for robustness testing of
#   “Gun Violence and the Political Economy of Forced Labor Investigations in Brazil.”
#   Builds shift–share instruments with 2- and 3-period lags and baseline (initial)
#   exposure shares, producing regression estimates and coefficient plots.
#
# OUTPUTS:
#   • shift_share_regression_plots_short_titles.pdf
#   • All estimates saved under:
#       C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Results
#
# NOTES:
#   • Fully compatible with main replication pipeline
#   • Uses ASCII-safe titles for cross-platform LaTeX/HTML export
#   • Ensures consistent lag logic and exposure construction across specifications
################################################################################

#### 0) CLEAN WORKSPACE AND SETUP ###############################################
rm(list = ls())

setwd("C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Results")
if (!dir.exists(getwd())) stop("Working directory does not exist: ", getwd())
cat("Working directory is now:", getwd(), "/n")

#### 1) PACKAGE INSTALLATION & LOADING ##########################################
required_packages <- c(
  "haven", "fixest", "dplyr", "remotes", "ritest", "ggplot2", "broom",
  "patchwork", "tidyverse", "Hmisc", "psych", "xtable", "broom.mixed",
  "coefplot", "gridExtra", "purrr", "conleyreg", "etwfe", "cowplot", "zoo"
)
new_packages <- setdiff(required_packages, installed.packages()[, "Package"])
if (length(new_packages)) install.packages(new_packages)
if (!requireNamespace("ritest", quietly = TRUE)) remotes::install_github("grantmcdermott/ritest")
invisible(lapply(required_packages, library, character.only = TRUE))

#### 2) LOAD RAW DATA ###########################################################
df_reg <- readRDS("C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Data/dfGunsPaperNovember24.rds") %>%
  select(
    ibge_cod, year, homicidios, suicidios, slave, operations,
    taurusRevenue, values_imp, mean_uf_cod,
    mean_harvested_area_interp, mean_pib_mun_interp, mean_population_interp,
    mean_number_cattle_interp, mean_revenue_interp,
    pop_rural_2000, latitude, longitude
  ) %>%
  filter(mean_revenue_interp > 1000) %>%
  na.omit()

#### 3) LEVEL SHARES & SHIFT–SHARE ##############################################
df_clean <- df_reg %>%
  mutate(mean_uf_cod = as.numeric(substr(ibge_cod, 1, 2))) %>%
  filter(mean_uf_cod != 28, year >= 2000) %>%
  group_by(year) %>%
  mutate(
    rel_hom = homicidios / sum(homicidios, na.rm = TRUE),
    rel_sui = suicidios / sum(suicidios, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  group_by(ibge_cod) %>%
  mutate(
    exposure             = mean(rel_hom, na.rm = TRUE),
    exp_suicide          = mean(rel_sui, na.rm = TRUE),
    initial_exposure     = rel_hom[year == min(year)],
    initial_suicide      = rel_sui[year == min(year)],
    imports_scaled       = values_imp / max(values_imp, na.rm = TRUE),
    taurusRevenue_scaled = taurusRevenue / max(taurusRevenue, na.rm = TRUE),
    ss_imports           = exposure * imports_scaled,
    ss_taurus            = exposure * taurusRevenue_scaled,
    ss_imports_suic      = exp_suicide * imports_scaled,
    ss_taurus_suic       = exp_suicide * taurusRevenue_scaled,
    gdp_per_capita       = mean_pib_mun_interp / mean_population_interp
  ) %>%
  ungroup()

#### 4) CREATE 2-3 PERIOD LAGS & INITIAL SHARES #################################
df_clean <- df_clean %>%
  arrange(ibge_cod, year) %>%
  group_by(ibge_cod) %>%
  mutate(
    lag2_exposure        = lag(exposure, 2),
    lag3_exposure        = lag(exposure, 3),
    lag2_exp_suicide     = lag(exp_suicide, 2),
    lag3_exp_suicide     = lag(exp_suicide, 3),

    lag2_ss_imports      = lag2_exposure * imports_scaled,
    lag3_ss_imports      = lag3_exposure * imports_scaled,
    lag2_ss_taurus       = lag2_exposure * taurusRevenue_scaled,
    lag3_ss_taurus       = lag3_exposure * taurusRevenue_scaled,

    lag2_ss_imports_suic = lag2_exp_suicide * imports_scaled,
    lag3_ss_imports_suic = lag3_exp_suicide * imports_scaled,
    lag2_ss_taurus_suic  = lag2_exp_suicide * taurusRevenue_scaled,
    lag3_ss_taurus_suic  = lag3_exp_suicide * taurusRevenue_scaled,

    ss_imports_init      = initial_exposure * imports_scaled,
    ss_taurus_init       = initial_exposure * taurusRevenue_scaled
  ) %>%
  ungroup()

#### MAIN REGRESSIONS ###########################################################
controls <- c(
  "log(mean_harvested_area_interp + 1)",
  "log(gdp_per_capita + 1)",
  "log(mean_number_cattle_interp + 1)",
  "log(mean_revenue_interp + 1)"
)

formulas <- list(
  reg_lag2_imp  = paste("operations ~ lag2_ss_imports +", paste(controls, collapse = " + "), "| ibge_cod + year"),
  reg_lag3_imp  = paste("operations ~ lag3_ss_imports +", paste(controls, collapse = " + "), "| ibge_cod + year"),
  reg_init_imp  = paste("operations ~ ss_imports_init +", paste(controls, collapse = " + "), "| ibge_cod + year"),
  reg_lag2_tau  = paste("operations ~ lag2_ss_taurus +", paste(controls, collapse = " + "), "| ibge_cod + year"),
  reg_lag3_tau  = paste("operations ~ lag3_ss_taurus +", paste(controls, collapse = " + "), "| ibge_cod + year"),
  reg_init_tau  = paste("operations ~ ss_taurus_init +", paste(controls, collapse = " + "), "| ibge_cod + year")
)

models <- lapply(formulas, function(f) {
  feols(as.formula(f), data = df_clean, cluster = ~ibge_cod)
})

#### COEFFICIENT PLOTS ##########################################################
plot_coef <- function(m, title) {
  df <- broom::tidy(m) %>%
    filter(grepl("ss_", term)) %>%
    mutate(term = case_when(
      grepl("lag2_ss_imports", term) ~ "Lag2 Imports",
      grepl("lag3_ss_imports", term) ~ "Lag3 Imports",
      grepl("ss_imports_init", term) ~ "Init Imports",
      grepl("lag2_ss_taurus", term) ~ "Lag2 Taurus",
      grepl("lag3_ss_taurus", term) ~ "Lag3 Taurus",
      grepl("ss_taurus_init", term) ~ "Init Taurus",
      TRUE ~ term
    ))
  ggplot(df, aes(term, estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = title, y = "Estimate", x = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

titles <- c(
  "Lag2 Imports",
  "Lag3 Imports",
  "Init Imports",
  "Lag2 Taurus",
  "Lag3 Taurus",
  "Init Taurus"
)

plots <- mapply(plot_coef, models, titles, SIMPLIFY = FALSE)
combined <- cowplot::plot_grid(plotlist = plots, ncol = 3)

ggsave("shift_share_regression_plots_short_titles.pdf", combined, width = 14, height = 8)
