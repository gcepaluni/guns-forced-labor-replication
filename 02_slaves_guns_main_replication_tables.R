################################################################################
# GUN VIOLENCE AND THE POLITICAL ECONOMY OF FORCED LABOR INVESTIGATIONS IN BRAZIL
# COMPREHENSIVE ROBUSTNESS PIPELINE — REPLICATION FOR REVIEWER RESPONSE
# Author: Gabriel Cepaluni
# Date:   2025
################################################################################
# DESCRIPTION:
#   Generates all robustness outputs for
#   “Gun Violence and the Political Economy of Forced Labor Investigations in Brazil.”
#   Specifically builds:
#     • Main and template regression tables
#     • Leave-one-out (LOO) robustness checks
#     • Extended and reviewer-requested robustness tables
#
# OUTPUTS:
#   All LaTeX (.tex) and CSV (.csv) tables are saved under:
#     C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Results
#
# NOTES:
#   • Fully reproducible with R ≥ 4.3 and fixest ≥ 0.11
#   • Ensures consistency with the main analysis pipeline
#   • Uses ASCII-safe export for journal submission
################################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

# Load required libraries
suppressPackageStartupMessages({
  library(haven)
  library(fixest)
  library(dplyr)
  library(xtable)
  library(broom)
  library(tidyr)
  library(purrr)
})

#####################################################################################################################################
######################################## OUTPUT DIRECTORY & HELPERS #################################################################
#####################################################################################################################################

# All outputs will be written here (Windows-safe forward slashes)
out_dir <- "C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Track only the files produced in THIS run
created_files <- character(0)

save_tex <- function(xt, name) {
  path <- file.path(out_dir, name)
  print(xt,
        file = path,
        floating = FALSE, include.rownames = FALSE, type = "latex")
  created_files <<- c(created_files, path)
}
save_csv <- function(df, name) {
  path <- file.path(out_dir, name)
  write.csv(df, path, row.names = FALSE)
  created_files <<- c(created_files, path)
}

# R² (within) safely
safe_r2_within <- function(model) {
  tryCatch(as.numeric(fixest::fitstat(model, "r2_within")), error = function(e) NA_real_)
}

# Significance stars
stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.010) return("**")
  if (p < 0.050) return("*")
  if (p < 0.100) return("+")
  ""
}

#####################################################################################################################################
######################################## DATA PREPARATION ###########################################################################
#####################################################################################################################################

# Load and prepare the dataset
df_reg <- readRDS("C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Data/dfGunsPaperNovember24.rds") %>%
  select(
    ibge_cod, year, homicidios, suicidios, slave, operations, taurusRevenue, values_imp,
    mean_uf_cod, mean_harvested_area_interp, mean_pib_mun_interp, mean_population_interp,
    mean_number_cattle_interp, mean_revenue_interp, pop_rural_2000, latitude, longitude
  ) %>%
  filter(mean_revenue_interp > 1000) %>%
  tidyr::drop_na()

df_clean <- df_reg %>%
  mutate(
    ibge_cod   = as.character(ibge_cod),
    state_code = as.numeric(substr(ibge_cod, 1, 2))
  ) %>%
  filter(state_code != 28, year >= 2000) %>%
  group_by(year) %>%
  mutate(
    rel_hom = homicidios / sum(homicidios, na.rm = TRUE),
    rel_sui = suicidios  / sum(suicidios,  na.rm = TRUE)
  ) %>%
  ungroup() %>%
  group_by(ibge_cod) %>%
  mutate(
    exposure             = mean(rel_hom, na.rm = TRUE),
    exp_suicide          = mean(rel_sui,  na.rm = TRUE),
    imports_scaled       = values_imp       / max(values_imp,       na.rm = TRUE),
    taurusRevenue_scaled = taurusRevenue    / max(taurusRevenue,    na.rm = TRUE),
    ss_imports           = exposure    * imports_scaled,
    ss_taurus            = exposure    * taurusRevenue_scaled,
    ss_imports_suic      = exp_suicide * imports_scaled,
    ss_taurus_suic       = exp_suicide * taurusRevenue_scaled,
    gdp_per_capita       = mean_pib_mun_interp / mean_population_interp
  ) %>%
  ungroup()

cat("Dataset loaded:", nrow(df_clean), "observations,", n_distinct(df_clean$ibge_cod), "municipalities\n")
cat("Time period:", min(df_clean$year), "to", max(df_clean$year), "\n")

#####################################################################################################################################
######################################## STATS EXTRACTION ###########################################################################
#####################################################################################################################################

extract_reg_stats <- function(reg_model, outcome_var, data, instrument_var) {
  if (is.null(reg_model)) return(NULL)
  cs <- broom::tidy(reg_model)
  mc <- cs[cs$term == instrument_var, , drop = FALSE]
  if (nrow(mc) == 0) return(NULL)
  
  coef_est  <- mc$estimate[1]
  std_error <- mc$std.error[1]
  t_stat    <- coef_est / std_error
  p_val     <- mc$p.value[1]
  
  data.frame(
    coefficient  = round(coef_est, 4),
    std_error    = round(std_error, 4),
    t_statistic  = round(t_stat, 3),
    p_value      = round(p_val, 4),
    r_squared    = round(safe_r2_within(reg_model), 4),
    n_obs        = reg_model$nobs,
    n_clusters   = dplyr::n_distinct(data$ibge_cod),
    outcome_mean = round(mean(data[[outcome_var]], na.rm = TRUE), 3)
  )
}

#####################################################################################################################################
######################################## MAIN RESULTS — STEPWISE CONTROLS ###########################################################
#####################################################################################################################################

create_main_results_table <- function(outcome_var, instrument_var, data, table_title) {
  
  # (1) Municipality & Year FE
  f1 <- as.formula(paste(outcome_var, "~", instrument_var, "| ibge_cod + year"))
  m1 <- feols(f1, panel.id = ~ ibge_cod + year, cluster = ~ ibge_cod, data = data)
  
  # (2) + Economic controls
  f2 <- as.formula(paste0(
    outcome_var, " ~ ", instrument_var, " + ",
    "log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + ",
    "log(mean_number_cattle_interp+1) + log(mean_revenue_interp+1) | ibge_cod + year"
  ))
  m2 <- feols(f2, panel.id = ~ ibge_cod + year, cluster = ~ ibge_cod, data = data)
  
  # (3) + State linear trend
  f3 <- as.formula(paste0(
    outcome_var, " ~ ", instrument_var, " + ",
    "log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + ",
    "log(mean_number_cattle_interp+1) + log(mean_revenue_interp+1) | ",
    "ibge_cod + year + state_code[year]"
  ))
  m3 <- feols(f3, panel.id = ~ ibge_cod + year, cluster = ~ ibge_cod, data = data)
  
  s1 <- extract_reg_stats(m1, outcome_var, data, instrument_var)
  s2 <- extract_reg_stats(m2, outcome_var, data, instrument_var)
  s3 <- extract_reg_stats(m3, outcome_var, data, instrument_var)
  
  tab <- data.frame(
    Specification = c("Municipality + Year FE", "Previous + Controls", "Previous + State linear trend"),
    Coefficient   = c(s1$coefficient,  s2$coefficient,  s3$coefficient),
    Std_Error     = c(s1$std_error,    s2$std_error,    s3$std_error),
    T_Statistic   = c(s1$t_statistic,  s2$t_statistic,  s3$t_statistic),
    P_Value       = c(s1$p_value,      s2$p_value,      s3$p_value),
    R_Squared     = c(s1$r_squared,    s2$r_squared,    s3$r_squared),
    N_Obs         = c(s1$n_obs,        s2$n_obs,        s3$n_obs),
    N_Clusters    = c(s1$n_clusters,   s2$n_clusters,   s3$n_clusters),
    Outcome_Mean  = c(s1$outcome_mean, s2$outcome_mean, s3$outcome_mean)
  )
  
  latex <- xtable(
    tab,
    caption = table_title,
    label   = paste0("tab:", gsub("[^A-Za-z0-9]", "_", tolower(table_title))),
    digits  = c(0, 0, 4, 4, 3, 4, 4, 0, 0, 3)
  )
  list(table = tab, latex = latex, models = list(m1, m2, m3))
}

cat("\n================================================================================\n")
cat("MAIN RESULTS TABLES — STEPWISE ADDITION OF CONTROLS\n")
cat("================================================================================\n\n")

t_ops_imp  <- create_main_results_table("operations", "ss_imports", df_clean,
                                        "Effect of Gun Imports Shift-Share on Labor Audits")
t_ops_tau  <- create_main_results_table("operations", "ss_taurus",  df_clean,
                                        "Effect of Taurus Revenue Shift-Share on Labor Audits")
t_slv_imp  <- create_main_results_table("slave",      "ss_imports", df_clean,
                                        "Effect of Gun Imports Shift-Share on Slavery Cases")
t_slv_tau  <- create_main_results_table("slave",      "ss_taurus",  df_clean,
                                        "Effect of Taurus Revenue Shift-Share on Slavery Cases")

# Save Tables 1–4 (LaTeX + CSV)
save_tex(t_ops_imp$latex,  "table01_main_audits_imports.tex");  save_csv(t_ops_imp$table,  "table01_main_audits_imports.csv")
save_tex(t_ops_tau$latex,  "table02_main_audits_taurus.tex");   save_csv(t_ops_tau$table,  "table02_main_audits_taurus.csv")
save_tex(t_slv_imp$latex,  "table03_main_slavery_imports.tex"); save_csv(t_slv_imp$table,  "table03_main_slavery_imports.csv")
save_tex(t_slv_tau$latex,  "table04_main_slavery_taurus.tex");  save_csv(t_slv_tau$table,  "table04_main_slavery_taurus.csv")

#####################################################################################################################################
######################################## TEMPLATE-STYLE TABLES WITH ALL CONTROLS ####################################################
#####################################################################################################################################

# Function to get actual variable names from a model
get_model_terms <- function(model) {
  if (is.null(model)) return(character(0))
  tt <- broom::tidy(model)
  return(tt$term)
}

# Pretty labels for variables - comprehensive mapping
pretty_var <- function(term) {
  # Create a comprehensive mapping of potential variable names
  recode_map <- c(
    # Instrument variables
    "ss_imports" = "Shift-share (Imports)",
    "ss_taurus" = "Shift-share (Taurus)",
    
    # Control variables - various possible formats
    "log(mean_harvested_area_interp+1)" = "Log harvested area",
    "log(mean_harvested_area_interp + 1)" = "Log harvested area",
    "mean_harvested_area_interp" = "Log harvested area",
    
    "log(gdp_per_capita+1)" = "Log GDP per capita",
    "log(gdp_per_capita + 1)" = "Log GDP per capita", 
    "gdp_per_capita" = "Log GDP per capita",
    
    "log(mean_number_cattle_interp+1)" = "Log cattle herd",
    "log(mean_number_cattle_interp + 1)" = "Log cattle herd",
    "mean_number_cattle_interp" = "Log cattle herd",
    
    "log(mean_revenue_interp+1)" = "Log municipal revenue",
    "log(mean_revenue_interp + 1)" = "Log municipal revenue",
    "mean_revenue_interp" = "Log municipal revenue"
  )
  
  # Return mapped name if exists, otherwise return original
  ifelse(term %in% names(recode_map), recode_map[term], term)
}

# Enhanced function to extract coefficient and standard error for a term
term_cells <- function(model, term, digits = 3) {
  if (is.null(model)) return(c("", ""))
  
  tt <- broom::tidy(model)
  
  # Try exact match first
  r <- tt[tt$term == term, , drop = FALSE]
  
  # If no exact match, try variations
  if (nrow(r) == 0) {
    # Try with different spacing around +
    alt_terms <- c(
      gsub("\\+", " + ", term),  # Add spaces around +
      gsub(" \\+ ", "+", term),  # Remove spaces around +
      gsub("log\\((.*)\\+1\\)", "\\1", term),  # Remove log wrapper
      gsub("log\\((.*)\\)", "\\1", term)  # Remove log wrapper (alternative)
    )
    
    for (alt_term in alt_terms) {
      r <- tt[tt$term == alt_term, , drop = FALSE]
      if (nrow(r) > 0) break
    }
  }
  
  if (nrow(r) == 0) return(c("", ""))  # Return empty if still no match
  
  est <- sprintf(paste0("%.", digits, "f"), r$estimate[1])
  se  <- sprintf(paste0("(%.", digits, "f)"), r$std.error[1])
  c(paste0(est, stars(r$p.value[1])), se)
}

# Write template-style LaTeX table
write_template_table <- function(models, dep_label, caption, label, instrument_term, control_terms, file_name, data_for_counts) {
  K <- length(models)
  
  # Debug: Print actual terms from first model with controls
  cat("\n=== DEBUG: Actual terms in model 2 ===\n")
  if (length(models) >= 2 && !is.null(models[[2]])) {
    actual_terms <- get_model_terms(models[[2]])
    cat("All terms:", paste(actual_terms, collapse = ", "), "\n")
    
    # Filter out fixed effects and intercept - be more specific
    control_candidates <- actual_terms[!grepl("\\(Intercept\\)", actual_terms)]
    control_candidates <- control_candidates[control_candidates != instrument_term]
    cat("Control candidates after filtering:", paste(control_candidates, collapse = ", "), "\n")
    
    # Keep the original control_terms but show what we found
    cat("Original control_terms:", paste(control_terms, collapse = ", "), "\n")
    
    # Try to match each original control term with actual terms
    matched_controls <- character()
    for (orig_term in control_terms) {
      # Try exact match first
      if (orig_term %in% actual_terms) {
        matched_controls <- c(matched_controls, orig_term)
        cat("Exact match for:", orig_term, "\n")
      } else {
        # Try variations
        variations <- c(
          gsub("\\+", " + ", orig_term),  # Add spaces around +
          gsub(" \\+ ", "+", orig_term),  # Remove spaces around +
          gsub("log\\((.*)\\+1\\)", "\\1", orig_term)  # Remove log wrapper
        )
        
        found_match <- FALSE
        for (var_term in variations) {
          if (var_term %in% actual_terms) {
            matched_controls <- c(matched_controls, var_term)
            cat("Variation match for:", orig_term, "->", var_term, "\n")
            found_match <- TRUE
            break
          }
        }
        
        if (!found_match) {
          cat("No match found for:", orig_term, "\n")
        }
      }
    }
    
    # Use matched controls if we found any, otherwise use all candidates
    if (length(matched_controls) > 0) {
      control_terms <- matched_controls
      cat("Using matched controls:", paste(control_terms, collapse = ", "), "\n")
    } else {
      control_terms <- control_candidates
      cat("Using all candidates:", paste(control_terms, collapse = ", "), "\n")
    }
  }
  
  headers <- paste0(" & ", paste0("(", seq_len(K), ")", collapse = " & "), " \\\\")
  top <- c("\\begin{table}[ht]",
           "\\begin{center}",
           paste0("\\caption{", caption, "}"),
           "\\resizebox{12cm}{!} {",
           paste0("\\begin{tabular}{l", paste(rep("c", K), collapse = ""), "}"),
           "\\toprule",
           headers,
           "\\midrule",
           paste0("\\multicolumn{", K+1, "}{l}{\\textbf{\\emph{Dependent Variable: ", dep_label, "}}} \\\\"),
           "\\midrule")
  
  vars <- c(instrument_term, control_terms)
  
  body <- character()
  for (term in vars) {
    cat("Processing term:", term, "\n")
    row_coef <- c(pretty_var(term))
    row_se   <- c("")
    for (i in seq_along(models)) {
      m <- models[[i]]
      cells <- term_cells(m, term, digits = 3)
      cat("  Model", i, "- cells:", cells[1], cells[2], "\n")
      row_coef <- c(row_coef, cells[1])
      row_se   <- c(row_se,   cells[2])
    }
    body <- c(body,
              paste(row_coef, collapse = " & "), " \\\\",
              paste(row_se,   collapse = " & "), " \\\\",
              "")
  }
  
  r2s   <- vapply(models, safe_r2_within, numeric(1))
  r2s   <- ifelse(is.na(r2s), "NA", sprintf("%.3f", r2s))
  Ns    <- vapply(models, function(m) m$nobs, numeric(1))
  Clust <- rep(dplyr::n_distinct(data_for_counts$ibge_cod), K)
  
  summ <- c("\\hline",
            paste("R-squared", paste(c("", r2s), collapse = " & "), "\\\\"),
            paste("Municipality clusters", paste(c("", Clust), collapse = " & "), "\\\\"),
            paste("N", paste(c("", Ns), collapse = " & "), "\\\\"),
            "\\hline",
            paste("Year fixed effects", paste(c("", rep("YES", K)), collapse = " & "), "\\\\"),
            paste("Municipality fixed effects", paste(c("", rep("YES", K)), collapse = " & "), "\\\\"),
            paste("State linear trend", paste(c("", rep("", K-1), "YES"), collapse = " & "), "\\\\"),
            "\\bottomrule")
  
  notes <- c(paste0("\\multicolumn{", K+1, "}{p{12cm}}{\\footnotesize \\emph{Notes}: ",
                    "All specifications include municipality and year fixed effects. ",
                    "Robust standard errors clustered at the municipality level in parentheses. ",
                    "$ + p < .1, * p < .05, ** p < .01, *** p < .001. $}"),
             "\\end{tabular}%",
             "}",
             paste0("\\label{", label, "}"),
             "\\end{center}",
             "\\end{table}%")
  
  lines <- c(top, body, summ, notes)
  out_path <- file.path(out_dir, file_name)
  writeLines(lines, con = out_path)
  created_files <<- c(created_files, out_path)
}

# Initial control terms (will be updated automatically in write_template_table)
controls_initial <- c(
  "log(mean_harvested_area_interp+1)",
  "log(gdp_per_capita+1)",
  "log(mean_number_cattle_interp+1)",
  "log(mean_revenue_interp+1)"
)

# Four template-style tables with all controls visible
write_template_table(t_ops_imp$models, "Labor Audits",
                     "Effect of Gun Imports Shift-Share on Labor Audits",
                     "tab:template_ops_imports", "ss_imports", controls_initial,
                     "table01_main_audits_imports_template.tex", df_clean)

write_template_table(t_ops_tau$models, "Labor Audits",
                     "Effect of Taurus Revenue Shift-Share on Labor Audits",
                     "tab:template_ops_taurus", "ss_taurus", controls_initial,
                     "table02_main_audits_taurus_template.tex", df_clean)

write_template_table(t_slv_imp$models, "Slavery Cases",
                     "Effect of Gun Imports Shift-Share on Slavery Cases",
                     "tab:template_slavery_imports", "ss_imports", controls_initial,
                     "table03_main_slavery_imports_template.tex", df_clean)

write_template_table(t_slv_tau$models, "Slavery Cases",
                     "Effect of Taurus Revenue Shift-Share on Slavery Cases",
                     "tab:template_slavery_taurus", "ss_taurus", controls_initial,
                     "table04_main_slavery_taurus_template.tex", df_clean)

#####################################################################################################################################
######################################## STATE-BY-STATE LEAVE-ONE-OUT (UF) ##########################################################
#####################################################################################################################################

leave_one_out <- function(outcome_var, instrument_var, data) {
  states <- sort(unique(data$state_code))
  results <- vector("list", length(states))
  names(results) <- states
  
  for (i in seq_along(states)) {
    uf <- states[i]
    cat("LOO:", outcome_var, "~", instrument_var, "— excluding state", uf, "\n")
    dt  <- data[data$state_code != uf, ]
    fml <- as.formula(paste0(
      outcome_var, " ~ ", instrument_var, " + ",
      "log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + ",
      "log(mean_number_cattle_interp+1) + log(mean_revenue_interp+1) | ibge_cod + year"
    ))
    m <- tryCatch(feols(fml, data = dt, panel.id = ~ ibge_cod + year, cluster = ~ ibge_cod),
                  error = function(e) NULL)
    if (is.null(m)) next
    cs <- broom::tidy(m)
    mc <- cs[cs$term == instrument_var, , drop = FALSE]
    if (!nrow(mc)) next
    results[[i]] <- data.frame(
      excluded_state = uf,
      coefficient = round(mc$estimate[1], 4),
      std_error  = round(mc$std.error[1], 4),
      p_value    = round(mc$p.value[1], 4),
      n_obs      = m$nobs,
      r_squared  = round(safe_r2_within(m), 4)
    )
  }
  
  out <- dplyr::bind_rows(results)
  if (nrow(out)) {
    out$ci_lower <- round(out$coefficient - 1.96 * out$std_error, 4)
    out$ci_upper <- round(out$coefficient + 1.96 * out$std_error, 4)
  }
  out
}

cat("\n================================================================================\n")
cat("LEAVE-ONE-OUT ROBUSTNESS (STATE-BY-STATE)\n")
cat("================================================================================\n\n")

loo_ops_imp  <- leave_one_out("operations", "ss_imports", df_clean)
loo_ops_tau  <- leave_one_out("operations", "ss_taurus",  df_clean)
loo_slv_imp  <- leave_one_out("slave",      "ss_imports", df_clean)
loo_slv_tau  <- leave_one_out("slave",      "ss_taurus",  df_clean)

if (nrow(loo_ops_imp)) {
  xt <- xtable(loo_ops_imp,
               caption = "Leave-One-State-Out: Labor Audits Response to Gun Imports",
               label = "tab:loo_audits_imports",
               digits = c(0,0,4,4,4,0,4,4,4))
  save_tex(xt, "table05_loo_audits_imports.tex")
  save_csv(loo_ops_imp, "table05_loo_audits_imports.csv")
}
if (nrow(loo_ops_tau)) {
  xt <- xtable(loo_ops_tau,
               caption = "Leave-One-State-Out: Labor Audits Response to Taurus Revenue",
               label = "tab:loo_audits_taurus",
               digits = c(0,0,4,4,4,0,4,4,4))
  save_tex(xt, "table06_loo_audits_taurus.tex")
  save_csv(loo_ops_tau, "table06_loo_audits_taurus.csv")
}
if (nrow(loo_slv_imp)) {
  xt <- xtable(loo_slv_imp,
               caption = "Leave-One-State-Out: Slavery Cases Response to Gun Imports",
               label = "tab:loo_slavery_imports",
               digits = c(0,0,4,4,4,0,4,4,4))
  save_tex(xt, "table07_loo_slavery_imports.tex")
  save_csv(loo_slv_imp, "table07_loo_slavery_imports.csv")
}
if (nrow(loo_slv_tau)) {
  xt <- xtable(loo_slv_tau,
               caption = "Leave-One-State-Out: Slavery Cases Response to Taurus Revenue",
               label = "tab:loo_slavery_taurus",
               digits = c(0,0,4,4,4,0,4,4,4))
  save_tex(xt, "table08_loo_slavery_taurus.tex")
  save_csv(loo_slv_tau, "table08_loo_slavery_taurus.csv")
}

#####################################################################################################################################
######################################## COMPREHENSIVE ROBUSTNESS CHECKS ############################################################
#####################################################################################################################################

create_robustness_table <- function(outcome_var, instrument_var, data) {
  
  base_fml <- as.formula(paste0(
    outcome_var, " ~ ", instrument_var, " + ",
    "log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + ",
    "log(mean_number_cattle_interp+1) + log(mean_revenue_interp+1) | ibge_cod + year"
  ))
  
  # (1) Base + Controls
  m_base <- feols(base_fml, panel.id = ~ ibge_cod + year, cluster = ~ ibge_cod, data = data)
  
  # (2) + State linear trend
  trend_fml <- as.formula(paste0(
    outcome_var, " ~ ", instrument_var, " + ",
    "log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + ",
    "log(mean_number_cattle_interp+1) + log(mean_revenue_interp+1) | ibge_cod + year + state_code[year]"
  ))
  m_trend <- feols(trend_fml, panel.id = ~ ibge_cod + year, cluster = ~ ibge_cod, data = data)
  
  # (3) Conley SE (100 km) via fixest::vcov_conley
  m_conley <- tryCatch({
    feols(base_fml,
          panel.id = ~ ibge_cod + year,
          data = data,
          vcov = fixest::vcov_conley(spatial = ~ longitude + latitude, time = ~ year, dist_cutoff = 100))
  }, error = function(e) NULL)
  
  # (4) Drop Rural Top 10%
  thr90 <- quantile(data$pop_rural_2000, 0.9, na.rm = TRUE)
  m_r90 <- feols(base_fml, panel.id = ~ ibge_cod + year, cluster = ~ ibge_cod,
                 data = dplyr::filter(data, pop_rural_2000 <= thr90))
  
  # (5) Drop Rural Top 20%
  thr80 <- quantile(data$pop_rural_2000, 0.8, na.rm = TRUE)
  m_r80 <- feols(base_fml, panel.id = ~ ibge_cod + year, cluster = ~ ibge_cod,
                 data = dplyr::filter(data, pop_rural_2000 <= thr80))
  
  # (6) Poisson (PPML-FE)
  m_ppml <- tryCatch({
    if (outcome_var %in% c("operations", "slave")) {
      fepois(base_fml, panel.id = ~ ibge_cod + year, cluster = ~ ibge_cod, data = data)
    } else NULL
  }, error = function(e) NULL)
  
  grab <- function(m) {
    if (is.null(m)) return(c(NA,NA,NA,NA,NA))
    cs <- broom::tidy(m)
    mc <- cs[cs$term == instrument_var, , drop = FALSE]
    if (!nrow(mc)) return(c(NA,NA,NA,NA,NA))
    c(round(mc$estimate[1], 4),
      round(mc$std.error[1], 4),
      round(mc$p.value[1],  4),
      m$nobs,
      round(safe_r2_within(m), 4))
  }
  
  res <- rbind(
    grab(m_base),
    grab(m_trend),
    grab(m_conley),
    grab(m_r90),
    grab(m_r80),
    grab(m_ppml)
  )
  
  data.frame(
    Specification = c("Base + Controls",
                      "Base + State linear trend",
                      "Base + Conley SE (100km)",
                      "Drop Rural Top 10%",
                      "Drop Rural Top 20%",
                      "Poisson"),
    Coefficient = res[,1],
    Std_Error  = res[,2],
    P_Value    = res[,3],
    N_Obs      = res[,4],
    R_Squared  = res[,5]
  )
}

cat("\n================================================================================\n")
cat("COMPREHENSIVE ROBUSTNESS CHECKS\n")
cat("================================================================================\n\n")

rob_ops_imp  <- create_robustness_table("operations", "ss_imports", df_clean)
rob_ops_tau  <- create_robustness_table("operations", "ss_taurus",  df_clean)
rob_slv_imp  <- create_robustness_table("slave",      "ss_imports", df_clean)
rob_slv_tau  <- create_robustness_table("slave",      "ss_taurus",  df_clean)

# Save Tables 9–12
save_tex(xtable(rob_ops_imp,
                caption="Robustness: Labor Audits Response to Gun Imports Shift-Share",
                label="tab:rob_audits_imports",
                digits=c(0,0,4,4,4,0,4)),
         "table09_robustness_audits_imports.tex")
save_csv(rob_ops_imp, "table09_robustness_audits_imports.csv")

save_tex(xtable(rob_ops_tau,
                caption="Robustness: Labor Audits Response to Taurus Revenue Shift-Share",
                label="tab:rob_audits_taurus",
                digits=c(0,0,4,4,4,0,4)),
         "table10_robustness_audits_taurus.tex")
save_csv(rob_ops_tau, "table10_robustness_audits_taurus.csv")

save_tex(xtable(rob_slv_imp,
                caption="Robustness: Slavery Cases Response to Gun Imports Shift-Share",
                label="tab:rob_slavery_imports",
                digits=c(0,0,4,4,4,0,4)),
         "table11_robustness_slavery_imports.tex")
save_csv(rob_slv_imp, "table11_robustness_slavery_imports.csv")

save_tex(xtable(rob_slv_tau,
                caption="Robustness: Slavery Cases Response to Taurus Revenue Shift-Share",
                label="tab:rob_slavery_taurus",
                digits=c(0,0,4,4,4,0,4)),
         "table12_robustness_slavery_taurus.tex")
save_csv(rob_slv_tau, "table12_robustness_slavery_taurus.csv")

#####################################################################################################################################
######################################## FINAL SUMMARY ##############################################################################
#####################################################################################################################################

cat("\n================================================================================\n")
cat("FILES WRITTEN TO:\n"); cat(out_dir, "\n")
cat("================================================================================\n")
# Show only files produced in this run (that still exist)
created_files <- unique(created_files[file.exists(created_files)])
print(created_files)
