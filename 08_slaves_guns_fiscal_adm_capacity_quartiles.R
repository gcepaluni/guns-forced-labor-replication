##########################################################################################################
######################## Revenue/GDP Quartiles: Shift-share Imports ######################################
##########################################################################################################

library(dplyr)
library(fixest)
library(broom)
library(writexl)
library(tidylog)

# --------------------------------------------------------------------------------
# 0. Load panel data
# --------------------------------------------------------------------------------
panel <- readRDS("C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Data/slavery_guns_census_panel_with_statecapacity_cleaned.rds")

# ----------------------------------------------------------------------------------------
# 1. Create quartiles for Revenue/GDP ratio (baseline 2000)
# ----------------------------------------------------------------------------------------
baseline <- panel %>% 
  filter(year == 2000 & !is.na(mean_revenue_interp) & !is.na(gdp_per_capita)) %>%
  mutate(
    total_gdp = gdp_per_capita * 1000,  # Adjust multiplier if needed
    revenue_gdp_ratio = mean_revenue_interp / total_gdp
  ) %>%
  select(ibge_cod, revenue_gdp_ratio) %>%
  rename(rev_gdp_2000 = revenue_gdp_ratio)

quartile_cutoffs <- quantile(baseline$rev_gdp_2000,
                             probs = c(0, 0.25, 0.5, 0.75, 1),
                             na.rm = TRUE)

panel <- panel %>%
  left_join(baseline, by = 'ibge_cod') %>% 
  mutate(
    revgdp_quartile = cut(
      rev_gdp_2000,
      breaks = quartile_cutoffs,
      include.lowest = TRUE,
      labels = c("Low", "Medium-Low", "Medium-High", "High")
    )
  )

cat("\nQuartile distribution (municipalities by baseline Revenue/GDP ratio):\n")
print(table(panel$revgdp_quartile, useNA = "ifany"))

# --------------------------------------------------------------------------------
# 2. Run regressions by quartile
# --------------------------------------------------------------------------------
results_list <- list()

for (q in levels(panel$revgdp_quartile)) {
  subdata <- panel %>% filter(revgdp_quartile == q)
  
  if (nrow(subdata) < 200) {
    warning(paste("Skipping", q, "- too few observations"))
    next
  }
  
  mod_audits <- feols(
    operations ~ ss_imports +
      log(mean_harvested_area_interp + 1) +
      log(gdp_per_capita + 1) +
      log(mean_number_cattle_interp + 1) +
      log(mean_revenue_interp + 1) | ibge_cod + year,
    data = subdata, panel.id = ~ibge_cod + year, cluster = ~ibge_cod
  )
  
  mod_slaves <- feols(
    slave ~ ss_imports +
      log(mean_harvested_area_interp + 1) +
      log(gdp_per_capita + 1) +
      log(mean_number_cattle_interp + 1) +
      log(mean_revenue_interp + 1) | ibge_cod + year,
    data = subdata, panel.id = ~ibge_cod + year, cluster = ~ibge_cod
  )
  
  results_list[[q]] <- list(audits = mod_audits, slaves = mod_slaves)
}

# --------------------------------------------------------------------------------
# 3. Extract tidy results
# --------------------------------------------------------------------------------
results_table <- bind_rows(
  lapply(names(results_list), function(q) {
    bind_rows(
      tidy(results_list[[q]]$audits, conf.int = TRUE) %>%
        filter(term == "ss_imports") %>%
        mutate(outcome = "Audits", quartile = q,
               n_obs = nobs(results_list[[q]]$audits),
               r2 = as.numeric(fitstat(results_list[[q]]$audits, "wr2"))),
      tidy(results_list[[q]]$slaves, conf.int = TRUE) %>%
        filter(term == "ss_imports") %>%
        mutate(outcome = "Slaves", quartile = q,
               n_obs = nobs(results_list[[q]]$slaves),
               r2 = as.numeric(fitstat(results_list[[q]]$slaves, "wr2")))
    )
  })
)

# --------------------------------------------------------------------------------
# 4. Export results
# --------------------------------------------------------------------------------
out_dir <- "C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Results"
csv_file <- file.path(out_dir, "revgdp_quartile_imports_results.csv")
xlsx_file <- file.path(out_dir, "revgdp_quartile_imports_results.xlsx")

write.csv(results_table, csv_file, row.names = FALSE)
write_xlsx(list("Coefficient Estimates" = results_table), xlsx_file)

cat(sprintf("\n✓ Results saved: %s\n", csv_file))
cat(sprintf("✓ Results saved: %s\n", xlsx_file))

# --------------------------------------------------------------------------------
# 5. LaTeX table (exact template format)
# --------------------------------------------------------------------------------
tex_file <- file.path(out_dir, "revgdp_quartile_imports_results.tex")

# Helper function for significance stars
add_stars <- function(pval) {
  ifelse(pval < 0.001, "$^{***}$",
         ifelse(pval < 0.01, "$^{**}$",
                ifelse(pval < 0.05, "$^{*}$",
                       ifelse(pval < 0.1, "$^{+}$", ""))))
}

# Build rows for each outcome
build_table_rows <- function(outcome_label, results_df) {
  subset <- results_df %>% 
    filter(outcome == outcome_label) %>% 
    arrange(factor(quartile, levels = c("Low", "Medium-Low", "Medium-High", "High")))
  
  coefs <- sprintf("%.3f", subset$estimate)
  ses <- sprintf("(%.3f)", subset$std.error)
  stars <- add_stars(subset$p.value)
  n_obs <- subset$n_obs
  r2_vals <- sprintf("%.3f", subset$r2)
  
  coef_row <- paste0("Shift--Share (Imports) & ", 
                     paste0(coefs, stars, collapse = " & "), " \\\\")
  se_row <- paste0(" & ", paste(ses, collapse = " & "), " \\\\")
  n_row <- paste0("Observations & ", paste(n_obs, collapse = " & "), " \\\\")
  r2_row <- paste0("$R^2$ & ", paste(r2_vals, collapse = " & "), " \\\\")
  ctrl_row <- paste0("Controls & ", paste(rep("Yes", 4), collapse = " & "), " \\\\")
  
  return(list(
    header = paste0("\\multicolumn{5}{l}{\\textbf{", outcome_label, "}} \\\\"),
    coef = coef_row,
    se = se_row,
    n = n_row,
    r2 = r2_row,
    ctrl = ctrl_row
  ))
}

audits <- build_table_rows("Audits", results_table)
slaves <- build_table_rows("Slaves", results_table)

tex_content <- sprintf("\\begin{table}[!htbp]
\\begin{center}
\\small
\\caption{Heterogeneous Effects of Gun Imports by Municipal Revenue/GDP Ratio Quartiles}
\\label{tab:revgdp_quartiles}
\\begin{tabular}{lcccc}
\\toprule
 & Low & Medium-Low & Medium-High & High \\\\
\\midrule
%s
%s
%s
%s
%s
%s
\\midrule
%s
%s
%s
%s
%s
%s
\\bottomrule
\\end{tabular}
\\end{center}
\\footnotesize \\emph{Notes}: Municipalities are grouped into quartiles of baseline (year 2000) ratio of municipal revenue to GDP. All specifications include municipality and year fixed effects, with controls for log harvested area, log GDP per capita, log cattle herd, and log municipal revenue. Standard errors clustered at the municipal level are reported in parentheses. $^{+}p<0.10$, $^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.
\\end{table}
",
                       audits$header, audits$coef, audits$se, audits$n, audits$r2, audits$ctrl,
                       slaves$header, slaves$coef, slaves$se, slaves$n, slaves$r2, slaves$ctrl
)

writeLines(tex_content, tex_file)
cat(sprintf("✓ LaTeX table saved: %s\n", tex_file))


##########################################################################################################
######################## HDI Quartiles: Shift-share Imports ##############################################
##########################################################################################################

library(dplyr)
library(fixest)
library(broom)
library(writexl)
library(tidylog)

# --------------------------------------------------------------------------------
# 0. Load panel data
# --------------------------------------------------------------------------------
panel <- readRDS("C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Data/slavery_guns_census_panel_with_statecapacity_cleaned.rds")

# ----------------------------------------------------------------------------------------
# 1. Create quartiles for HDI (baseline 2000)
# ----------------------------------------------------------------------------------------
baseline <- panel %>% 
  filter(year == 2000 & !is.na(p_idhm_interp)) %>%
  select(ibge_cod, p_idhm_interp) %>%
  rename(hdi_2000 = p_idhm_interp)

quartile_cutoffs <- quantile(baseline$hdi_2000,
                             probs = c(0, 0.25, 0.5, 0.75, 1),
                             na.rm = TRUE)

panel <- panel %>%
  left_join(baseline, by = 'ibge_cod') %>% 
  mutate(
    hdi_quartile = cut(
      hdi_2000,
      breaks = quartile_cutoffs,
      include.lowest = TRUE,
      labels = c("Low", "Medium-Low", "Medium-High", "High")
    )
  )

cat("\nQuartile distribution (municipalities by baseline HDI):\n")
print(table(panel$hdi_quartile, useNA = "ifany"))

# --------------------------------------------------------------------------------
# 2. Run regressions by quartile
# --------------------------------------------------------------------------------
results_list <- list()

for (q in levels(panel$hdi_quartile)) {
  subdata <- panel %>% filter(hdi_quartile == q)
  
  if (nrow(subdata) < 200) {
    warning(paste("Skipping", q, "- too few observations"))
    next
  }
  
  mod_audits <- feols(
    operations ~ ss_imports +
      log(mean_harvested_area_interp + 1) +
      log(gdp_per_capita + 1) +
      log(mean_number_cattle_interp + 1) +
      log(mean_revenue_interp + 1) | ibge_cod + year,
    data = subdata, panel.id = ~ibge_cod + year, cluster = ~ibge_cod
  )
  
  mod_slaves <- feols(
    slave ~ ss_imports +
      log(mean_harvested_area_interp + 1) +
      log(gdp_per_capita + 1) +
      log(mean_number_cattle_interp + 1) +
      log(mean_revenue_interp + 1) | ibge_cod + year,
    data = subdata, panel.id = ~ibge_cod + year, cluster = ~ibge_cod
  )
  
  results_list[[q]] <- list(audits = mod_audits, slaves = mod_slaves)
}

# --------------------------------------------------------------------------------
# 3. Extract tidy results
# --------------------------------------------------------------------------------
results_table <- bind_rows(
  lapply(names(results_list), function(q) {
    bind_rows(
      tidy(results_list[[q]]$audits, conf.int = TRUE) %>%
        filter(term == "ss_imports") %>%
        mutate(outcome = "Audits", quartile = q,
               n_obs = nobs(results_list[[q]]$audits),
               r2 = as.numeric(fitstat(results_list[[q]]$audits, "wr2"))),
      tidy(results_list[[q]]$slaves, conf.int = TRUE) %>%
        filter(term == "ss_imports") %>%
        mutate(outcome = "Slaves", quartile = q,
               n_obs = nobs(results_list[[q]]$slaves),
               r2 = as.numeric(fitstat(results_list[[q]]$slaves, "wr2")))
    )
  })
)

# --------------------------------------------------------------------------------
# 4. Export results
# --------------------------------------------------------------------------------
out_dir <- "C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Results"
csv_file <- file.path(out_dir, "hdi_quartile_imports_results.csv")
xlsx_file <- file.path(out_dir, "hdi_quartile_imports_results.xlsx")

write.csv(results_table, csv_file, row.names = FALSE)
write_xlsx(list("Coefficient Estimates" = results_table), xlsx_file)

cat(sprintf("\n✓ Results saved: %s\n", csv_file))
cat(sprintf("✓ Results saved: %s\n", xlsx_file))

# --------------------------------------------------------------------------------
# 5. LaTeX table (exact template format)
# --------------------------------------------------------------------------------
tex_file <- file.path(out_dir, "hdi_quartile_imports_results.tex")

# Helper function for significance stars
add_stars <- function(pval) {
  ifelse(pval < 0.001, "$^{***}$",
         ifelse(pval < 0.01, "$^{**}$",
                ifelse(pval < 0.05, "$^{*}$",
                       ifelse(pval < 0.1, "$^{+}$", ""))))
}

# Build rows for each outcome
build_table_rows <- function(outcome_label, results_df) {
  subset <- results_df %>% 
    filter(outcome == outcome_label) %>% 
    arrange(factor(quartile, levels = c("Low", "Medium-Low", "Medium-High", "High")))
  
  coefs <- sprintf("%.3f", subset$estimate)
  ses <- sprintf("(%.3f)", subset$std.error)
  stars <- add_stars(subset$p.value)
  n_obs <- subset$n_obs
  r2_vals <- sprintf("%.3f", subset$r2)
  
  coef_row <- paste0("Shift--Share (Imports) & ", 
                     paste0(coefs, stars, collapse = " & "), " \\\\")
  se_row <- paste0(" & ", paste(ses, collapse = " & "), " \\\\")
  n_row <- paste0("Observations & ", paste(n_obs, collapse = " & "), " \\\\")
  r2_row <- paste0("$R^2$ & ", paste(r2_vals, collapse = " & "), " \\\\")
  ctrl_row <- paste0("Controls & ", paste(rep("Yes", 4), collapse = " & "), " \\\\")
  
  return(list(
    header = paste0("\\multicolumn{5}{l}{\\textbf{", outcome_label, "}} \\\\"),
    coef = coef_row,
    se = se_row,
    n = n_row,
    r2 = r2_row,
    ctrl = ctrl_row
  ))
}

audits <- build_table_rows("Audits", results_table)
slaves <- build_table_rows("Slaves", results_table)

tex_content <- sprintf("\\begin{table}[!htbp]
\\begin{center}
\\small
\\caption{Heterogeneous Effects of Gun Imports by Municipal Human Development Index Quartiles}
\\label{tab:hdi_quartiles}
\\begin{tabular}{lcccc}
\\toprule
 & Low & Medium-Low & Medium-High & High \\\\
\\midrule
%s
%s
%s
%s
%s
%s
\\midrule
%s
%s
%s
%s
%s
%s
\\bottomrule
\\end{tabular}
\\end{center}
\\footnotesize \\emph{Notes}: Municipalities are grouped into quartiles of baseline (year 2000) Municipal Human Development Index (IDHM). All specifications include municipality and year fixed effects, with controls for log harvested area, log GDP per capita, log cattle herd, and log municipal revenue. Standard errors clustered at the municipal level are reported in parentheses. $^{+}p<0.10$, $^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.
\\end{table}
",
                       audits$header, audits$coef, audits$se, audits$n, audits$r2, audits$ctrl,
                       slaves$header, slaves$coef, slaves$se, slaves$n, slaves$r2, slaves$ctrl
)

writeLines(tex_content, tex_file)
cat(sprintf("✓ LaTeX table saved: %s\n", tex_file))