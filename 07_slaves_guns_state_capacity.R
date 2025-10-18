# =========================================================
# Streamlined Replication in R — Reviewer Response
# Focus: Municipal Capacity & "False Leads" Tests
# Data:   C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Data/slavery_guns_census_panel.rds
# Output: C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Results
# Models: Linear only (feols: LPM)
# Table:  One big LaTeX + HTML
# =========================================================

## ---- 0) Packages ----
req <- c("data.table","fixest","modelsummary","dplyr","stringr")
new <- req[!req %in% installed.packages()[,"Package"]]
if (length(new)) install.packages(new, dependencies = TRUE, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))
fixest::setFixest_nthreads(1)

## ---- 1) Paths ----
path_data <- "C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Data/slavery_guns_census_panel.rds"
path_out  <- "C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Results"
if (!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)

tex_file  <- file.path(path_out, "main_results_summary.tex")
html_file <- file.path(path_out, "main_results_summary.html")

## ---- 2) Load & prepare data ----
dt <- readRDS(path_data)
data.table::setDT(dt)

# Ensure required columns
need <- c("ibge_cod","year","mean_uf_cod",
          "operations","slave","ss_taurus",
          "mean_harvested_area_interp","gdp_per_capita",
          "mean_number_cattle_interp","mean_revenue_interp",
          "p_desp_funcao_seguranca_mun","p_idhm")
miss <- setdiff(need, names(dt))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

# Analysis variables
dt <- dt |>
  mutate(
    l_harv    = log(mean_harvested_area_interp + 1),
    l_gdppc   = log(gdp_per_capita + 1),
    l_cattle  = log(mean_number_cattle_interp + 1),
    l_revenue = log(mean_revenue_interp + 1),
    audited   = as.integer(operations > 0),
    hit_any   = ifelse(operations > 0, as.integer(slave > 0), NA_integer_)
  )

# Capacity proxies (standardized)
cap_vars <- c("p_desp_funcao_seguranca_mun","p_idhm","mean_revenue_interp")
cap_vars <- intersect(cap_vars, names(dt))
for (v in cap_vars) dt[[paste0(v,"_z")]] <- as.numeric(scale(dt[[v]]))

# Proxy choice
cap_choice <- dplyr::case_when(
  "p_desp_funcao_seguranca_mun_z" %in% names(dt) ~ "p_desp_funcao_seguranca_mun_z",
  "p_idhm_z"                        %in% names(dt) ~ "p_idhm_z",
  "mean_revenue_interp_z"           %in% names(dt) ~ "mean_revenue_interp_z",
  TRUE ~ NA_character_
)

## ---- 3) Models ----
# Baseline
reg8 <- feols(
  slave ~ ss_taurus + l_harv + l_gdppc + l_cattle + l_revenue | ibge_cod + year,
  cluster = ~ ibge_cod, data = dt
)

# False-leads (audited only)
dt_ops <- dt[operations > 0]

m_hit_lpm <- feols(
  hit_any ~ ss_taurus + l_harv + l_gdppc + l_cattle + l_revenue | ibge_cod + year,
  cluster = ~ ibge_cod, data = dt_ops
)

# UF × year FE
m_hit_lpm_ufy <- feols(
  hit_any ~ ss_taurus + l_harv + l_gdppc + l_cattle + l_revenue |
    ibge_cod + year + mean_uf_cod^year,
  cluster = ~ ibge_cod, data = dt_ops
)

# Capacity interactions
fit_hit_int <- NULL
if (!is.na(cap_choice)) {
  f_hit  <- as.formula(paste0(
    "hit_any ~ ss_taurus * ", cap_choice,
    " + l_harv + l_gdppc + l_cattle + l_revenue | ibge_cod + year"
  ))
  fit_hit_int <- feols(f_hit, cluster = ~ ibge_cod, data = dt_ops)
}

# Audits outcome
m_ops <- feols(
  operations ~ ss_taurus + l_harv + l_gdppc + l_cattle + l_revenue | ibge_cod + year,
  cluster = ~ ibge_cod, data = dt
)
m_ops_ufy <- feols(
  operations ~ ss_taurus + l_harv + l_gdppc + l_cattle + l_revenue |
    ibge_cod + year + mean_uf_cod^year,
  cluster = ~ ibge_cod, data = dt
)

## ---- 4) One big table (LaTeX + HTML) ----
options(modelsummary_format_numeric_latex = "plain")

mods <- list(
  "Slavery (baseline FE)"           = reg8,
  "Hit (LPM | audited)"             = m_hit_lpm,
  "Hit (LPM | audited, UF×year FE)" = m_hit_lpm_ufy,
  "Audits (count, FE)"              = m_ops,
  "Audits (count, UF×year FE)"      = m_ops_ufy
)
if (!is.null(fit_hit_int)) mods[["Hit (LPM) × Capacity"]] <- fit_hit_int

# Flags for FE
ufy_names <- c("Hit (LPM | audited, UF×year FE)","Audits (count, UF×year FE)")
fe_flags <- data.frame(rowname = c("Municipality FE","Year FE","UF×Year FE","Cluster SE"),
                       check.names = FALSE)
for (nm in names(mods)) {
  fe_flags[[nm]] <- c("Yes","Yes", ifelse(nm %in% ufy_names, "Yes", "No"), "Municipality")
}

# Coefficient labels
coef_map <- c(
  "ss_taurus" = "Shift–Share (Taurus)",
  "l_harv"    = "log(Harvested area + 1)",
  "l_gdppc"   = "log(GDP per capita + 1)",
  "l_cattle"  = "log(Cattle + 1)",
  "l_revenue" = "log(Municipal revenue + 1)"
)

# LaTeX
modelsummary(
  mods, coef_map = coef_map,
  statistic = "({std.error})", stars = TRUE,
  gof_omit  = "IC|Log|Within|Pseudo|Std|F|Theta",
  add_rows  = fe_flags, output = tex_file
)

# HTML
modelsummary(
  mods, coef_map = coef_map,
  statistic = "({std.error})", stars = TRUE,
  gof_omit  = "IC|Log|Within|Pseudo|Std|F|Theta",
  add_rows  = fe_flags, output = html_file
)

cat("Tables saved to:\n", tex_file, "\n", html_file, "\n")
