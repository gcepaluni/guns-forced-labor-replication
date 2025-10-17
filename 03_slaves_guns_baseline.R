########################################################################################################
############ Reviewer Request: Baseline(2000) trends & optional UF×Year fixed effects (compact) ########
########################################################################################################
# This script is robust to missing/renamed columns and prints ONLY the two compact tables it builds:
#   (1) Outcome = operations   (2) Outcome = slave
# Each table shows the coefficient (SE) for both shift-share regressors under:
#   - Municipality + Year FE + baseline(2000)×t trends
#   - Municipality + UF×Year FE + baseline(2000)×t trends
########################################################################################################

rm(list = ls()); options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(fixest)
  library(xtable)
})

# --------------------------------------------------------------------------------
# PATHS
# --------------------------------------------------------------------------------
in_path <- "C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Data/slavery_guns_census_panel.rds"
out_dir <- "C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Track only files created in this run
created_files <- character(0)
save_tex <- function(xt, name) {
  path <- file.path(out_dir, name)
  print(xt, file = path, floating = FALSE, include.rownames = FALSE, type = "latex")
  created_files <<- c(created_files, path)
}
save_csv <- function(df, name) {
  path <- file.path(out_dir, name)
  write.csv(df, path, row.names = FALSE)
  created_files <<- c(created_files, path)
}

# --------------------------------------------------------------------------------
# Load
# --------------------------------------------------------------------------------
DT <- as.data.table(readRDS(in_path))
stopifnot(all(c("ibge_cod","year") %in% names(DT)))
DT[, ibge_cod := as.character(ibge_cod)]
DT[, year := as.integer(year)]
if (!"state_code" %in% names(DT)) DT[, state_code := as.integer(substr(ibge_cod, 1, 2))]

# --------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------
first_present <- function(candidates, in_names) {
  x <- intersect(candidates, in_names)
  if (length(x)) x[1] else NA_character_
}
pick_baseline_year <- function(x_years, target = 2000L) {
  x <- sort(unique(x_years[!is.na(x_years)]))
  if (!length(x)) return(NA_integer_)
  if (target %in% x) return(target)
  x[which.min(abs(x - target))]
}

# --------------------------------------------------------------------------------
# Baseline year per municipality (target 2000; else closest)
# --------------------------------------------------------------------------------
by_base <- DT[, .(base_year = pick_baseline_year(year, 2000L)), by = ibge_cod]
DT <- merge(DT, by_base, by = "ibge_cod", all.x = TRUE, allow.cartesian = TRUE)

# --------------------------------------------------------------------------------
# Locate likely columns present in this dataset (robust to naming)
# --------------------------------------------------------------------------------
nm <- names(DT)

# Population (level)
pop_col <- first_present(c("population","mean_population","pop_total","pop","mean_population_interp"), nm)

# GDP per capita (direct) OR GDP level + Population
gdppc_col <- first_present(c("gdppc","gdp_pc","pib_pc","pib_per_capita","gdp_per_capita"), nm)
gdp_col   <- first_present(c("gdp","pib","pib_mun","mean_pib_mun_interp"), nm)

# Rural share (direct) OR rural pop / total pop
rurshare_col  <- first_present(c("rural_share_2000","rural_share","share_rural","pop_rural_share"), nm)
rural_pop_col <- first_present(c("pop_rural","rural_pop","pop_rural_2000"), nm)

# Agricultural intensity proxy: harvested area per cap (preferred) else cattle per cap
harv_col   <- first_present(c("harvested_area","mean_harvested_area_interp","agri_area","area_agri"), nm)
cattle_col <- first_present(c("number_cattle","mean_number_cattle_interp","cattle"), nm)

# --------------------------------------------------------------------------------
# Build municipality-specific baseline values (at base_year)
# --------------------------------------------------------------------------------
B <- DT[year == base_year,
        .(ibge_cod,
          rural_share_2000 = {
            if (!is.na(rurshare_col)) {
              as.numeric(get(rurshare_col))
            } else if (!is.na(rural_pop_col) && !is.na(pop_col)) {
              rp <- as.numeric(get(rural_pop_col)); tp <- as.numeric(get(pop_col))
              fifelse(tp > 0, rp / tp, NA_real_)
            } else NA_real_
          },
          pop_2000 = if (!is.na(pop_col)) as.numeric(get(pop_col)) else NA_real_,
          gdppc_2000 = {
            if (!is.na(gdppc_col)) {
              as.numeric(get(gdppc_col))
            } else if (!is.na(gdp_col) && !is.na(pop_col)) {
              g <- as.numeric(get(gdp_col)); p <- as.numeric(get(pop_col))
              fifelse(p > 0, g / p, NA_real_)
            } else NA_real_
          },
          agint_2000 = {
            if (!is.na(harv_col) && !is.na(pop_col)) {
              ha <- as.numeric(get(harv_col)); p <- as.numeric(get(pop_col))
              fifelse(p > 0, ha / p, NA_real_)
            } else if (!is.na(cattle_col) && !is.na(pop_col)) {
              c <- as.numeric(get(cattle_col)); p <- as.numeric(get(pop_col))
              fifelse(p > 0, c / p, NA_real_)
            } else NA_real_
          }
        )
]

DT <- merge(DT, B, by = "ibge_cod", all.x = TRUE, allow.cartesian = TRUE)

# --------------------------------------------------------------------------------
# Construct municipality-specific baseline trends
# --------------------------------------------------------------------------------
DT[, t := year - 2000L]
DT[, tr_rural := rural_share_2000 * t]
DT[, tr_pop   := log1p(pop_2000)   * t]
DT[, tr_gdppc := log1p(gdppc_2000) * t]
DT[, tr_agint := agint_2000        * t]

trend_vars <- c("tr_rural","tr_pop","tr_gdppc","tr_agint")
trend_vars <- trend_vars[trend_vars %in% names(DT)]
trend_vars <- trend_vars[vapply(trend_vars, function(v) any(is.finite(DT[[v]]), na.rm = TRUE), logical(1))]

if (!length(trend_vars)) {
  stop("Baseline trend variables could not be constructed (all missing/constant). 
Check availability of rural share, population, GDP per capita, and agricultural intensity in the panel.")
}

cat("Baseline trend variables used:", paste(trend_vars, collapse = ", "), "\n")

# --------------------------------------------------------------------------------
# Outcomes & regressors (use only those present)
# --------------------------------------------------------------------------------
ys <- intersect(c("operations","slave"), names(DT))
xs <- intersect(c("ss_taurus","ss_imports"), names(DT))

if (!length(ys)) stop("No outcome available among: operations, slave.")
if (!length(xs)) stop("No shift-share regressor available among: ss_taurus, ss_imports.")

cat("Outcomes:", paste(ys, collapse = ", "), "\n")
cat("Regressors:", paste(xs, collapse = ", "), "\n")

# --------------------------------------------------------------------------------
# Estimation wrappers
#   Spec A: Municipality + Year FE + baseline(2000)×t trends
#   Spec B: Municipality + UF×Year FE + baseline(2000)×t trends
# --------------------------------------------------------------------------------
mk_trend_rhs <- paste(trend_vars, collapse = " + ")

run_one <- function(y, x) {
  fA <- as.formula(paste0(y, " ~ ", x, " + ", mk_trend_rhs, " | ibge_cod + year"))
  fB <- as.formula(paste0(y, " ~ ", x, " + ", mk_trend_rhs, " | ibge_cod + state_code^year"))
  mA <- feols(fA, data = DT, cluster = ~ ibge_cod)
  mB <- feols(fB, data = DT, cluster = ~ ibge_cod)
  list(A = mA, B = mB)
}

fmt_cell <- function(mod, x) {
  ct <- summary(mod)$coeftable
  if (x %in% rownames(ct)) {
    est <- unname(ct[x, "Estimate"])
    se  <- unname(ct[x, "Std. Error"])
    return(sprintf("%.4f (%.4f)", est, se))
  } else {
    return("")
  }
}

# --------------------------------------------------------------------------------
# Build TWO compact tables: one per outcome, columns = both regressors
# --------------------------------------------------------------------------------
for (y in ys) {
  # Fit once per regressor
  fits <- lapply(xs, function(x) run_one(y, x))
  names(fits) <- xs
  
  # Assemble rows by spec
  row_specA <- c("Municipality + Year FE + baseline(2000)×t",
                 vapply(xs, function(x) fmt_cell(fits[[x]][["A"]], x), character(1)))
  row_specB <- c("Municipality + UF×Year FE + baseline(2000)×t",
                 vapply(xs, function(x) fmt_cell(fits[[x]][["B"]], x), character(1)))
  
  # Make data.frame; dynamic column names based on available regressors
  col_names <- c("Specification", paste0(xs, " (Coef SE)"))
  out_df <- as.data.frame(rbind(row_specA, row_specB), stringsAsFactors = FALSE)
  names(out_df) <- col_names
  
  # Save LaTeX + CSV (one table per outcome)
  cap  <- sprintf("Effect of shift-share instruments on %s with municipality-specific baseline(2000) trends", y)
  lab  <- paste0("tab:", y, "_baseline_trends_compact")
  digits_vec <- rep(0, ncol(out_df) + 1)  # xtable requires length = ncol + 1
  xt   <- xtable(out_df, caption = cap, label = lab, digits = digits_vec)
  stub <- paste0("revresp_", y, "_baseline_trends_compact")
  save_tex(xt, paste0(stub, ".tex"))
  save_csv(out_df, paste0(stub, ".csv"))
}

# --------------------------------------------------------------------------------
# Print only files created in this run
# --------------------------------------------------------------------------------
created_files <- unique(created_files[file.exists(created_files)])
cat("\nFiles created in this run:\n")
print(created_files)
