# ====================================================================
# Reproducible Pipeline to Produce:
#   1) execution_ratios_trend.png  (2019–2025 budget execution ratios)
#   2) national_trend_appointments.png (2000–2019 total AFT appointments)
# --------------------------------------------------------------------
# Inputs expected:
#   - slavery_guns_budget_2000_2025.txt  (budget for forced labor actions)
#   - AFT_99_a_2025_por_UF.csv           (AFT appointments by UF and year)
# Outputs:
#   - execution_ratios_trend.png
#   - national_trend_appointments.png
# ====================================================================

## 0) Install/load packages ----------------------------------------------------
need <- c("data.table","dplyr","tidyr","ggplot2","scales","stringi","readr","grid")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(need, require, character.only = TRUE))
fixest::setFixest_nthreads(1)

## 1) Paths --------------------------------------------------------------------
budget_file <- "C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Data/slavery_guns_budget_2000_2025.txt"
aft_file    <- "C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/New_Data/LAI/Dados/AFT_99_a_2025_por_UF.csv"
out_dir     <- "C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/WD_replication_files/Results"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## 2) Helpers ------------------------------------------------------------------
convert_br_num <- function(x) {
  x <- gsub("\\.", "", x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

norm_names <- function(x) {
  x2 <- stringi::stri_trans_general(x, "Latin-ASCII")
  x2 <- tolower(x2)
  x2 <- gsub("[^a-z0-9]+", "_", x2)
  x2
}

find_col <- function(patterns, nm_vec) {
  for (p in patterns) {
    m <- grep(p, nm_vec, perl = TRUE)
    if (length(m)) return(nm_vec[m[1]])
  }
  return(NA_character_)
}

## 3) Figure 1: execution_ratios_trend.png -------------------------------------
orc <- tryCatch(fread(budget_file, sep = ",", header = TRUE),
                error = function(e) fread(budget_file, sep = "\t", header = TRUE))

setnames(orc, old = names(orc), new = norm_names(names(orc)))
nm <- names(orc)

col_year  <- find_col(c("^(ano|year)$"), nm)
col_act   <- find_col(c("^(acao|acao_.*|action)$"), nm)
col_alloc <- find_col(c("(dotacao_atualizada|dotacao|updated_allocation|^dot$)"), nm)
col_emp   <- find_col(c("^(empenhado|committed)$"), nm)
col_liq   <- find_col(c("^(liquidado|liquidated)$"), nm)
col_pago  <- find_col(c("^(pago|paid|^desp$)"), nm)

req <- c(col_year,col_act,col_alloc,col_emp,col_liq,col_pago)
if (anyNA(req) || any(!nzchar(req))) {
  stop("Could not auto-detect all required columns. Found names:\n",
       paste(nm, collapse = ", "))
}

orc <- orc |>
  dplyr::select(
    Year = all_of(col_year),
    Action = all_of(col_act),
    Updated_Allocation = all_of(col_alloc),
    Committed = all_of(col_emp),
    Liquidated = all_of(col_liq),
    Paid = all_of(col_pago)
  )

orc$Year <- as.integer(orc$Year)
for (v in c("Updated_Allocation","Committed","Liquidated","Paid")) {
  if (is.character(orc[[v]])) orc[[v]] <- convert_br_num(orc[[v]])
}

act_norm <- norm_names(orc$Action)
orc$Action <- ifelse(grepl("erradicacao.*escr|combate_ao_trabalho_escravo|slave", act_norm, ignore.case = TRUE),
                     "Inspection to Eradicate Slave Labor", orc$Action)
orc_sl <- dplyr::filter(orc, Action == "Inspection to Eradicate Slave Labor")

orc_plot <- orc_sl |>
  dplyr::select(Year, Updated_Allocation, Committed, Liquidated, Paid) |>
  tidyr::complete(Year = 2019:2025) |>
  dplyr::mutate(
    dplyr::across(c(Updated_Allocation,Committed,Liquidated,Paid), ~tidyr::replace_na(.x, 0)),
    committed_share  = dplyr::if_else(Updated_Allocation > 0, Committed  / Updated_Allocation, NA_real_),
    liquidated_share = dplyr::if_else(Updated_Allocation > 0, Liquidated / Updated_Allocation, NA_real_),
    paid_share       = dplyr::if_else(Updated_Allocation > 0, Paid       / Updated_Allocation, NA_real_)
  )

orc_long <- orc_plot |>
  dplyr::select(Year, committed_share, liquidated_share, paid_share) |>
  tidyr::pivot_longer(cols = -Year, names_to = "Type", values_to = "Share") |>
  dplyr::mutate(
    Type = factor(Type,
                  levels = c("committed_share","liquidated_share","paid_share"),
                  labels = c("Committed","Liquidated","Paid"))
  )

p_exec <- ggplot(orc_long, aes(x = Year, y = Share, color = Type)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1)) +
  labs(x = "Year", y = "Execution share (of updated allocation)", color = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.position   = "top",
    legend.key.width  = unit(1.6, "lines"),
    panel.border      = element_rect(colour = "grey60", fill = NA, linewidth = 0.6)
  )

ggsave(file.path(out_dir, "execution_ratios_trend.png"),
       p_exec, width = 8, height = 5, bg = "white", dpi = 300)

## 4) Figure 2: national_trend_appointments.png --------------------------------
aft <- readr::read_csv(aft_file, show_col_types = FALSE)
names(aft) <- norm_names(names(aft))

if (!all(c("ano","total_geral") %in% names(aft))) {
  stop("Could not find columns for year ('ano') and national total ('total_geral'). Found: ",
       paste(names(aft), collapse = ", "))
}

df <- aft |>
  dplyr::filter(!is.na(ano)) |>
  dplyr::mutate(ano = as.integer(ano)) |>
  dplyr::filter(ano >= 2000, ano <= 2019)

p_nat <- ggplot(df, aes(x = ano, y = total_geral)) +
  geom_line(linewidth = 1.1, color = "#2c7fb8") +
  geom_point(size = 2, color = "#2c7fb8") +
  geom_vline(xintercept = 2019, linetype = "dashed", color = "grey60") +
  labs(x = "Year", y = "Total appointments (national)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(out_dir, "national_trend_appointments.png"),
       p_nat, width = 8, height = 5, bg = "white", dpi = 300)

message("✅ Done. PNGs saved to: ", normalizePath(out_dir))
