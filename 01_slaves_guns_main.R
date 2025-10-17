#####################################################################################################################################
############################### CODE ################################################################################################
#####################################################################################################################################

rm(list = ls())

setwd("C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/Resultados")

# Install packages if they are not already installed
required_packages <- c("haven", "fixest", "dplyr", "remotes", 
                       "ritest", "ggplot2", "broom", "patchwork", 
                       "tidyverse", "Hmisc", "psych", "xtable", 
                       "broom.mixed", "coefplot", "gridExtra", 
                       "purrr", "conleyreg", "etwfe", "cowplot")

# Check which packages are not installed
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

# Install any missing packages
if(length(new_packages)) {
  install.packages(new_packages)
}

# Install ritest from GitHub if it's not installed
if (!requireNamespace("ritest", quietly = TRUE)) {
  remotes::install_github("grantmcdermott/ritest")
}

# Load all packages
lapply(required_packages, library, character.only = TRUE)

library(haven)
library(fixest)
library(dplyr)
# install.packages("remotes")
# remotes::install_github("grantmcdermott/ritest")
library(ritest)
library(ggplot2)
library(broom)
library(patchwork)
library(tidyverse)
library(Hmisc)
library(psych)
library(xtable)
library(broom.mixed)
library(coefplot)
library(gridExtra)
library(purrr)
library(dplyr)
#install.packages('conleyreg')
library(conleyreg)
library(etwfe)
library(cowplot)

### Cepaluni Data

df_reg <- readRDS("C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/Datasets/dfGunsPaperNovember24.rds") %>% 
  select(ibge_cod, year, homicidios, slave, operations, taurusRevenue, values_imp, 
        mean_uf_cod, mean_harvested_area_interp, mean_pib_mun_interp, mean_population_interp,
        mean_number_cattle_interp, mean_revenue_interp,
        pop_rural_2000, latitude, longitude) %>% 
  filter(mean_revenue_interp > 1000)  %>%
  na.omit()

df_clean <- df_reg %>%
  mutate(mean_uf_cod = as.numeric(substr(ibge_cod, 1, 2))) %>%
  filter(mean_uf_cod != 28, year >= 2000) %>%
  group_by(year) %>%
  mutate(rel_hom = homicidios / sum(homicidios, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(ibge_cod) %>%
  mutate(
    exposure = mean(rel_hom, na.rm = TRUE),
    imports_scaled = values_imp / max(values_imp, na.rm = TRUE),
    taurusRevenue_scaled = taurusRevenue / max(taurusRevenue, na.rm = TRUE),
    ss_imports = exposure * imports_scaled,
    ss_taurus = exposure * taurusRevenue_scaled,
    gdp_per_capita = mean_pib_mun_interp / mean_population_interp
  ) %>%
  ungroup()

check <- df_clean %>% 
  group_by(ibge_cod, year) %>% 
  filter(n()>1) %>%
  arrange(ibge_cod) %>% 
  ungroup()
check

names(df_clean)

# Check number of unique municipalities
num_unique_ibgecode <- df_clean %>%
  summarise(unique_ibgecodes = n_distinct(ibge_cod))

print(num_unique_ibgecode)

# Check unique municipalities for each year
unique_municipalities_per_year <- df_clean %>%
  group_by(year) %>%
  summarise(unique_municipalities = n_distinct(ibge_cod))

print(unique_municipalities_per_year)

# Check unique municipalities for each year
duplicates_check <- df_clean %>%
  group_by(year, ibge_cod) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

if(nrow(duplicates_check) > 0) {
  print("There are duplicates")
  print(duplicates_check)
} else {
  print("Each municipality is unique for each year")
}

zero_homicides_count <- sum(df_clean$exposure == 0)

#####################################################################################################################################
############################### Save the dataset ####################################################################################
#####################################################################################################################################

# Set your working directory to where you want to save the file
# It's good practice to use setwd() to ensure your file paths are managed correctly
# setwd("C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/Datasets")

# Now, save the dataframe to a .csv file in the specified directory
# write.csv(df_clean, "shift_shareMarch2024.csv", row.names = FALSE)


#####################################################################################################################################
############################### Descriptive Stats ###################################################################################
#####################################################################################################################################

vars <- c('slave', 'operations', 'values_imp', 'taurusRevenue_scaled', 'homicidios',
          "mean_harvested_area_interp", "mean_pib_mun_interp","mean_population_interp", 
          'gdp_per_capita', "mean_number_cattle_interp", "mean_revenue_interp"
)

# Calculate descriptive statistics for the specified variables
descriptive_stats <- data.frame(
  Mean = sapply(df_clean[, vars], mean, na.rm = TRUE),
  Min = sapply(df_clean[, vars], min, na.rm = TRUE),
  Max = sapply(df_clean[, vars], max, na.rm = TRUE),
  Std_Dev = sapply(df_clean[, vars], sd, na.rm = TRUE),
  N_Obs = sapply(df_clean[, vars], function(x) sum(!is.na(x)))
)

# Rename the columns
descriptive_stats <- descriptive_stats %>%
  rename(
    Mean = Mean,
    Min = Min,
    Max = Max,
    Std.Dev. = Std_Dev,
    N.Obs. = N_Obs
  )

# Print the descriptive statistics table in LaTeX format
print(xtable(descriptive_stats, digits = 3, caption = "Descriptive Statistics", label = "tab:descriptive_stats"), floating = FALSE)

#####################################################################################################################################
############################### Descriptive Stats: Quantiles ########################################################################
#####################################################################################################################################

vars <- c('slave', 'operations', 'values_imp', 'taurusRevenue_scaled', 'homicidios',
          "mean_harvested_area_interp", "mean_population_interp", 
          'gdp_per_capita', "mean_number_cattle_interp", "mean_revenue_interp"
)

# Calculate descriptive statistics for the specified variables
descriptive_stats <- data.frame(
  Mean = sapply(df_clean[, vars], mean, na.rm = TRUE),
  Q1 = sapply(df_clean[, vars], quantile, probs = 0.25, na.rm = TRUE),
  Median = sapply(df_clean[, vars], median, na.rm = TRUE),
  Q3 = sapply(df_clean[, vars], quantile, probs = 0.75, na.rm = TRUE),
  Min = sapply(df_clean[, vars], min, na.rm = TRUE),
  Max = sapply(df_clean[, vars], max, na.rm = TRUE),
  Std_Dev = sapply(df_clean[, vars], sd, na.rm = TRUE),
  N_Obs = sapply(df_clean[, vars], function(x) sum(!is.na(x)))
)

# You don't need to rename the columns if they are already correctly named
# The rename operation was redundant in the previous code

# Print the descriptive statistics table in LaTeX format
print(xtable(descriptive_stats, digits = 3, caption = "Descriptive Statistics", label = "tab:descriptive_stats"), floating = FALSE)

#####################################################################################################################################
############################### Main Regressions ####################################################################################
#####################################################################################################################################

# Central Results - Operations and Slavery Significant and in the Right Direction
reg1 <- feols(operations ~ ss_imports | ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg1)

reg2 <- feols(operations ~ ss_taurus | ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg2)

reg3 <- feols(slave ~ ss_imports | ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg3)

reg4 <- feols(slave ~ ss_taurus | ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg4)

# Regressions with Controls
# Central Results - Operations and Slavery Significant and in the Right Direction
reg5 <- feols(operations ~ ss_imports + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log(mean_revenue_interp+1) | ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg5)

reg6 <- feols(operations ~ ss_taurus + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log( mean_revenue_interp+1)| ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg6)

reg7 <- feols(slave ~ ss_imports + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log(mean_revenue_interp+1)| ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg7)

reg8 <- feols(slave ~ ss_taurus + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log( mean_revenue_interp+1)| ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg8)

# Initial Robustness - Unit-Specific Time Trends - Results are robust
reg9 <- feols(operations ~ ss_imports + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log( mean_revenue_interp+1)| ibge_cod+year + mean_uf_cod[year], panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg9)

reg10 <- feols(operations ~ ss_taurus + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log( mean_revenue_interp+1) | ibge_cod+year + mean_uf_cod[year], panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
               data= df_clean)
summary(reg10)

reg11 <- feols(slave ~ ss_imports + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log( mean_revenue_interp+1)| ibge_cod+year + mean_uf_cod[year], panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
               data= df_clean)
summary(reg11)

reg12 <- feols(slave ~ ss_taurus + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log( mean_revenue_interp+1) | ibge_cod+year + mean_uf_cod[year], panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
               data= df_clean)
summary(reg12)

# Calculating Impact on Audits
reg1$coefficients[1]*(df_clean %>% select(ibge_cod, exposure) %>% distinct %>% ungroup %>% select(exposure) %>% as.matrix %>% quantile(0.5, na.rm = T)) / (df_clean %>% select(ibge_cod, operations) %>% ungroup %>% select(operations) %>% as.matrix %>% mean(na.rm = T))
reg2$coefficients[1]*(df_clean %>% select(ibge_cod, exposure) %>% distinct %>% ungroup %>% select(exposure) %>% as.matrix %>% quantile(0.5, na.rm = T)) / (df_clean %>% select(ibge_cod, operations) %>% ungroup %>% select(operations) %>% as.matrix %>% mean(na.rm = T))
# Calculating Impacts on Slavery
reg3$coefficients[1]*(df_clean %>% select(ibge_cod, exposure) %>% distinct %>% ungroup %>% select(exposure) %>% as.matrix %>% quantile(0.5, na.rm = T)) / (df_clean %>% select(ibge_cod, slave) %>% ungroup %>% select(slave) %>% as.matrix %>% mean(na.rm = T))
reg4$coefficients[1]*(df_clean %>% select(ibge_cod, exposure) %>% distinct %>% ungroup %>% select(exposure) %>% as.matrix %>% quantile(0.5, na.rm = T)) / (df_clean %>% select(ibge_cod, slave) %>% ungroup %>% select(slave) %>% as.matrix %>% mean(na.rm = T))

#####################################################################################################################################
############################### Coef. Plots #########################################################################################
#####################################################################################################################################

# Define a function to create the coefficient plot for a given regression model
create_coeff_plot <- function(reg_model, model_name) {
  coef_df <- try(tidy(reg_model), silent = TRUE)
  
  if (inherits(coef_df, "try-error")) {
    # If tidy() fails, create an empty plot
    return(ggplot() + ggtitle(paste("Error in", model_name)) + theme_void())
  }
  
  coef_df <- coef_df %>%
    mutate(term = case_when(
      term == "ss_imports" ~ "Imports",
      term == "ss_taurus" ~ "Taurus",
      TRUE ~ term
    )) %>%
    filter(term %in% c("Imports", "Taurus"))
  
  if (nrow(coef_df) == 0) {
    # If no relevant coefficients, create an empty plot
    return(ggplot() + ggtitle(paste("No coefficients in", model_name)) + theme_void())
  }
  
  p <- ggplot(coef_df, aes(x = term, y = estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = " ", y = "Intention to treat") +
    ggtitle(model_name) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# List of regression models and their names
reg_models <- list(reg1, reg2, reg3, reg4, reg5, reg6, reg7, reg8, reg9, reg10, reg11, reg12)
model_names <- c("Audits (Bivariate)", "Audits (Bivariate)", "Slaves (Bivariate)", "Slaves (Bivariate)", 
                 "Audits (C)", "Audits (C)", "Slaves (C)", "Slaves (C)", 
                 "Audits (C + SSTT)", "Audits (C + SSTT)", "Slaves (C + SSTT)", "Slaves (C + SSTT)")

# Create a list of coefficient plots for each model
coef_plots <- lapply(seq_along(reg_models), function(i) {
  create_coeff_plot(reg_models[[i]], model_names[i])
})

# Combine the plots using cowplot by converting the list into individual grobs
combined_plots <- plot_grid(plotlist = coef_plots, ncol = 4, nrow = ceiling(length(coef_plots) / 4))

print(combined_plots)

# Save the combined plots
ggsave("main_results.pdf", combined_plots, width = 12, height = 10)  # Adjust width and height as needed

#####################################################################################################################################
############################### Common Trend Plots: Shift-share ###################################################################
#####################################################################################################################################

# Calculate the mean of bartik_taurus
mean_bartik_taurus <- mean(df_clean$exposure, na.rm = TRUE)

# Create new dummy variables for above and below the mean bartik_taurus
df_clean <- df_clean %>%
  mutate(dummy_above_mean = ifelse(exposure > mean_bartik_taurus, 1, 0),
         dummy_below_mean = ifelse(exposure <= mean_bartik_taurus, 1, 0))

# Group by the new dummies in separate datasets
collapsed_data_above_mean <- df_clean %>%
  filter(dummy_above_mean == 1) %>%
  group_by(year) %>%
  summarise(mean_audits_above = mean(operations, na.rm = TRUE),
            mean_slave_above = mean(slave, na.rm = TRUE), .groups = "drop")

collapsed_data_below_mean <- df_clean %>%
  filter(dummy_below_mean == 1) %>%
  group_by(year) %>%
  summarise(mean_audits_below = mean(operations, na.rm = TRUE),
            mean_slave_below = mean(slave, na.rm = TRUE), .groups = "drop")

# Join the datasets to make plotting easier
collapsed_data_bartik_taurus <- full_join(collapsed_data_above_mean, collapsed_data_below_mean, by = "year")

# Adjust your plots to use the new dataset.
plot_3 <- ggplot() +
  geom_line(data = collapsed_data_bartik_taurus, aes(x = year, y = mean_audits_above, color = "Exposure Values Above the Mean")) +
  geom_line(data = collapsed_data_bartik_taurus, aes(x = year, y = mean_audits_below, color = "Exposure Values Below the Mean")) +
  labs(x = "Year", y = "Mean Audits", color = "Shift-share Taurus") +
  scale_x_continuous(breaks = c(2001:2019)) +
  theme_minimal()

plot_4 <- ggplot() +
  geom_line(data = collapsed_data_bartik_taurus, aes(x = year, y = mean_slave_above, color = "Exposure Values Above the Mean")) +
  geom_line(data = collapsed_data_bartik_taurus, aes(x = year, y = mean_slave_below, color = "Exposure Values Below the Mean")) +
  labs(x = "Year", y = "Mean Slaves", color = "Shift-share Taurus") +
  scale_x_continuous(breaks = c(2001:2019)) +
  theme_minimal()

# Print or save the plots
print(plot_3)
print(plot_4)

# Combine and save the plot using cowplot
combined_plots2 <- cowplot::plot_grid(plot_3, plot_4, ncol = 1)
print(combined_plots2)

ggsave("bartik_common_trends_plots.pdf", combined_plots2, width = 10, height = 10)

#####################################################################################################################################
############################### Randomization Inference #############################################################################
#####################################################################################################################################

est_ri_reg1 <- ritest(reg1, "ss_imports", cluster='ibge_cod', reps=500, seed=4321)
est_ri_reg1
# Show the plot
plot(est_ri_reg1)

est_ri_reg2 <- ritest(reg2, "ss_taurus", cluster='ibge_cod', reps=500, seed=4321)
est_ri_reg2
# Show the plot
plot(est_ri_reg2)

# Run ritest for reg3 (bartik_imports)
est_ri_reg3 <- ritest(reg3, "ss_imports", cluster = 'ibge_cod', reps = 500, seed = 4321)
est_ri_reg3
# Show the plot
plot(est_ri_reg3)

# Run ritest for reg4 (bartik_taurus)
est_ri_reg4 <- ritest(reg4, "ss_taurus", cluster = 'ibge_cod', reps = 500, seed = 4321)
est_ri_reg4
# Show the plot
plot(est_ri_reg4)

# Run ritest for reg5 (bartik_imports)
est_ri_reg5 <- ritest(reg5, "ss_imports", cluster = 'ibge_cod', reps = 500, seed = 4321)
est_ri_reg5
# Show the plot
plot(est_ri_reg5)

# Run ritest for reg6 (bartik_taurus)
est_ri_reg6 <- ritest(reg6, "ss_taurus", cluster = 'ibge_cod', reps = 500, seed = 4321)
est_ri_reg6
# Show the plot
plot(est_ri_reg6)

# Run ritest for reg7 (bartik_imports)
est_ri_reg7 <- ritest(reg7, "ss_imports", cluster = 'ibge_cod', reps = 500, seed = 4321)
est_ri_reg7
# Show the plot
plot(est_ri_reg7)

# Run ritest for reg8 (bartik_taurus)
est_ri_reg8 <- ritest(reg8, "ss_taurus", cluster = 'ibge_cod', reps = 500, seed = 4321)
est_ri_reg8
# Show the plot
plot(est_ri_reg8)

# Run ritest for reg9 (bartik_imports)
est_ri_reg9 <- ritest(reg9, "ss_imports", cluster = 'ibge_cod', reps = 500, seed = 4321)
est_ri_reg9
# Show the plot
plot(est_ri_reg9)

# Run ritest for reg10 (bartik_taurus)
est_ri_reg10 <- ritest(reg10, "ss_taurus", cluster = 'ibge_cod', reps = 500, seed = 4321)
est_ri_reg10
# Show the plot
plot(est_ri_reg10)

# Run ritest for reg11 (bartik_imports)
est_ri_reg11 <- ritest(reg11, "ss_imports", cluster = 'ibge_cod', reps = 500, seed = 4321)
est_ri_reg11
# Show the plot
plot(est_ri_reg11)

# Run ritest for reg12 (bartik_taurus)
est_ri_reg12 <- ritest(reg12, "ss_taurus", cluster = 'ibge_cod', reps = 500, seed = 4321)
est_ri_reg12
# Show the plot
plot(est_ri_reg12)

# Set the layout for the combined plot (2 columns)
par(mfrow = c(6, 2)) # 6 rows and 2 columns

# Set the margins to adjust the plotting area
# The default margin values are c(bottom, left, top, right) = c(5.1, 4.1, 4.1, 2.1)
# You can modify the values based on your preference
par(mar = c(1, 1, 1, 1)) # Increased the bottom margin to add space between figures

# Plot for est_ri_reg1
plot(est_ri_reg1)

# Plot for est_ri_reg2
plot(est_ri_reg2)

# Plot for est_ri_reg3
plot(est_ri_reg3)

# Plot for est_ri_reg4
plot(est_ri_reg4)

# Plot for est_ri_reg5
plot(est_ri_reg5)

# Plot for est_ri_reg6
plot(est_ri_reg6)

# Plot for est_ri_reg7
plot(est_ri_reg7)

# Plot for est_ri_reg8
plot(est_ri_reg8)

# Plot for est_ri_reg9
plot(est_ri_reg9)

# Plot for est_ri_reg10
plot(est_ri_reg10)

# Plot for est_ri_reg11
plot(est_ri_reg11)

# Plot for est_ri_reg12
plot(est_ri_reg12)

# Assuming you have already run est_ri_reg1 to est_ri_reg12 as mentioned earlier
# Open a PDF device to save the plots
png("ritest_plots.png", width = 800, height = 1200) # Adjust width and height as needed

# Close the PDF device
dev.off()

#####################################################################################################################################
############################### Rural Population: Audits ############################################################################
#####################################################################################################################################

# Calculate deciles of rural population
rural_pop_deciles <- quantile(df_clean$pop_rural_2000, probs = seq(0, 1, by = 0.1))
rural_pop_deciles

# Filter the data frame to keep only rows where pop_rural_2000 <= 12477
df_filtered_drop_high_10 <- df_clean %>%
  filter(pop_rural_2000 <= 12558)

# Run the regression on the filtered dataset
reg_rob1 <- feols(operations ~ ss_imports + log(mean_harvested_area_interp+1) +
                    log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +
                    log(mean_revenue_interp+1) |
                    ibge_cod + year,
                  panel.id = ~ibge_cod + year,
                  cluster = ~ibge_cod,
                  data = df_filtered_drop_high_10)

# Summary of the regression
summary(reg_rob1)

# Filter the data frame to keep only rows where pop_rural_2000 <= 8305
df_filtered_drop_high_20 <- df_clean %>%
  filter(pop_rural_2000 <= 8384)

# Run the regression on the filtered dataset
reg_rob2 <- feols(operations ~ ss_imports + log(mean_harvested_area_interp+1) +
                    log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +
                    log(mean_revenue_interp+1) |
                    ibge_cod + year,
                  panel.id = ~ibge_cod + year,
                  cluster = ~ibge_cod,
                  data = df_filtered_drop_high_20)

# Summary of the regression
summary(reg_rob2)

# Run the regression on the filtered dataset
reg_rob3 <- feols(slave ~ ss_imports + log(mean_harvested_area_interp+1) +
                    log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +
                    log(mean_revenue_interp+1) |
                    ibge_cod + year,
                  panel.id = ~ibge_cod + year,
                  cluster = ~ibge_cod,
                  data = df_filtered_drop_high_10)

# Summary of the regression
summary(reg_rob3)

# Run the regression on the filtered dataset
reg_rob4 <- feols(slave ~ ss_imports + log(mean_harvested_area_interp+1) +
                    log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +
                    log(mean_revenue_interp+1) |
                    ibge_cod + year,
                  panel.id = ~ibge_cod + year,
                  cluster = ~ibge_cod,
                  data = df_filtered_drop_high_20)

# Summary of the regression
summary(reg_rob4)


#####################################################################################################################################
############################### Coefficient Plot ####################################################################################
#####################################################################################################################################

# Define a function to create coefficient plots
create_coeff_plot <- function(reg_model, model_name) {
  coef_df <- tidy(reg_model) %>%
    filter(term %in% c("ss_imports",
                       "log(mean_harvested_area_interp+1)",
                       "log(gdp_per_capita+1)",
                       "log(mean_number_cattle_interp+1)",
                       "log(mean_revenue_interp+1)")) %>%
    mutate(term = ifelse(term == "ss_imports", "Shift-share Imports", term))  # Rename "bartik_imports" to "Shift-share 1"
  
  p <- ggplot(coef_df, aes(x = term, y = estimate)) +
    geom_point(shape = 16, size = 2.777, color = "red", fill = "blue") +  # Custom point aesthetics
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(aes(ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error), width = 0.2, size = 1, color = "black") +  
    labs(x = " ", y = "Intention to treat") +  # Fixed missing quote
    ggtitle(model_name) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# Correctly define model names with swapped terms
model_names <- c("Drop 10% with high 2000 rural population (Audits)", "Drop 20% with high 2000 rural population (Audits)",
                 "Drop 10% with high 2000 rural population (Slaves)", "Drop 20% with high 2000 rural population (Slaves)")

# List of regression models
reg_models <- list(reg_rob1, reg_rob2, reg_rob3, reg_rob4)

# Create coefficient plots
coef_plots <- lapply(seq_along(reg_models), function(i) {
  create_coeff_plot(reg_models[[i]], model_names[i])
})

# Combine the plots using cowplot
combined_plots <- cowplot::plot_grid(plotlist = coef_plots, ncol = 2)

# Print the combined plots
print(combined_plots)

# Save the combined plots as a PDF
ggsave("rural_pop_coef_plot.pdf", combined_plots, width = 10, height = 7)

#####################################################################################################################################
############################### Shift-share (Imports and Taurus): Conley std. errors ################################################
#####################################################################################################################################
#####################################################################################################################################
############################### Coef. Plots #########################################################################################
#####################################################################################################################################

# Define the function to create coefficient plots
create_coeff_plot <- function(model, title) {
  # Tidy the model summary to get coefficients and std. error for the Shift-share variables
  coef_df <- tidy(model) %>%
    mutate(term = ifelse(term == "ss_imports", "Shift-share (Imports)",
                         ifelse(term == "ss_taurus", "Shift-share (Taurus)", term))) %>%
    filter(term %in% c("Shift-share (Imports)", "Shift-share (Taurus)"))
  
  # Create the plot
  p <- ggplot(coef_df, aes(x = term, y = estimate)) +
    geom_errorbar(aes(ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error),
                  width = 0.3, color = "black") +  # Changed width to 0.3 (you can adjust this value)
    geom_point(size = 3, shape = 16, color = "blue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = title, x = NULL, y = "Intention to treat") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}


# regression 1: without additional controls
reg_conley1 <- feols(slave ~ ss_imports + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log(mean_revenue_interp+1)| ibge_cod+year, panel.id = ~ibge_cod+year, vcov = conley(100),
                     data= df_clean)
reg_conley1

# regression 2: with additional controls
reg_conley2 <- feols(slave ~ ss_imports + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log(mean_revenue_interp+1)| ibge_cod+year, panel.id = ~ibge_cod+year, vcov = conley(100),
                     data= df_clean)
reg_conley2

# regression 3: Different dependent variable, without additional controls
reg_conley3 <- feols(operations ~ ss_imports + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log(mean_revenue_interp+1)| ibge_cod+year, panel.id = ~ibge_cod+year, vcov = conley(100),
                     data= df_clean)
reg_conley3

# regression 4: Different dependent variable, with additional controls
reg_conley4 <- feols(operations ~ ss_imports + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log(mean_revenue_interp+1)| ibge_cod+year, panel.id = ~ibge_cod+year, vcov = conley(100),
                     data= df_clean)
reg_conley4

# regression 5: without additional controls
reg_conley5 <- feols(slave ~ ss_taurus + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log(mean_revenue_interp+1)| ibge_cod+year, panel.id = ~ibge_cod+year, vcov = conley(100),
                     data= df_clean)
reg_conley5

# regression 6: with additional controls
reg_conley6 <- feols(slave ~ ss_taurus + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log(mean_revenue_interp+1)| ibge_cod+year, panel.id = ~ibge_cod+year, vcov = conley(100),
                     data= df_clean)
reg_conley6

# regression 7: Different dependent variable, without additional controls
reg_conley7 <- feols(operations ~ ss_taurus + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log(mean_revenue_interp+1)| ibge_cod+year, panel.id = ~ibge_cod+year, vcov = conley(100),
                     data= df_clean)
reg_conley7

# regression 8: Different dependent variable, with additional controls
reg_conley8 <- feols(operations ~ ss_taurus + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log(mean_revenue_interp+1)| ibge_cod+year, panel.id = ~ibge_cod+year, vcov = conley(100),
                     data= df_clean)
reg_conley8

# Create coefficient plots for each regression
plot1 <- create_coeff_plot(reg_conley1, "Audits (Bivariate)")
plot2 <- create_coeff_plot(reg_conley2, "Audits (Controls)")
plot3 <- create_coeff_plot(reg_conley3, "Slaves (Bivariate)")
plot4 <- create_coeff_plot(reg_conley4, "Slaves (Controls)")
plot5 <- create_coeff_plot(reg_conley5, "Audits (Bivariate)")
plot6 <- create_coeff_plot(reg_conley6, "Audits (Controls)")
plot7 <- create_coeff_plot(reg_conley7, "Slaves (Bivariate)")
plot8 <- create_coeff_plot(reg_conley8, "Slaves (Controls)")

# Combine the plots using cowplot, arranging them in a 2x4 grid
combined_plot <- cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, ncol = 2)

# Display the combined plot
print(combined_plot)

# Save the combined plot to a file
ggsave("shift_share_conley_cowplot.pdf", combined_plot, width = 7, height = 11)

#####################################################################################################################################
############################### Main regressions (Poisson) ##########################################################################
#####################################################################################################################################

# Central Results - Slavery and Operations Significant and in the Right Direction
reg_poisson1 <- fepois(slave ~ ss_imports | ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg_poisson1)

reg_poisson2 <- fepois(slave ~ ss_taurus | ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg_poisson2)

reg_poisson3 <- fepois(operations ~ ss_imports | ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg_poisson3)

reg_poisson4 <- fepois(operations ~ ss_taurus | ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg_poisson4)

# reg_poissonressions with Controls
# Central Results - Slavery and Operations Significant and in the Right Direction
reg_poisson5 <- fepois(slave ~ ss_imports + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log(mean_revenue_interp+1) | ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg_poisson5)

reg_poisson6 <- fepois(slave ~ ss_taurus + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log( mean_revenue_interp+1)| ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg_poisson6)

reg_poisson7 <- fepois(operations ~ ss_imports + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log(mean_revenue_interp+1)| ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg_poisson7)

reg_poisson8 <- fepois(operations ~ ss_taurus + log(mean_harvested_area_interp+1) + log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +log( mean_revenue_interp+1)| ibge_cod+year, panel.id = ~ibge_cod+year, cluster = ~ibge_cod,
              data= df_clean)
summary(reg_poisson8)

#####################################################################################################################################
############################### Coef. Plots #########################################################################################
#####################################################################################################################################

# Corrected function to create the coefficient plot
create_coeff_plot <- function(reg_poisson_model, model_name) {
  coef_df <- tidy(reg_poisson_model) %>%
    mutate(term = ifelse(term == "ss_imports", "Imports", ifelse(term == "ss_taurus", "Taurus", term))) %>%
    filter(term %in% c("Imports", "Taurus"))
  
  p <- ggplot(coef_df, aes(x = term, y = estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = " ", y = "Intention to treat") + # Corrected: added y-axis label
    ggtitle(model_name) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# List of regression models and their names, replacing Audits by Slaves and Slaves by Audits
reg_poisson_models <- list(reg_poisson1, reg_poisson2, reg_poisson3, reg_poisson4, reg_poisson5, reg_poisson6, reg_poisson7, reg_poisson8) 
model_names <- c("Audits (Bivariate)", "Audits (Bivariate)", "Slaves (Bivariate)", "Slaves (Bivariate)", 
                 "Audits (Controls)", "Audits (Controls)", "Slaves (Controls)", "Slaves (Controls)")

# Create a list of coefficient plots for each model
coef_plots <- lapply(seq_along(reg_poisson_models), function(i) {
  create_coeff_plot(reg_poisson_models[[i]], model_names[i])
})

# Combine the plots using cowplot
combined_plots <- cowplot::plot_grid(plotlist = coef_plots, ncol = 2)

# Print the combined plots
print(combined_plots)

# Save the combined plots with specified dimensions
ggsave("main_results_poisson.pdf", combined_plots, width = 10, height = 10)


#####################################################################################################################################
############################### Homicidies Quantiles ################################################################################
#####################################################################################################################################

# Calculate the quantiles for homicidios
quantiles <- quantile(df_clean$homicidios, probs = c(0.25, 0.5, 0.75))

# Create a new variable in df_clean that categorizes based on the quantiles
df_clean$homicide_quantile <- cut(df_clean$homicidios, 
                                  breaks = c(-Inf, quantiles, Inf), 
                                  labels = c("Low", "Medium-Low", "Medium-High", "High"))

# Make sure the levels of 'homicide_quantile' are set in the correct order for plotting
df_clean$homicide_quantile <- factor(df_clean$homicide_quantile, levels = c("Low", "Medium-Low", "Medium-High", "High"))

# Run the regression for each quantile group
reg_low <- feols(operations ~ ss_imports + log(mean_harvested_area_interp+1) + 
                   log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +
                   log(mean_revenue_interp+1) | ibge_cod + year, 
                 panel.id = ~ibge_cod + year, 
                 cluster = ~ibge_cod, 
                 data = df_clean[df_clean$homicide_quantile == "Low", ])

reg_medium_low <- feols(operations ~ ss_imports + log(mean_harvested_area_interp+1) + 
                          log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +
                          log(mean_revenue_interp+1) | ibge_cod + year, 
                        panel.id = ~ibge_cod + year, 
                        cluster = ~ibge_cod, 
                        data = df_clean[df_clean$homicide_quantile == "Medium-Low", ])

reg_medium_high <- feols(operations ~ ss_imports + log(mean_harvested_area_interp+1) + 
                           log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +
                           log(mean_revenue_interp+1) | ibge_cod + year, 
                         panel.id = ~ibge_cod + year, 
                         cluster = ~ibge_cod, 
                         data = df_clean[df_clean$homicide_quantile == "Medium-High", ])

reg_high <- feols(operations ~ ss_imports + log(mean_harvested_area_interp+1) + 
                    log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +
                    log(mean_revenue_interp+1) | ibge_cod + year, 
                  panel.id = ~ibge_cod + year, 
                  cluster = ~ibge_cod, 
                  data = df_clean[df_clean$homicide_quantile == "High", ])

# Combine the regression results into a tidy data frame
tidy_low <- tidy(reg_low) %>% mutate(quantile = "Low")
tidy_medium_low <- tidy(reg_medium_low) %>% mutate(quantile = "Medium-Low")
tidy_medium_high <- tidy(reg_medium_high) %>% mutate(quantile = "Medium-High")
tidy_high <- tidy(reg_high) %>% mutate(quantile = "High")

# Combine all the tidy results into one data frame
tidy_all <- bind_rows(tidy_low, tidy_medium_low, tidy_medium_high, tidy_high)

# Filter out only the coefficient for bartik_imports and rename it
tidy_all <- tidy_all %>% 
  filter(term == "ss_imports") %>% 
  mutate(term = "Shift-share 1")

# Create the coefficient plot using ggplot, ensuring the correct order of quantiles
coef_plot <- ggplot(tidy_all, aes(x = factor(quantile, levels = c("Low", "Medium-Low", "Medium-High", "High")), 
                                  y = estimate, color = quantile)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) +  # Make points larger
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                position = position_dodge(width = 0.5), width = 0.25, size = 1) +  # Make error bars thicker
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.2) +  # Add a thick horizontal line at zero
  labs(title = "Coefficient Plot for Gun Imports's Impact on Audits by Homicide Quantiles", 
       x = "Homicide Quantiles", 
       y = "Estimate (Shift-share Imports)") +
  scale_x_discrete(limits = c("Low", "Medium-Low", "Medium-High", "High")) +  # Explicitly set the order
  theme_minimal()

# Save the plot as a PNG file
ggsave("coefficient_plot_shift_share.png", plot = coef_plot, width = 8, height = 6)

# Print the plot
print(coef_plot)

#####################################################################################################################################
########################### Homicidies Quantiles: Shift-share 2 #####################################################################
#####################################################################################################################################

# Calculate the quantiles for homicidios
quantiles <- quantile(df_clean$homicidios, probs = c(0.25, 0.5, 0.75))

# Create a new variable in df_clean that categorizes based on the quantiles
df_clean$homicide_quantile <- cut(df_clean$homicidios, 
                                  breaks = c(-Inf, quantiles, Inf), 
                                  labels = c("Low", "Medium-Low", "Medium-High", "High"))

# Make sure the levels of 'homicide_quantile' are set in the correct order for plotting
df_clean$homicide_quantile <- factor(df_clean$homicide_quantile, levels = c("Low", "Medium-Low", "Medium-High", "High"))

# Run the regression for each quantile group
reg_low <- feols(operations ~ ss_taurus + log(mean_harvested_area_interp+1) + 
                   log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +
                   log(mean_revenue_interp+1) | ibge_cod + year, 
                 panel.id = ~ibge_cod + year, 
                 cluster = ~ibge_cod, 
                 data = df_clean[df_clean$homicide_quantile == "Low", ])

reg_medium_low <- feols(operations ~ ss_taurus + log(mean_harvested_area_interp+1) + 
                          log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +
                          log(mean_revenue_interp+1) | ibge_cod + year, 
                        panel.id = ~ibge_cod + year, 
                        cluster = ~ibge_cod, 
                        data = df_clean[df_clean$homicide_quantile == "Medium-Low", ])

reg_medium_high <- feols(operations ~ ss_taurus + log(mean_harvested_area_interp+1) + 
                           log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +
                           log(mean_revenue_interp+1) | ibge_cod + year, 
                         panel.id = ~ibge_cod + year, 
                         cluster = ~ibge_cod, 
                         data = df_clean[df_clean$homicide_quantile == "Medium-High", ])

reg_high <- feols(operations ~ ss_taurus + log(mean_harvested_area_interp+1) + 
                    log(gdp_per_capita+1) + log(mean_number_cattle_interp+1) +
                    log(mean_revenue_interp+1) | ibge_cod + year, 
                  panel.id = ~ibge_cod + year, 
                  cluster = ~ibge_cod, 
                  data = df_clean[df_clean$homicide_quantile == "High", ])

# Combine the regression results into a tidy data frame
tidy_low <- tidy(reg_low) %>% mutate(quantile = "Low")
tidy_medium_low <- tidy(reg_medium_low) %>% mutate(quantile = "Medium-Low")
tidy_medium_high <- tidy(reg_medium_high) %>% mutate(quantile = "Medium-High")
tidy_high <- tidy(reg_high) %>% mutate(quantile = "High")

# Combine all the tidy results into one data frame
tidy_all <- bind_rows(tidy_low, tidy_medium_low, tidy_medium_high, tidy_high)

# Filter out only the coefficient for bartik_taurus and rename it
tidy_all <- tidy_all %>% 
  filter(term == "ss_taurus") %>% 
  mutate(term = "Shift-share 2")

# Create the coefficient plot using ggplot, ensuring the correct order of quantiles
coef_plot <- ggplot(tidy_all, aes(x = factor(quantile, levels = c("Low", "Medium-Low", "Medium-High", "High")), 
                                  y = estimate, color = quantile)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) +  # Make points larger
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), 
                position = position_dodge(width = 0.5), width = 0.25, size = 1) +  # Make error bars thicker
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.2) +  # Add a thick horizontal line at zero
  labs(title = "Coefficient Plot for Taurus Weapons' Impact on Audits by Homicide Quantiles", 
       x = "Homicide Quantiles", 
       y = "Estimate (Shift-share Taurus)") +
  scale_x_discrete(limits = c("Low", "Medium-Low", "Medium-High", "High")) +  # Explicitly set the order
  theme_minimal()

# Save the plot as a PNG file
ggsave("coefficient_plot_shift_share_2.png", plot = coef_plot, width = 8, height = 6)

# Print the plot
print(coef_plot)

#####################################################################################################################################
############################### Regions #############################################################################################
#####################################################################################################################################

# Get the unique values of the 'mean_uf_cod' column
unique_values <- unique(df_clean$mean_uf_cod)

# Print the unique values
print(unique_values)

# Assuming 'df_clean' is the name of your dataset, and 'mean_uf_cod' is the column with state codes
# Create a new variable 'region' representing the regions

# Define a function to map states to their respective regions
state_to_region <- function(state_code) {
  if (state_code %in% c(11, 12, 13, 14, 15, 16, 17)) {
    return("North")
  } else if (state_code %in% c(21, 22, 23, 24, 25, 26, 27, 29)) {
    return("Northeast")
  } else if (state_code %in% c(31, 32, 33, 35)) {
    return("Southeast")
  } else if (state_code %in% c(41, 42, 43)) {
    return("South")
  } else if (state_code %in% c(50, 51, 52, 53)) {
    return("Central-West")
  } else {
    return("Unknown") # In case there are any state codes not matching the defined regions
  }
}

# Use the state_to_region function to create the 'region' variable based on 'mean_uf_cod'
df_clean$region <- sapply(df_clean$mean_uf_cod, state_to_region)

# Define a function to run the regressions and get summary
run_regression_summary <- function(data, response_var, independent_vars) {
  formula <- paste(response_var, "~", paste(independent_vars, collapse = " + "), "| ibge_cod + year")
  reg <- fixest::feols(as.formula(formula), data = data, panel.id = ~ibge_cod + year, cluster = ~ibge_cod)
  summary(reg)
}

# List of response variables and independent variables for each regression
response_vars <- c("slave", "operations", "slave", "operations")
independent_vars_list <- list(
  c("ss_imports", "log(mean_harvested_area_interp + 1)", "log(gdp_per_capita + 1)", "log(mean_number_cattle_interp + 1)", "log(mean_revenue_interp + 1)"),
  c("ss_imports", "log(mean_harvested_area_interp + 1)", "log(gdp_per_capita + 1)", "log(mean_number_cattle_interp + 1)", "log(mean_revenue_interp + 1)"),
  c("ss_taurus", "log(mean_harvested_area_interp + 1)", "log(gdp_per_capita + 1)", "log(mean_number_cattle_interp + 1)", "log(mean_revenue_interp + 1)"),
  c("ss_taurus", "log(mean_harvested_area_interp + 1)", "log(gdp_per_capita + 1)", "log(mean_number_cattle_interp + 1)", "log(mean_revenue_interp + 1)")
)

# Create an empty list to store the regression summaries
reg_summaries <- list()

# Iterate over the combinations of response and independent variables
for (i in seq_along(response_vars)) {
  response_var <- response_vars[i]
  independent_vars <- independent_vars_list[[i]]
  
  # Group the data by 'region'
  df_grouped <- df_clean %>%
    group_by(region)
  
  # Run the regression for each group and get the summary
  reg_summary <- df_grouped %>%
    do(regression_summary = run_regression_summary(., response_var, independent_vars))
  
  # Extract the regression summaries from the list and add them to the 'reg_summaries' list
  reg_summaries[[paste("Regression", i)]] <- reg_summary$regression_summary
}

# Print the regression summaries for each region
for (i in seq_along(reg_summaries)) {
  print(paste("Summary for Regression", i, ":", sep = " "))
  print(reg_summaries[[i]])
}

#####################################################################################################################################
############################### Region Plots: Slaves / Shift-share 1 (imports) ######################################################
#####################################################################################################################################

# Assuming 'df_clean' is the name of your dataset, and 'region' is the column with Brazilian region information

# Define a function to run the regression and get the coefficients
run_regression_coefficients <- function(data, response_var, independent_var, control_vars) {
  formula <- paste(response_var, "~", independent_var, "+", paste(control_vars, collapse = " + "), "| ibge_cod + year")
  reg <- by(data, data$region, function(df) fixest::feols(as.formula(formula), data = df, panel.id = ~ibge_cod + year, cluster = ~ibge_cod))
  coef_point_estimate <- sapply(reg, function(model) {
    coef_summary <- tidy(model)
    coef_summary$estimate[coef_summary$term == independent_var]
  })
  coef_se <- sapply(reg, function(model) {
    coef_summary <- tidy(model)
    coef_summary$std.error[coef_summary$term == independent_var]
  })
  coef_ci_low <- coef_point_estimate - 1.96 * coef_se
  coef_ci_high <- coef_point_estimate + 1.96 * coef_se
  return(data.frame(region = names(reg), coef_point_estimate, coef_ci_low, coef_ci_high))
}

# Response variable and independent variable (bartik_imports)
response_var <- "slave"
independent_var <- "ss_imports"

# Additional control variables
control_vars <- c("log(mean_harvested_area_interp+1)", "log(gdp_per_capita+1)", "log(mean_number_cattle_interp+1)", "log(mean_revenue_interp+1)")

# Create a data frame with the coefficients for each region
coefficients_df <- run_regression_coefficients(df_clean, response_var, independent_var, control_vars)

# Assuming 'coefficients_df' is correctly created from the 'run_regression_coefficients' function

# Print the summary of the regression for each region
for (i in 1:nrow(coefficients_df)) {
  region <- coefficients_df$region[i]
  cat("Summary for Region:", region, "\n")
  cat("Coefficient Intention to treat:", coefficients_df$coef_point_estimate[i], "\n") # Corrected line
  cat("95% CI Low:", coefficients_df$coef_ci_low[i], "\n")
  cat("95% CI High:", coefficients_df$coef_ci_high[i], "\n\n")
}

sl_bart1 <- ggplot(coefficients_df, aes(x = region, y = coef_point_estimate)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = coef_ci_low, ymax = coef_ci_high), width = 0.2, color = "blue") +
  geom_hline(yintercept = 0, lty = 2, lwd = 1, colour = "grey50") +
  labs(x = "Brazilian Regions", y = "Intention to treat", title = "Coefficient Plot for Shift-share 1 (imports) on Slavery") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the slaves_bartik_imports plot
ggsave(filename = "Slavesbartik_imports.pdf",  plot = sl_bart1, width = 12, height = 6, units = "in")

#####################################################################################################################################
############################### Region Plots: Audits / Shift-share 1 (imports) ######################################################
#####################################################################################################################################

# Assuming 'df_clean' is the name of your dataset, and 'region' is the column with Brazilian region information

# Define a function to run the regression and get the coefficients
run_regression_coefficients <- function(data, response_var, independent_var, control_vars) {
  formula <- paste(response_var, "~", independent_var, "+", paste(control_vars, collapse = " + "), "| ibge_cod + year")
  reg <- by(data, data$region, function(df) fixest::feols(as.formula(formula), data = df, panel.id = ~ibge_cod + year, cluster = ~ibge_cod))
  coef_point_estimate <- sapply(reg, function(model) {
    coef_summary <- tidy(model)
    coef_summary$estimate[coef_summary$term == independent_var]
  })
  coef_se <- sapply(reg, function(model) {
    coef_summary <- tidy(model)
    coef_summary$std.error[coef_summary$term == independent_var]
  })
  coef_ci_low <- coef_point_estimate - 1.96 * coef_se
  coef_ci_high <- coef_point_estimate + 1.96 * coef_se
  return(data.frame(region = names(reg), coef_point_estimate, coef_ci_low, coef_ci_high))
}

# Response variable and independent variable (operations and bartik_imports)
response_var <- "operations"
independent_var <- "ss_imports"

# Additional control variables
control_vars <- c("log(mean_harvested_area_interp+1)", "log(gdp_per_capita+1)", "log(mean_number_cattle_interp+1)", "log(mean_revenue_interp+1)")

# Create a data frame with the coefficients for each region
coefficients_df <- run_regression_coefficients(df_clean, response_var, independent_var, control_vars)

# Print the summary of the regression for each region
for (i in 1:nrow(coefficients_df)) {
  region <- coefficients_df$region[i]
  cat("Summary for Region:", region, "\n")
  cat("Coefficient Estimate:", coefficients_df$coef_point_estimate[i], "\n")  # Corrected this line
  cat("95% CI Low:", coefficients_df$coef_ci_low[i], "\n")
  cat("95% CI High:", coefficients_df$coef_ci_high[i], "\n\n")
}

audit_bart1 <- ggplot(coefficients_df, aes(x = region, y = coef_point_estimate)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = coef_ci_low, ymax = coef_ci_high), width = 0.2, color = "blue") +
  geom_hline(yintercept = 0, lty = 2, lwd = 1, colour = "grey50") +
  labs(x = "Brazilian Regions", y = "Intention to treat", title = "Coefficient Plot for Shift-share 1 (imports) on Audits") + # Corrected 'y' label and syntax
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the slaves_bartik_imports plot
ggsave(filename = "Auditbartik_imports.pdf",  plot = audit_bart1, width = 12, height = 6, units = "in")

#####################################################################################################################################
############################### Region Plots: Slaves / Shift-share 2 (Taurus) #######################################################
#####################################################################################################################################

# Define a function to run the regression and get the coefficients
run_regression_coefficients <- function(data, response_var, independent_var, control_vars) {
  formula <- paste(response_var, "~", independent_var, "+", paste(control_vars, collapse = " + "), "| ibge_cod + year")
  reg <- by(data, data$region, function(df) fixest::feols(as.formula(formula), data = df, panel.id = ~ibge_cod + year, cluster = ~ibge_cod))
  coef_point_estimate <- sapply(reg, function(model) {
    coef_summary <- tidy(model)
    coef_summary$estimate[coef_summary$term == independent_var]
  })
  coef_se <- sapply(reg, function(model) {
    coef_summary <- tidy(model)
    coef_summary$std.error[coef_summary$term == independent_var]
  })
  coef_ci_low <- coef_point_estimate - 1.96 * coef_se
  coef_ci_high <- coef_point_estimate + 1.96 * coef_se
  return(data.frame(region = names(reg), coef_point_estimate, coef_ci_low, coef_ci_high))
}

# Response variable and independent variable (operations and bartik_taurus)
response_var <- "slave"
independent_var <- "ss_taurus"

# Additional control variables
control_vars <- c("log(mean_harvested_area_interp+1)", "log(gdp_per_capita+1)", "log(mean_number_cattle_interp+1)", "log(mean_revenue_interp+1)")

# Create a data frame with the coefficients for each region
coefficients_df <- run_regression_coefficients(df_clean, response_var, independent_var, control_vars)

# Print the summary of the regression for each region
for (i in 1:nrow(coefficients_df)) {
  region <- coefficients_df$region[i]
  cat("Summary for Region:", region, "\n")
  cat("Intention to treat:", coefficients_df$coef_point_estimate[i], "\n") # Fixed this line
  cat("95% CI Low:", coefficients_df$coef_ci_low[i], "\n")
  cat("95% CI High:", coefficients_df$coef_ci_high[i], "\n\n")
}

sl_bart2 <- ggplot(coefficients_df, aes(x = region, y = coef_point_estimate)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = coef_ci_low, ymax = coef_ci_high), width = 0.2, color = "blue") +
  geom_hline(yintercept = 0, lty = 2, lwd = 1, colour = "grey50") +
  labs(x = "Brazilian Regions", y = "Intention to treat", title = "Coefficient Plot for Shift-share 2 (Taurus) on Slavery") + # Corrected this line
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the slaves_bartik_imports plot
ggsave(filename = "Slavesbartik_taurus.pdf",  plot = sl_bart2, width = 12, height = 6, units = "in")

#####################################################################################################################################
############################### Region Plots: Audits / Shift-share 2 (Taurus) #########################################################
#####################################################################################################################################

# Assuming 'df_clean' is the name of your dataset, and 'region' is the column with Brazilian region information

# Define a function to run the regression and get the coefficients
run_regression_coefficients <- function(data, response_var, independent_var, control_vars) {
  formula <- paste(response_var, "~", independent_var, "+", paste(control_vars, collapse = " + "), "| ibge_cod + year")
  reg <- by(data, data$region, function(df) fixest::feols(as.formula(formula), data = df, panel.id = ~ibge_cod + year, cluster = ~ibge_cod))
  coef_point_estimate <- sapply(reg, function(model) {
    coef_summary <- tidy(model)
    coef_summary$estimate[coef_summary$term == independent_var]
  })
  coef_se <- sapply(reg, function(model) {
    coef_summary <- tidy(model)
    coef_summary$std.error[coef_summary$term == independent_var]
  })
  coef_ci_low <- coef_point_estimate - 1.96 * coef_se
  coef_ci_high <- coef_point_estimate + 1.96 * coef_se
  return(data.frame(region = names(reg), coef_point_estimate, coef_ci_low, coef_ci_high))
}

# Response variable and independent variable (operations and bartik_taurus)
response_var <- "operations"
independent_var <- "ss_taurus"

# Additional control variables
control_vars <- c("log(mean_harvested_area_interp+1)", "log(gdp_per_capita+1)", "log(mean_number_cattle_interp+1)", "log(mean_revenue_interp+1)")

# Create a data frame with the coefficients for each region
coefficients_df <- run_regression_coefficients(df_clean, response_var, independent_var, control_vars)

# Print the summary of the regression for each region
for (i in 1:nrow(coefficients_df)) {
  region <- coefficients_df$region[i]
  cat("Summary for Region:", region, "\n")
  cat("Intention to treat:", coefficients_df$coef_point_estimate[i], "\n") # Corrected line
  cat("95% CI Low:", coefficients_df$coef_ci_low[i], "\n")
  cat("95% CI High:", coefficients_df$coef_ci_high[i], "\n\n")
}


audits_bartik_taurus <- ggplot(coefficients_df, aes(x = region, y = coef_point_estimate)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = coef_ci_low, ymax = coef_ci_high), width = 0.2, color = "blue") +
  geom_hline(yintercept = 0, lty = 2, lwd = 1, colour = "grey50") +
  labs(x = "Brazilian Regions", y = "Intention to treat", title = "Coefficient Plot for Shift-share 2 (Taurus) on Audits") + # Corrected line
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the audits_bartik_taurus plot
ggsave(filename = "Auditbartik_taurus.pdf",  plot = audits_bartik_taurus, width = 12, height = 6, units = "in")

#####################################################################################################################################
############################### Region Bar Plots ####################################################################################
#####################################################################################################################################

#### Audits ####

# Aggregate data for each region and calculate the sum of operations and slaves
df_agg_operations_slaves <- df_clean %>%
  group_by(region) %>%
  summarise(
    sum_operations = sum(operations),
    sum_slaves = sum(slave)
  )

# Create the bar plot for the sum of operations
plot_operations <- ggplot(df_agg_operations_slaves, aes(x = region)) +
  geom_col(aes(y = sum_operations, fill = "Audits"), width = 0.3, position = "dodge") +
  labs(title = "Audits per Region",
       x = "Region",
       y = "Frequency",
       fill = "Variable") +
  scale_fill_manual(values = '#DF0030') +
  theme_minimal() +
  geom_text(aes(y = sum_operations, label = round(sum_operations, 3)), vjust = -0.5, size = 3)

# Save the plots
ggsave("sum_of_operations_plot.pdf", plot = plot_operations, width = 10, height = 6, dpi = 300)

#### Slaves ####

# Create the bar plot for the sum of slaves
plot_slaves <- ggplot(df_agg_operations_slaves, aes(x = region)) +
  geom_col(aes(y = sum_slaves, fill = "Slaves"), width = 0.3, position = "dodge") +
  labs(title = "Slaves per Region",
       x = "Region",
       y = "Frequency",
       fill = "Variable") +
  scale_fill_manual(values = '#25D366') +
  theme_minimal() +
  geom_text(aes(y = sum_slaves, label = round(sum_slaves, 3)), vjust = -0.5, size = 3)

# Save the plots
ggsave("sum_of_slaves_plot.pdf", plot = plot_slaves, width = 10, height = 6, dpi = 300)


# Aggregate data for each region and calculate the sum of exposure
df_agg_exposure <- df_clean %>%
  group_by(region) %>%
  summarise(
    sum_exposure = sum(exposure)
  )

# Create the bar plot for the sum of exposure using Republican blue color
plot_exposure <- ggplot(df_agg_exposure, aes(x = region)) +
  geom_col(aes(y = sum_exposure, fill = "Exposure"), width = 0.3, position = "dodge") +
  labs(title = "Exposure per Region",
       x = "Region",
       y = "Frequency",
       fill = "Variable") +
  scale_fill_manual(values = '#007AFF') +  # Republican blue color
  theme_minimal() +
  geom_text(aes(y = sum_exposure, label = round(sum_exposure, 3)), vjust = -0.5, size = 3)

# Save the plot
ggsave("sum_of_exposure_plot.pdf", plot = plot_exposure, width = 10, height = 6, dpi = 300)



# Aggregate data for each region and calculate the sum of revenues
df_statecapacity <- df_clean %>%
  group_by(region) %>%
  summarise(
    sum_statecapacity = sum(mean_revenue_interp/mean_population_interp)
  )

# Create the bar plot for the sum of exposure using Republican blue color
plot_exposure <- ggplot(df_statecapacity, aes(x = region)) +
  geom_col(aes(y = sum_statecapacity, fill = " Average Revenue per Capita"), width = 0.3, position = "dodge") +
  labs(title = "State Capacity per Region",
       x = "Region",
       y = "Frequency",
       fill = "Variable") +
  scale_fill_manual(values = '#007AFF') +  # Republican blue color
  theme_minimal() +
  geom_text(aes(y = sum_statecapacity, label = round(sum_statecapacity, 3)), vjust = -0.5, size = 3)

# Save the plot
ggsave("state_capacity_plot.pdf", plot = plot_exposure, width = 10, height = 6, dpi = 300)

#####################################################################################################################################
############################# Budget Execution and Institutional Capacity ###########################################################
#####################################################################################################################################

# ====================================================================
# Reproducible Pipeline to Produce:
#   1) execution_ratios_trend.png  (2019–2025 budget execution ratios)
#   2) national_trend_appointments.png (2000–2019 total AFT appointments)
# --------------------------------------------------------------------
# Inputs expected (either in working directory OR at the Windows paths):
#   - Dados_Orcamentarios_2000_2025.txt   (budget for forced labor actions)
#   - AFT_99_a_2025_por_UF.csv            (AFT appointments by UF and year)
# Outputs (created in ./Resultados by default):
#   - execution_ratios_trend.png
#   - national_trend_appointments.png
# ====================================================================

## 0) Install/load packages ----------------------------------------------------
need <- c("data.table","dplyr","tidyr","ggplot2","scales","stringi","readr","grid")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(need, require, character.only = TRUE))

## 1) Paths --------------------------------------------------------------------
try_paths <- function(fn, win_path) {
  if (file.exists(win_path)) return(win_path)
  if (file.exists(fn))       return(fn)
  stop(paste0("File not found: '", fn, "'. Put it in the working dir or at: ", win_path))
}

budget_file <- try_paths(
  "Dados_Orcamentarios_2000_2025.txt",
  "C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/New_Data/LAI/Dados/Dados_Orcamentarios_2000_2025.txt"
)

aft_file <- try_paths(
  "AFT_99_a_2025_por_UF.csv",
  "C:/Users/gabic/Dropbox/Slaves_Shif-Share_IV/WD_RR/New_Data/LAI/Dados/AFT_99_a_2025_por_UF.csv"
)

out_dir <- "Resultados"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## 2) Helpers ------------------------------------------------------------------
# Convert Brazilian numeric strings like "1.234.567,89" -> 1234567.89
convert_br_num <- function(x) {
  x <- gsub("\\.", "", x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

# Normalize names: lowercase, remove accents, replace non-alphanumerics with "_"
norm_names <- function(x) {
  x2 <- stringi::stri_trans_general(x, "Latin-ASCII")
  x2 <- tolower(x2)
  x2 <- gsub("[^a-z0-9]+", "_", x2)
  x2
}

# Flexible column finder: tries each pattern until it finds a match
find_col <- function(patterns, nm_vec) {
  for (p in patterns) {
    m <- grep(p, nm_vec, perl = TRUE)
    if (length(m)) return(nm_vec[m[1]])
  }
  return(NA_character_)
}

## 3) Figure 1: execution_ratios_trend.png -------------------------------------
# Load budget data robustly (comma or tab separators)
orc <- tryCatch(fread(budget_file, sep = ",", header = TRUE),
                error = function(e) fread(budget_file, sep = "\t", header = TRUE))

# Standardize column names
setnames(orc, old = names(orc), new = norm_names(names(orc)))
nm <- names(orc)

# Identify columns (accepts PT/EN variants)
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

# Keep/rename minimal set
orc <- orc |>
  dplyr::select(
    Year = all_of(col_year),
    Action = all_of(col_act),
    Updated_Allocation = all_of(col_alloc),
    Committed = all_of(col_emp),
    Liquidated = all_of(col_liq),
    Paid = all_of(col_pago)
  )

# Coerce types
orc$Year <- as.integer(orc$Year)
for (v in c("Updated_Allocation","Committed","Liquidated","Paid")) {
  if (is.character(orc[[v]])) orc[[v]] <- convert_br_num(orc[[v]])
}

# Keep the forced-labor inspection action (regex tolerant)
act_norm <- norm_names(orc$Action)
orc$Action <- ifelse(grepl("erradicacao.*escr|combate_ao_trabalho_escravo|slave", act_norm, ignore.case = TRUE),
                     "Inspection to Eradicate Slave Labor", orc$Action)
orc_sl <- dplyr::filter(orc, Action == "Inspection to Eradicate Slave Labor")

# Build 2019–2025 series and compute ratios
orc_plot <- orc_sl |>
  dplyr::select(Year, Updated_Allocation, Committed, Liquidated, Paid) |>
  tidyr::complete(Year = 2019:2025) |>
  dplyr::mutate(
    dplyr::across(c(Updated_Allocation,Committed,Liquidated,Paid), ~tidyr::replace_na(.x, 0)),
    committed_share  = dplyr::if_else(Updated_Allocation > 0, Committed  / Updated_Allocation, NA_real_),
    liquidated_share = dplyr::if_else(Updated_Allocation > 0, Liquidated / Updated_Allocation, NA_real_),
    paid_share       = dplyr::if_else(Updated_Allocation > 0, Paid       / Updated_Allocation, NA_real_)
  )

# Long format with human-readable Type labels (no hyphens)
orc_long <- orc_plot |>
  dplyr::select(Year, committed_share, liquidated_share, paid_share) |>
  tidyr::pivot_longer(cols = -Year, names_to = "Type", values_to = "Share") |>
  dplyr::mutate(
    Type = factor(Type,
                  levels = c("committed_share","liquidated_share","paid_share"),
                  labels = c("Committed","Liquidated","Paid"))
  )

# Plot 1: framed panel, no title, same visual language as plot 2
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
    panel.border      = element_rect(colour = "grey60", fill = NA, linewidth = 0.6) # frame
  )

ggsave(file.path(out_dir, "execution_ratios_trend.png"),
       p_exec, width = 8, height = 5, bg = "white", dpi = 300)

## 4) Figure 2: national_trend_appointments.png --------------------------------
aft <- readr::read_csv(aft_file, show_col_types = FALSE)
names(aft) <- norm_names(names(aft))

# Expect columns like "ano" (year) and "total_geral" (national total)
if (!all(c("ano","total_geral") %in% names(aft))) {
  stop("Could not find columns for year ('ano') and national total ('total_geral'). Found: ",
       paste(names(aft), collapse = ", "))
}

df <- aft |>
  dplyr::filter(!is.na(ano)) |>
  dplyr::mutate(ano = as.integer(ano)) |>
  dplyr::filter(ano >= 2000, ano <= 2019)

# Plot 2: plain white panel, no title, harmonized axes/lines with plot 1
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

message("Done. PNGs saved to: ", normalizePath(out_dir))
