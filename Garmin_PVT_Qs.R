##----------------------------------------------------------------------------------------------- ##
#                                       COMPARE DATA UNITS
##----------------------------------------------------------------------------------------------- ##
#
# Description:  This script assesses the associations between Garmin, PVT, and questionnaire data.
#
# Date:         Aug 2025
# R.version:    4.4.0 (2024-04-24)
# Rstudio:      2025.05.0+496
#
## ---------------------------------------------------------------------------------------------- ##

# clear environment
rm(list=ls())

# ------------------------------------------------------------------------------------------------ #
#                                      Settings & Dependencies
# ------------------------------------------------------------------------------------------------ #

# numbers of external volumes
heronderzoek <- "heronderzoek"


# Import libraries
# ----------------------------------- #

# read files
library(readxl)
library(writexl)

# data manipulation
library(tidyverse)
library(dplyr)

# data imputation


# statistics
library(ltm)
library(psych)
library(rstatix)
library(lme4)
library(lmerTest)

# visualizations
library(ggplot2)
library(ggpubr)
library(GGally)
library(corrplot)


# Functions
# ----------------------------------- #



# Define working directory
# ----------------------------------- #

# external volume, change value if necessary 
#setwd(paste("~/Networkshares/", heronderzoek, "/Groep Geuze/25U-0078_PRESS/", sep = ""))
setwd("/Volumes/heronderzoek/Groep Geuze/25U-0078_PRESS/")


# ~/Networkshares/ of ~/Volumes/





# ------------------------------------------------------------------------------------------------ #
#                                   Data Collection and preparation
# ------------------------------------------------------------------------------------------------ #

# Read garmin
# ----------------------------------- #
df_garmin <- read_csv("E_ResearchData/2_ResearchData/1. Verwerkte data/Garmin/Garmin_data_clean2.csv", show_col_types = FALSE)

df_garmin_avg <- df_garmin %>%
  group_by(Castor.ID, measurement) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  filter(measurement != "x - Transitiedag")

# Read PVT
# ----------------------------------- #
df_pvt <- read_csv("E_ResearchData/2_ResearchData/1. Verwerkte data/PVT/PVT_readouts.csv", show_col_types = FALSE)

#change one entered non-existing PVT-ID at T2 (206), most likely is a typo for 906:
df_pvt$PVTid[which(df_pvt$PVTid==206)] <- 906

df_codelist <- data.frame(read_excel("C_PersonalData/2_CodeList/Sleutelbestand_PRESS.xlsx"))  %>%
  select("Castor.ID", "PVT.ID")  %>%
  mutate(PVT.ID = as.numeric(PVT.ID))

df_pvt <- full_join(df_pvt, df_codelist, by = c("PVTid" = "PVT.ID"))

df_pvt <- df_pvt %>%
  mutate(measurement = case_when(timepoint == "M1" ~ "T0",
                                 timepoint == "M2" ~ "T1",
                                 timepoint == "M3" ~ "T2",
                                 timepoint == "M4" ~ "T3",
                                 TRUE ~ NA_character_))


# Read questionnaires
# ----------------------------------- #
df_in_slaap_vallen <- read_csv("E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_In_slaap_vallen_export_20250514_clean.csv", show_col_types = FALSE)
df_spanning        <- read_csv("E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_Spanning_export_20250514_clean.csv", show_col_types = FALSE)
df_vermoeidheid    <- read_csv("E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_Vermoeidheid_export_20250514_clean.csv", show_col_types = FALSE)
df_vas             <- read_csv("E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_Alertheid_prikkelbaarheid_en_schrikachtigheid_export_20250514_clean.csv", show_col_types = FALSE)

df_in_slaap_vallen <- df_in_slaap_vallen %>%
  mutate(measurement = case_when(Survey.Package.Name == "Meting 1 vragenlijsten" ~ "T0",
                                 Survey.Package.Name == "Meting 2 vragenlijsten" ~ "T1",
                                 Survey.Package.Name == "Meting 3 vragenlijsten" ~ "T2",
                                 Survey.Package.Name == "Meting 4 vragenlijsten" ~ "T3",
                                 TRUE ~ NA_character_),
         Castor.ID = Castor.Participant.ID)

df_spanning <- df_spanning %>%
  mutate(measurement = case_when(Survey.Package.Name == "Meting 1 vragenlijsten" ~ "T0",
                                 Survey.Package.Name == "Meting 2 vragenlijsten" ~ "T1",
                                 Survey.Package.Name == "Meting 3 vragenlijsten" ~ "T2",
                                 Survey.Package.Name == "Meting 4 vragenlijsten" ~ "T3",
                                 TRUE ~ NA_character_),
         Castor.ID = Castor.Participant.ID)

df_vermoeidheid <- df_vermoeidheid %>%
  mutate(measurement = case_when(Survey.Package.Name == "Meting 1 vragenlijsten" ~ "T0",
                                 Survey.Package.Name == "Meting 2 vragenlijsten" ~ "T1",
                                 Survey.Package.Name == "Meting 3 vragenlijsten" ~ "T2",
                                 Survey.Package.Name == "Meting 4 vragenlijsten" ~ "T3",
                                 TRUE ~ NA_character_),
         Castor.ID = Castor.Participant.ID)

df_vas <- df_vas %>%
  mutate(measurement = case_when(Survey.Package.Name == "Meting 1 vragenlijsten" ~ "T0",
                                 Survey.Package.Name == "Meting 2 vragenlijsten" ~ "T1",
                                 Survey.Package.Name == "Meting 3 vragenlijsten" ~ "T2",
                                 Survey.Package.Name == "Meting 4 vragenlijsten" ~ "T3",
                                 TRUE ~ NA_character_),
         Castor.ID = Castor.Participant.ID)


# Read sleep diary
# ----------------------------------- #
df_slaapdagboek1 <- read_csv("E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_Slaapdagboek_1_nacht_export_20250514_clean.csv",show_col_types = FALSE)

df_slaapdagboek1 <- df_slaapdagboek1 %>%
  mutate(measurement = case_when(Survey.Package.Name == "Meting 1 vragenlijsten" ~ "T0",
                                 Survey.Package.Name == "Meting 2 vragenlijsten" ~ "T1",
                                 Survey.Package.Name == "Meting 3 vragenlijsten" ~ "T2",
                                 Survey.Package.Name == "Meting 4 vragenlijsten" ~ "T3",
                                 TRUE ~ NA_character_),
         Castor.ID = Castor.Participant.ID)

# Calculate total sleep time (TST)
df_slaapdagboek1$TST <- NA
df_slaapdagboek1 <- df_slaapdagboek1 %>%
  group_by( c(Castor.Participant.ID, Survey.Completed.On) ) %>%
  summarise(across(.cols = where(is.numeric), .fns = ~ if (all(is.na(.))) NA else mean(., na.rm = TRUE)), .groups = "drop") %>%
  
  


# ------------------------------------------------------------------------------------------------ #
#                                          Merge data frames
# ------------------------------------------------------------------------------------------------ #


convert_measurement <- function(df) {
  df %>%
    mutate(
      measurement = case_when(
        str_detect(Survey.Package.Name, "Meting 1") ~ "T0",
        str_detect(Survey.Package.Name, "Meting 2") ~ "T1",
        str_detect(Survey.Package.Name, "Meting 3") ~ "T2",
        str_detect(Survey.Package.Name, "Meting 4") ~ "T3",
        TRUE ~ NA_character_
      ),
      Castor.ID = Castor.Participant.ID
    ) %>%
    select(Castor.ID, measurement, sum)
}

df_pvt <- df_pvt %>%
  select(Castor.ID, measurement, rt_mean) %>%
  filter(measurement != "T2") #skip T2 due to missing data for almost all participants

df_in_slaap_vallen_sum <- convert_measurement(df_in_slaap_vallen) %>%
  rename(sum_slaap = sum)

df_spanning_sum <- convert_measurement(df_spanning) %>%
  rename(sum_spanning = sum)

df_vermoeidheid_sum <- convert_measurement(df_vermoeidheid) %>%
  rename(sum_vermoeidheid = sum)

df_vas <- df_vas %>%
  select("Castor.ID", "measurement", "alertness", "irritability", "skittish")
  

# Merge all with df_garmin_avg
df_merged <- df_garmin_avg %>%
  full_join(df_pvt, by = c("Castor.ID", "measurement")) %>%
  full_join(df_in_slaap_vallen_sum, by = c("Castor.ID", "measurement")) %>%
  full_join(df_spanning_sum, by = c("Castor.ID", "measurement")) %>%
  full_join(df_vermoeidheid_sum, by = c("Castor.ID", "measurement")) %>%
  full_join(df_vas, by = c("Castor.ID", "measurement"))  




# ------------------------------------------------------------------------------------------------ #
#                               Add T2/T3 and delta's
# ------------------------------------------------------------------------------------------------ #

# Combine T1 and T2 into one measurement
# -------------------------------------- #

df_t1_t2 <- df_merged %>%
  filter(measurement %in% c("T1", "T2")) %>%
  group_by(Castor.ID) %>%
  summarise(across(.cols = where(is.numeric), .fns = ~ if (all(is.na(.))) NA else mean(., na.rm = TRUE)), .groups = "drop") %>%
  mutate(measurement = "T1/2") %>%
  select(colnames(df_merged))

# Add to original
df_merged <- bind_rows(df_merged, df_t1_t2) %>%
  arrange(Castor.ID, measurement) %>%
  mutate(
    # Zet character "NA" en lege strings om in NA voor character-kolommen
    across(where(is.character), ~ na_if(., "NA")),
    across(where(is.character), ~ na_if(., "")),
    
    # Zet NaN om in NA voor numerieke kolommen
    across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))
rm(df_t1_t2)



# Calculate delta between measurement rounds
# --------------------------------------------- #

# Step 1: T1 vs T2-T3
df_diff1 <- df_merged %>%
  filter(measurement %in% c("T0", "T1/2")) %>%
  pivot_wider(
    names_from = measurement,
    values_from = c(HRV.Last.Night.Avg...ms., Resting.HR..bpm., Avg.Sleep.min, Nap.Duration.min,
                    HRV.Last.Night.High..ms., Min..Daily.HR..bpm., Avg.Sleep.Awake.Dur.min, Sleep.Score,
                    rt_mean, sum_slaap, sum_spanning, sum_vermoeidheid, alertness, irritability, skittish),
    names_sep = "_"
  ) %>%
  mutate(across(
    ends_with("_T1/2"),
    ~ ifelse(is.na(.x) | is.na(get(sub("_T1/2$", "_T1", cur_column()))), NA_real_, .x - get(sub("_T1/2$", "_T0", cur_column()))),
    .names = "{sub('_T1/2$', '', .col)}"
  )) %>%
  select(Castor.ID, all_of(c(
    "HRV.Last.Night.Avg...ms.", "Resting.HR..bpm.", "Avg.Sleep.min", "Nap.Duration.min",
    "HRV.Last.Night.High..ms.", "Min..Daily.HR..bpm.", "Avg.Sleep.Awake.Dur.min", "Sleep.Score",
    "rt_mean", "sum_slaap", "sum_spanning", "sum_vermoeidheid", "alertness", "irritability", "skittish"
  ))) %>%
  mutate(measurement = "delta_T0_T1/2") %>%
  select(Castor.ID, measurement, everything())


# Step 2: T1-T2 vs T3
df_diff2 <- df_merged %>%
  filter(measurement %in% c("T1/2", "T3")) %>%
  pivot_wider(
    names_from = measurement,
    values_from = c(HRV.Last.Night.Avg...ms., Resting.HR..bpm., Avg.Sleep.min, Nap.Duration.min,
                    HRV.Last.Night.High..ms., Min..Daily.HR..bpm., Avg.Sleep.Awake.Dur.min, Sleep.Score,
                    rt_mean, sum_slaap, sum_spanning, sum_vermoeidheid, alertness, irritability, skittish),
    names_sep = "_"
  ) %>%
  mutate(across(
    ends_with("_T1/2"),
    ~ ifelse(is.na(.x) | is.na(get(sub("_T1/2$", "_T3", cur_column()))), NA_real_, get(sub("_T1/2$", "_T3", cur_column())) - .x),
    .names = "{sub('_T1/2$', '', .col)}"
  )) %>%
  select(Castor.ID, all_of(c(
    "HRV.Last.Night.Avg...ms.", "Resting.HR..bpm.", "Avg.Sleep.min", "Nap.Duration.min",
    "HRV.Last.Night.High..ms.", "Min..Daily.HR..bpm.", "Avg.Sleep.Awake.Dur.min", "Sleep.Score",
    "rt_mean", "sum_slaap", "sum_spanning", "sum_vermoeidheid", "alertness", "irritability", "skittish"
  ))) %>%
  mutate(measurement = "delta_T1/2_T3") %>%
  select(Castor.ID, measurement, everything())

# Stap 4: Combine
df_final <- bind_rows(df_merged, df_diff1, df_diff2) %>%
  arrange(Castor.ID, factor(measurement, levels = c("T0", "T1", "T1/2", "T2", "T3", "T0_T1/2", "T1/2_T3")))
rm(df_diff1)
rm(df_diff2)





# remove castor ID's 110019 omdat Garmin data mist
df_final <- df_final %>%
  filter(!Castor.ID %in% c(110019))
# save data frame with all variables
write.csv(df_final,"E_ResearchData/2_ResearchData/1. Verwerkte data/PRESS_all_variables_clean.csv", row.names = FALSE)




# ------------------------------------------------------------------------------------------------ #
#                                   Calculate correlations
# ------------------------------------------------------------------------------------------------ #


# Create correlation matrix
# ----------------------------------- #
measurements <- c("T0","T1", "T2", "T3")

for (measurement_point in measurements) {
  
  # Filter and select numeric variables
  numeric_vars <- df_final %>%
    filter(measurement == measurement_point) %>%
    mutate(Castor.ID = as.character(Castor.ID)) %>%
    select(where(is.numeric)) %>%
    select(where(~ !all(is.na(.x))))
  
  # Skip if not enough variables
  if (ncol(numeric_vars) < 2) next
  
  # Compute correlation matrix
  cor_matrix <- cor(numeric_vars, use = "pairwise.complete.obs")
  
  fileName <- paste0("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/",
                     "corrPlot_Garmin_PVT_Qs_", measurement_point, ".png")
  
  # Open SVG device
  png(filename = fileName, width = 800, height = 600, res = 150)
  corrplot.mixed(cor_matrix,
                 upper = "ellipse",
                 lower = "number",
                 tl.cex = 0.5,
                 number.cex = 0.5,
                 tl.pos = "lt")
  dev.off()
}





# ------------------------------------------------------------------------------------------------ #
#                                   Visualize trajectories
# ------------------------------------------------------------------------------------------------ #

# Vars to plot
vars_to_plot <- c("Avg.Sleep.min", "sum_vermoeidheid", "alertness", "irritability", "skittish", "rt_mean")

df_plot <- df_merged %>%
  filter(measurement %in% c("T0", "T1/2", "T3")) %>%
  mutate(measurement = factor(measurement, levels = c("T0", "T1/2", "T3")))

for (varname in vars_to_plot) {
  # Plot met facets per deelnemer
  p_facet <- ggplot(df_plot, aes(x = measurement, y = .data[[varname]], group = Castor.ID)) +
    geom_line(alpha = 0.7) +
    geom_point(alpha = 0.7) +
    facet_wrap(~ Castor.ID, scales = "free_y") +
    labs(
      title = paste("Verloop van", varname, "per deelnemer (facets)"),
      x = "Measurement",
      y = varname
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  fileName <- paste0("verloop_", varname, ".svg")
  ggsave(paste0("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/", fileName),
         plot = p_facet, width = 260, height = 200, units = "mm")
  
  # Plot met alle deelnemers in 1 figuur zonder facets
  p_all <- ggplot(df_plot, aes(x = measurement, y = .data[[varname]], group = Castor.ID, color = factor(Castor.ID))) +
    geom_line(alpha = 0.7) +
    geom_point(alpha = 0.7) +
    labs(
      title = paste("Verloop van", varname, "voor alle deelnemers"),
      x = "Measurement",
      y = varname,
      color = "Deelnemer"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  fileName <- paste0("verloop_alleDeelnemers", varname, ".svg")
  ggsave(paste0("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/", fileName),
         plot = p_all, width = 260, height = 200, units = "mm")
}


# ------------------------------------------------------------------------------------------------ #
#                                   Assessment of associations
# ------------------------------------------------------------------------------------------------ #


# Correlation
# ----------------------------------- #

df_corr <- df_final %>%
  filter(measurement == "delta_T0_T1/2") %>%
  select("Avg.Sleep.min", "sum_vermoeidheid", "alertness", "irritability", "skittish", "rt_mean")

# Compute correlation matrix
cor_matrix <- cor(df_corr, use = "pairwise.complete.obs")
corrplot.mixed(cor_matrix,
               upper = "ellipse",
               lower = "number",
               tl.cex = 0.5,
               number.cex = 0.5,
               tl.pos = "lt")



# Regression models 
# ------------------------------------------------------ #

# Algemene functie om modellen te draaien
run_models <- function(predictor_var, predictor_name, outcomes, df_final) {
  
  # Data voor delta_T1_T2/T3
  df_delta <- df_final %>%
    filter(measurement == "delta_T0_T1/2") %>%
    select(Castor.ID, !!predictor_name := all_of(predictor_var), all_of(outcomes))
  
  # Data voor T2/T3
  df_T2T3 <- df_final %>%
    filter(measurement == "T1/2") %>%
    select(Castor.ID, all_of(outcomes)) %>%
    rename_with(~ paste0(.x, "_T1T2"), -Castor.ID)
  
  # Combineer in één dataframe
  df_combined <- df_delta %>%
    left_join(df_T2T3, by = "Castor.ID")
  
  # Functie voor één outcome
  run_both_models <- function(outcome) {
    # delta
    f1 <- as.formula(paste0(outcome, " ~ ", predictor_name))
    model_delta <- lm(f1, data = df_combined)
    res1 <- tidy(model_delta) %>%
      filter(term == predictor_name) %>%
      mutate(model = "delta_vs_delta",
             outcome = outcome,
             r_squared = summary(model_delta)$r.squared)
    
    # T2/T3
    outcome_T2T3 <- paste0(outcome, "_T1T2")
    f2 <- as.formula(paste0(outcome_T2T3, " ~ ", predictor_name))
    model_T2T3 <- lm(f2, data = df_combined)
    res2 <- tidy(model_T2T3) %>%
      filter(term == predictor_name) %>%
      mutate(model = "delta_vs_T1T2",
             outcome = outcome,
             r_squared = summary(model_T2T3)$r.squared)
    
    bind_rows(res1, res2)
  }
  
  # Run voor alle outcomes
  map_df(outcomes, run_both_models) %>%
    select(model, outcome, estimate, std.error, statistic, p.value, r_squared) %>%
    arrange(model, p.value)
}

# Outcomes
outcomes <- c("rt_mean", "sum_vermoeidheid", "alertness", "irritability", "skittish")

# Draai modellen voor elke voorspeller
result_sleep_min <- run_models("Avg.Sleep.min", "delta_sleep_min", outcomes, df_final)
result_sleep_score <- run_models("Sleep.Score", "delta_sleep_score", outcomes, df_final)
result_hrv <- run_models("HRV.Last.Night.Avg...ms.", "delta_hrv", outcomes, df_final)

# Printen
print(result_sleep_min)
print(result_sleep_score)
print(result_hrv)


# Chatgpt
# Welke resultaten zijn interessant/significant?
#   🔹 Slaapduur (result_sleep_min)
# alertness in delta_vs_T2T3:
#   p = 0.030, estimate = 0.0950, R² = 0.335 → significant positief effect, dus meer slaap hangt samen met meer alertheid bij T2/T3.
# 
# Overige resultaten zijn niet significant (p > 0.05), hoewel het effect op irritability in delta_vs_delta een trend laat zien (p = 0.124).
# 
# 🔹 Slaapscore (result_sleep_score)
# irritability in delta_vs_delta:
#   p = 0.0540, estimate = 0.700, R² = 0.275 → bijna significant, mogelijk relevant effect.
# 
# Overige resultaten zijn niet significant, al zijn sommige trends zichtbaar (zoals alertness bij delta_vs_T2T3, p = 0.104).
# 
# 🔹 HRV (result_hrv)
# irritability in delta_vs_delta:
#   p = 0.0482, estimate = 0.700, R² = 0.310 → significant positief effect → toename in HRV hangt samen met minder irritatie.
# 
# Andere resultaten zijn niet significant (maar skittish en alertness laten wel zwakke trends zien).













# andere manier
df_selected <- df_final %>%
  filter(measurement %in% c("T1", "T2/T3", "T4"))
model <- lmer(irritability ~ Avg.Sleep.min * measurement + (1 | Castor.ID), data = df_selected)
summary(model)












