##----------------------------------------------------------------------------------------------- ##
#                                 QUESTIONNAIRE DATA
##----------------------------------------------------------------------------------------------- ##
#
# Description:  This script cleans and describes the questionnaire data.
#
# Date:         May 2025
# R.version:    4.4.0 (2024-04-24)
# Rstudio:      2025.05.0+496
#
## ---------------------------------------------------------------------------------------------- ##


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

# visualizations
library(ggplot2)
library(ggpubr)

# Functions
# ----------------------------------- #



# Define working directory
# ----------------------------------- #

# external volume, change value if necessary 
setwd(paste("~/Networkshares/", heronderzoek, "/Groep Geuze/25U-0078_PRESS/", sep = ""))

# ~/Networkshares/ of ~/Volumes/

# ------------------------------------------------------------------------------------------------ #
#                                          Data Collection
# ------------------------------------------------------------------------------------------------ #



# ----------------------------------- #
# In slaap vallen
# ----------------------------------- #
df_in_slaap_vallen <- data.frame(read_delim("E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_In_slaap_vallen_export_20250514.csv", delim = ";", show_col_types = FALSE))
df_in_slaap_vallen <- df_in_slaap_vallen %>%
  filter(Survey.Progress == 100) %>%
  dplyr::select("Castor.Participant.ID", "Survey.Package.Name", starts_with("PSAS_")) %>%
  mutate(sum = rowSums(dplyr::select(., starts_with("PSAS_")), na.rm = TRUE))

write.csv(df_in_slaap_vallen,"E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_In_slaap_vallen_export_20250514_clean.csv", row.names = FALSE)


# ----------------------------------- #
# Spanning
# ----------------------------------- #
df_spanning <- data.frame(read_delim("E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_Spanning_export_20250514.csv", delim = ";", show_col_types = FALSE))
df_spanning <- df_spanning %>%
  filter(Survey.Progress == 100) %>%
  dplyr::select("Castor.Participant.ID", "Survey.Package.Name", starts_with("STAI_"))

# reverse omkeeritems en sla op
stai_items <- paste0("STAI_0", 1:6)
reverse_items <- c("STAI_01", "STAI_04", "STAI_05") # positive items need to be reversed (kalm, ontspannen en tevreden)

names(df_spanning)[names(df_spanning) %in% stai_items] <- paste0(stai_items, "_raw") #keep old column names

for (item in stai_items) {
  raw_col <- paste0(item, "_raw")
  new_col <- paste0(item, "_reversed")
  
  if (item %in% reverse_items) {
    df_spanning[[new_col]] <- 5 - df_spanning[[raw_col]]
  } else {
    df_spanning[[new_col]] <- df_spanning[[raw_col]]
  }
}
write.csv(df_spanning,"E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_Spanning_export_20250514_reversed.csv", row.names = FALSE)

df_spanning <- df_spanning %>%
  mutate(sum = rowSums(dplyr::select(., ends_with("_reversed")), na.rm = TRUE))

write.csv(df_spanning,"E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_Spanning_export_20250514_clean.csv", row.names = FALSE)

# ----------------------------------- #
# Vermoeidheid
# ----------------------------------- #
df_vermoeidheid <- data.frame(read_delim("E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_Vermoeidheid_export_20250514.csv", delim = ";", show_col_types = FALSE))
df_vermoeidheid <- df_vermoeidheid %>%
  filter(Survey.Progress == 100) %>%
  dplyr::select("Castor.Participant.ID", "Survey.Package.Name", starts_with("PROMIS_")) %>%
  mutate(sum = rowSums(dplyr::select(., starts_with("PROMIS_")), na.rm = TRUE))

write.csv(df_vermoeidheid,"E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_Vermoeidheid_export_20250514_clean.csv", row.names = FALSE)


# ----------------------------------- #
# VAS
# ----------------------------------- #
df_vas <- data.frame(read_delim("E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_Alertheid_prikkelbaarheid_en_schrikachtigheid_export_20250514.csv", delim = ";", show_col_types = FALSE))
df_vas <- df_vas %>%
  filter(Survey.Progress == 100) %>%
  dplyr::select("Castor.Participant.ID", "Survey.Package.Name", "alertness", "irritability", "skittish")

write.csv(df_vas,"E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_Alertheid_prikkelbaarheid_en_schrikachtigheid_export_20250514_clean.csv", row.names = FALSE)


# ------------------------------------------------------------------------------------------------ #
#                                     Cronbach's alpha
# ------------------------------------------------------------------------------------------------ #

# In slaap vallen
# ----------------------------------- #
ca_in_slaap_vallen <- cronbach.alpha(df_in_slaap_vallen %>% dplyr::select(starts_with("PSAS_")))[["alpha"]]
ca_m1_in_slaap_vallen <- cronbach.alpha(df_in_slaap_vallen %>% 
                          dplyr::filter(Survey.Package.Name == "Meting 1 vragenlijsten") %>%
                          dplyr::select(starts_with("PSAS_")))[["alpha"]]
ca_m2_in_slaap_vallen <- cronbach.alpha(df_in_slaap_vallen %>% 
                          dplyr::filter(Survey.Package.Name == "Meting 2 vragenlijsten") %>%
                          dplyr::select(starts_with("PSAS_")))[["alpha"]]
ca_m3_in_slaap_vallen <- cronbach.alpha(df_in_slaap_vallen %>% 
                          dplyr::filter(Survey.Package.Name == "Meting 3 vragenlijsten") %>%
                          dplyr::select(starts_with("PSAS_")))[["alpha"]]
ca_m4_in_slaap_vallen <- cronbach.alpha(df_in_slaap_vallen %>%
                          dplyr::filter(Survey.Package.Name == "Meting 4 vragenlijsten") %>%
                          dplyr::select(starts_with("PSAS_")))[["alpha"]]



# Spanning
# ----------------------------------- #
df_spanning <- na.omit(df_spanning) # remove the row with one missing value
ca_spanning <- cronbach.alpha(df_spanning %>% dplyr::select(ends_with("_reversed")))[["alpha"]]
ca_m1_spanning <- cronbach.alpha(df_spanning %>% 
                          dplyr::filter(Survey.Package.Name == "Meting 1 vragenlijsten") %>%
                          dplyr::select(ends_with("_reversed")))[["alpha"]]
ca_m2_spanning <- cronbach.alpha(df_spanning %>% 
                          dplyr::filter(Survey.Package.Name == "Meting 2 vragenlijsten") %>%
                          dplyr::select(ends_with("_reversed")))[["alpha"]]
ca_m3_spanning <- cronbach.alpha(df_spanning %>% 
                          dplyr::filter(Survey.Package.Name == "Meting 3 vragenlijsten") %>%
                          dplyr::select(ends_with("_reversed")))[["alpha"]]
ca_m4_spanning <- cronbach.alpha(df_spanning %>%
                          dplyr::filter(Survey.Package.Name == "Meting 4 vragenlijsten") %>%
                          dplyr::select(ends_with("_reversed")))[["alpha"]]



# Vermoeidheid
# ----------------------------------- #
ca_vermoeidheid <- cronbach.alpha(df_vermoeidheid %>% dplyr::select(starts_with("PROMIS_")))[["alpha"]]
ca_m1_vermoeidheid <- cronbach.alpha(df_vermoeidheid %>% 
                                       dplyr::filter(Survey.Package.Name == "Meting 1 vragenlijsten") %>%
                                       dplyr::select(starts_with("PROMIS_")))[["alpha"]]
ca_m2_vermoeidheid <- cronbach.alpha(df_vermoeidheid %>% 
                                       dplyr::filter(Survey.Package.Name == "Meting 2 vragenlijsten") %>%
                                       dplyr::select(starts_with("PROMIS_")))[["alpha"]]
ca_m3_vermoeidheid <- cronbach.alpha(df_vermoeidheid %>% 
                                       dplyr::filter(Survey.Package.Name == "Meting 3 vragenlijsten") %>%
                                       dplyr::select(starts_with("PROMIS_")))[["alpha"]]
ca_m4_vermoeidheid <- cronbach.alpha(df_vermoeidheid %>%
                                       dplyr::filter(Survey.Package.Name == "Meting 4 vragenlijsten") %>%
                                       dplyr::select(starts_with("PROMIS_")))[["alpha"]]



# ------------------------------------------------------------------------------------------------ #
#                     Test-hertest reliability + repeated measures ANOVA
# ------------------------------------------------------------------------------------------------ #

# Stap 1: Kies vragenlijst
vragenlijst <- "vermoeidheid"
df <- df_vermoeidheid

df$Survey.Package.Name <- factor(df$Survey.Package.Name,
                                 levels = c("Meting 1 vragenlijsten", "Meting 2 vragenlijsten", "Meting 3 vragenlijsten", "Meting 4 vragenlijsten"))

# Stap 2: Filter originele data op alleen die deelnemers
df_volledige_metingen <- df %>%
  group_by(Castor.Participant.ID) %>%
  filter(n_distinct(Survey.Package.Name) == 4) %>%
  ungroup()

# Stap 3: Bereken reliability 
cor(df_volledige_metingen[df_volledige_metingen$Survey.Package.Name == "Meting 1 vragenlijsten", "sum"], 
    df_volledige_metingen[df_volledige_metingen$Survey.Package.Name == "Meting 2 vragenlijsten", "sum"], 
    use = "complete.obs", method = "pearson")

cor(df_volledige_metingen[df_volledige_metingen$Survey.Package.Name == "Meting 1 vragenlijsten", "sum"], 
    df_volledige_metingen[df_volledige_metingen$Survey.Package.Name == "Meting 3 vragenlijsten", "sum"], 
    use = "complete.obs", method = "pearson")

cor(df_volledige_metingen[df_volledige_metingen$Survey.Package.Name == "Meting 1 vragenlijsten", "sum"], 
    df_volledige_metingen[df_volledige_metingen$Survey.Package.Name == "Meting 4 vragenlijsten", "sum"], 
    use = "complete.obs", method = "pearson")

cor(df_volledige_metingen[df_volledige_metingen$Survey.Package.Name == "Meting 2 vragenlijsten", "sum"], 
    df_volledige_metingen[df_volledige_metingen$Survey.Package.Name == "Meting 3 vragenlijsten", "sum"], 
    use = "complete.obs", method = "pearson")

cor(df_volledige_metingen[df_volledige_metingen$Survey.Package.Name == "Meting 2 vragenlijsten", "sum"], 
    df_volledige_metingen[df_volledige_metingen$Survey.Package.Name == "Meting 4 vragenlijsten", "sum"], 
    use = "complete.obs", method = "pearson")

cor(df_volledige_metingen[df_volledige_metingen$Survey.Package.Name == "Meting 3 vragenlijsten", "sum"], 
    df_volledige_metingen[df_volledige_metingen$Survey.Package.Name == "Meting 4 vragenlijsten", "sum"], 
    use = "complete.obs", method = "pearson")






# Combine dataframes
# ----------------------------------- #

extract_measurement <- function(x) {
  str_extract(x, "Meting \\d") %>% str_replace("Meting ", "T")
}

df_in_slaap_vallen <- df_in_slaap_vallen %>%
  mutate(measurement = extract_measurement(Survey.Package.Name)) %>%
  select(Castor.ID = Castor.Participant.ID, measurement, sum_in_slaap_vallen = sum)

df_spanning <- df_spanning %>%
  mutate(measurement = extract_measurement(Survey.Package.Name)) %>%
  select(Castor.ID = Castor.Participant.ID, measurement, sum_spanning = sum)

df_vermoeidheid <- df_vermoeidheid %>%
  mutate(measurement = extract_measurement(Survey.Package.Name)) %>%
  select(Castor.ID = Castor.Participant.ID, measurement, sum_vermoeidheid = sum)

df_vas <- df_vas %>%
  mutate(measurement = extract_measurement(Survey.Package.Name)) %>%
  select(Castor.ID = Castor.Participant.ID, measurement, alertness, irritability, skittish)

# Merge
df_sumscores <- df_in_slaap_vallen %>%
  full_join(df_spanning, by = c("Castor.ID", "measurement")) %>%
  full_join(df_vermoeidheid, by = c("Castor.ID", "measurement")) %>%
  full_join(df_vas, by = c("Castor.ID", "measurement"))



# Repeated measures ANOVA/Friedman
# --------------------------------------- #
cols_of_interest <- c("sum_in_slaap_vallen", "sum_spanning", "sum_vermoeidheid",
                      "alertness", "irritability", "skittish")

for (col in cols_of_interest) {
  print(col)
  
  # Filter complete cases per variable
  df_var <- df_sumscores %>%
    select(Castor.ID, measurement, all_of(col)) %>%
    filter(!is.na(.data[[col]]))
  
  # Set factor levels for measurement
  df_var$measurement <- factor(df_var$measurement,
                               levels = c("T1", "T2", "T3", "T4"))
  
  # Rename column to 'sum'
  df_var <- df_var %>%
    rename(sum = all_of(col))
  
  # Average per participant per measurement round
  df_avg <- df_var %>%
    group_by(Castor.ID, measurement) %>%
    summarise(sum = mean(sum, na.rm = TRUE), .groups = "drop")
  
  # Summary statistics
  df_avg %>%
    group_by(measurement) %>%
    get_summary_stats(sum, type = "mean_sd") %>%
    print()
  
  # Boxplot
  bxp <- ggboxplot(df_avg, x = "measurement", y = "sum", add = "point")
  print(bxp)
  
  # Shapiro-Wilk normality test
  shapiro_res <- df_avg %>%
    group_by(measurement) %>%
    rstatix::shapiro_test(sum)
  
  print(shapiro_res)
  
  # Check data sufficiency
  if (n_distinct(df_avg$measurement) > 1 && n_distinct(df_avg$Castor.ID) > 1) {
    
    # Check if any measurement group violates normality
    if (any(shapiro_res$p < 0.05)) {
      print("Normality violated, running Friedman and Wilcoxon instead of ANOVA.")
      
      # Keep only participants with all 5 measurements
      df_avg_posthoc <- df_avg %>%
        group_by(Castor.ID) %>%
        filter(n_distinct(measurement) == 4) %>%
        ungroup()
      
      # Friedman test
      friedman_res <- rstatix::friedman_test(data = df_avg_posthoc, sum ~ measurement | Castor.ID)
      print(friedman_res)
      
      if (friedman_res$p > 0.05) {
        print("ANOVA not significant")
        plot_title <- paste("FRIEDMAN NOT SIG - Meting:", col) 
      } else {
        plot_title <- paste("FRIEDMAN SIG - Meting:", col) 
      }
      
      
      # Wilcoxon post-hoc
      pwc <- rstatix::pairwise_wilcox_test(
        data = df_avg_posthoc,
        formula = sum ~ measurement,
        paired = TRUE,
        p.adjust.method = "bonferroni"
      )
      
      pwc <- pwc %>% add_xy_position(x = "measurement")
      
      final_plot <- bxp +
        stat_pvalue_manual(pwc, step.increase = 0.07, tip.length = 0.01) +
        labs(title = plot_title,
             subtitle = rstatix::get_test_label(friedman_res, detailed = TRUE),
             caption = rstatix::get_pwc_label(pwc)) +
        xlab(NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      print(final_plot)
      
      fileName <- paste0("boxplot_Friedman_", col, ".svg")
      ggsave(paste0("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/Questionnaires/", fileName),
             plot = final_plot, width = 260, height = 200, units = "mm")
      
    } else {
      print("Normality not violated, running repeated measures ANOVA.")
      
      # Run repeated measures ANOVA
      res.aov <- rstatix::anova_test(data = df_avg, dv = sum, wid = Castor.ID, within = measurement)
      
      print(paste("ANOVA results for:", col))
      print(get_anova_table(res.aov))
      
      if (get_anova_table(res.aov)$p > 0.05) {
        print("ANOVA not significant")
        plot_title <- paste("ANOVA NOT SIG - Meting:", col) 
      } else {
        plot_title <- paste("ANOVA SIG - Meting:", col) 
      }
      
      df_avg_posthoc <- df_avg %>%
        group_by(Castor.ID) %>%
        filter(n_distinct(measurement) == 4) %>%
        ungroup()
      
      pwc <- rstatix::pairwise_t_test(data = df_avg_posthoc,
                                      formula = sum ~ measurement,
                                      paired = TRUE,
                                      p.adjust.method = "bonferroni")
      
      pwc <- pwc %>% add_xy_position(x = "measurement")
      
      final_plot <- bxp +
        stat_pvalue_manual(pwc, step.increase = 0.07, tip.length = 0.01) +
        labs(title = plot_title,
             subtitle = get_test_label(res.aov, detailed = TRUE),
             caption = get_pwc_label(pwc)) +
        xlab(NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      print(final_plot)
      
      fileName <- paste0("boxplot_RepeatedMeasuresAnova_", col, ".svg")
      ggsave(paste0("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/Questionnaires/", fileName),
             plot = final_plot, width = 260, height = 200, units = "mm")
    }
    
  } else {
    print(paste("Not enough data for:", col))
  }
}





# Repeat with T2 and T3 combined
# --------------------------------------- #

df_3measures <- df_sumscores %>%
  filter(measurement %in% c("T1", "T2", "T3", "T4")) %>%
  mutate(measurement_grouped = case_when(
    measurement == "T1" ~ "T1",
    measurement %in% c("T2", "T3") ~ "T2-T3",
    measurement == "T4" ~ "T4"
  )) %>%
  group_by(Castor.ID, measurement_grouped) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  rename(measurement = measurement_grouped)

for (col in cols_of_interest) {
  print(col)
  
  # Filter complete cases per variable
  df_var <- df_3measures %>%
    select(Castor.ID, measurement, all_of(col)) %>%
    filter(!is.na(.data[[col]]))
  
  # Set factor levels for measurement
  df_var$measurement <- factor(df_var$measurement,
                               levels = c("T1", "T2-T3", "T4"))
  
  # Rename column to 'sum'
  df_var <- df_var %>%
    rename(sum = all_of(col))
  
  # Summary statistics
  df_var %>%
    group_by(measurement) %>%
    get_summary_stats(sum, type = "mean_sd") %>%
    print()
  
  # Boxplot
  bxp <- ggboxplot(df_var, x = "measurement", y = "sum", add = "point")
  print(bxp)
  
  # Shapiro-Wilk normality test
  shapiro_res <- df_var %>%
    group_by(measurement) %>%
    rstatix::shapiro_test(sum)
  
  print(shapiro_res)
  
  # Check data sufficiency
  if (n_distinct(df_var$measurement) > 1 && n_distinct(df_var$Castor.ID) > 1) {
    
    # Check if any measurement group violates normality
    if (any(shapiro_res$p < 0.05)) {
      print("Normality violated, running Friedman and Wilcoxon instead of ANOVA.")
      
      # Keep only participants with all 5 measurements
      df_var_posthoc <- df_var %>%
        group_by(Castor.ID) %>%
        filter(n_distinct(measurement) == 3) %>%
        ungroup()
      
      # Friedman test
      friedman_res <- rstatix::friedman_test(data = df_var_posthoc, sum ~ measurement | Castor.ID)
      print(friedman_res)
      
      if (friedman_res$p > 0.05) {
        print("ANOVA not significant")
        plot_title <- paste("FRIEDMAN NOT SIG - Meting:", col) 
      } else {
        plot_title <- paste("FRIEDMAN SIG - Meting:", col) 
      }
      
      
      # Wilcoxon post-hoc
      pwc <- rstatix::pairwise_wilcox_test(
        data = df_var_posthoc,
        formula = sum ~ measurement,
        paired = TRUE,
        p.adjust.method = "bonferroni"
      )
      
      pwc <- pwc %>% add_xy_position(x = "measurement")
      
      final_plot <- bxp +
        stat_pvalue_manual(pwc, step.increase = 0.07, tip.length = 0.01) +
        labs(title = plot_title,
             subtitle = rstatix::get_test_label(friedman_res, detailed = TRUE),
             caption = rstatix::get_pwc_label(pwc)) +
        xlab(NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      print(final_plot)
      
      fileName <- paste0("boxplot_Friedman_", col, "_3measures.svg")
      ggsave(paste0("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/Questionnaires/", fileName),
             plot = final_plot, width = 260, height = 200, units = "mm")
      
    } else {
      print("Normality not violated, running repeated measures ANOVA.")
      
      # Run repeated measures ANOVA
      res.aov <- rstatix::anova_test(data = df_var, dv = sum, wid = Castor.ID, within = measurement)
      
      print(paste("ANOVA results for:", col))
      print(get_anova_table(res.aov))
      
      if (get_anova_table(res.aov)$p > 0.05) {
        print("ANOVA not significant")
        plot_title <- paste("ANOVA NOT SIG - Meting:", col) 
      } else {
        plot_title <- paste("ANOVA SIG - Meting:", col) 
      }
      
      df_var_posthoc <- df_var %>%
        group_by(Castor.ID) %>%
        filter(n_distinct(measurement) == 3) %>%
        ungroup()
      
      pwc <- rstatix::pairwise_t_test(
        data = df_var_posthoc,
        formula = sum ~ measurement,
        paired = TRUE,
        p.adjust.method = "bonferroni"
      )
      
      pwc <- pwc %>% add_xy_position(x = "measurement")
      
      final_plot <- bxp +
        stat_pvalue_manual(pwc, step.increase = 0.07, tip.length = 0.01) +
        labs(title = plot_title,
             subtitle = get_test_label(res.aov, detailed = TRUE),
             caption = get_pwc_label(pwc)) +
        xlab(NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      print(final_plot)
      
      fileName <- paste0("boxplot_RepeatedMeasuresAnova_", col, "_3measures.svg")
      ggsave(paste0("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/Questionnaires/", fileName),
             plot = final_plot, width = 260, height = 200, units = "mm")
    }
    
  } else {
    print(paste("Not enough data for:", col))
  }
}




# Compare deltascores
# --------------------------------------- #

# Define comparisons
contrasts <- list("delta_1_2" = c("T1", "T2"),
                  "delta_3_4" = c("T3", "T4"),
                  "delta_1_2/3" = c("T1", "T2-T3"),
                  "delta_2/3_4" = c("T2-T3", "T4"))

all_variables_results <- list()

# add T2-T3
df_avg_T2T3 <- df_sumscores %>%
  filter(measurement %in% c("T2", "T3")) %>%
  group_by(Castor.ID) %>%
  filter(n_distinct(measurement) == 2) %>%
  summarise(across(cols_of_interest, mean, na.rm = TRUE)) %>%
  mutate(measurement = "T2-T3") %>%
  relocate(measurement, .after = Castor.ID)
df_sumscores_extended <- bind_rows(df_sumscores, df_avg_T2T3)


for (col in cols_of_interest) {
  print(paste("Analyzing variable:", col))
  
  # Filter complete cases for variable
  df_var <- df_sumscores_extended %>%
    select(Castor.ID, measurement, all_of(col)) %>%
    filter(!is.na(.data[[col]])) %>%
    rename(sum = all_of(col))
  
  # Create empty list for this variable's results
  results_list <- list()
  
  for (contrast_name in names(contrasts)) {
    meetings <- contrasts[[contrast_name]]
    
    # Filter to relevant measurements
    df_contrast <- df_var %>%
      filter(measurement %in% meetings)
    
    # Keep participants with both measurements
    df_complete <- df_contrast %>%
      group_by(Castor.ID) %>%
      filter(n_distinct(measurement) == 2) %>%
      ungroup()
    
    # Pivot to wide format
    df_wide <- df_complete %>%
      pivot_wider(names_from = measurement, values_from = sum)
    
    col1 <- meetings[1]
    col2 <- meetings[2]
    
    # Check if both measurement columns exist
    if (!all(c(col1, col2) %in% colnames(df_wide))) {
      warning(paste("Skipping", contrast_name, "for", col, "- not enough data"))
      next
    }
    
    # Calculate delta
    df_wide <- df_wide %>%
      mutate(delta = .data[[col2]] - .data[[col1]])
    
    # Shapiro-Wilk test
    if (nrow(df_wide) >= 3) {
      shap <- shapiro.test(df_wide$delta)
      
      if (shap$p.value > 0.05) {
        # Parametric
        test_res <- t.test(df_wide$delta, mu = 0, paired = FALSE)
        test_type <- "t-test"
      } else {
        # Non-parametric
        test_res <- wilcox.test(df_wide$delta, mu = 0, paired = FALSE, exact = FALSE)
        test_type <- "Wilcoxon"
      }
      
      result_df <- tidy(test_res) %>%
        mutate(
          contrast = contrast_name,
          variable = col,
          test = test_type
        )
      
      results_list[[contrast_name]] <- result_df
    }
  }
  
  # Combine and store per variable
  variable_results <- bind_rows(results_list)
  all_variables_results[[col]] <- variable_results
}

# Combine all variables
final_results <- bind_rows(all_variables_results)
final_results <- final_results %>%
  mutate(p.value.bonferroni = p.adjust(p.value, method = "bonferroni")) %>%
  select(variable, contrast, test, statistic, p.value, alternative, p.value.bonferroni, method)  %>%
  mutate(contrast = paste0(contrast, " vs. 0"))

# Save
write_xlsx(final_results, "F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/Questionnaires/delta_tests_all_variables.xlsx")
print(final_results)


