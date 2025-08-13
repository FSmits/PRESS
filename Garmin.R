##----------------------------------------------------------------------------------------------- ##
#                                   GARMIN DATA
##----------------------------------------------------------------------------------------------- ##
#
# Description:  This script cleans and describes the Garmin data.
#
# Date:         June 2025
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

# visualizations
library(ggplot2)
library(ggpubr)

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

# Read garmin code
# ----------------------------------- #
df_garmin_raw <- read_excel("E_ResearchData/2_ResearchData/1. Verwerkte data/Garmin/Garmin_data.xlsx")


# Read castor code
# ----------------------------------- #
df_codelist <- data.frame(read_excel("C_PersonalData/2_CodeList/Sleutelbestand_PRESS.xlsx"))  %>%
  select("Castor.ID", "Garmin.ID")  %>%
  mutate(Garmin.ID = as.numeric(Garmin.ID))


# Prepare data
# ----------------------------------- #
df_garmin <- data.frame(df_garmin_raw) %>%
  filter(Name != "All Athletes") %>%
  distinct(Date, Name, .keep_all = TRUE) %>% #remove duplicate 68_7CT006576, keep first line (second had NA)
  mutate(Name = str_replace(Name, "(?<!_)7[cC]", "_7c"),
         Code = as.numeric(str_extract(Name, "^\\d+")),
         Avg.Sleep..h.min. = as.numeric(difftime(Avg.Sleep..h.min., as.POSIXct("1899-12-31 00:00:00"), units = "mins")),
         Nap.Duration..h.min. = as.numeric(difftime(Nap.Duration..h.min., as.POSIXct("1899-12-31 00:00:00"), units = "mins")),
         Avg..Sleep.Awake.Dur...h.min. = as.numeric(difftime(Avg..Sleep.Awake.Dur...h.min., as.POSIXct("1899-12-31 00:00:00"), units = "mins"))) %>%
  rename (Avg.Sleep.min = Avg.Sleep..h.min.,
          Nap.Duration.min = Nap.Duration..h.min.,
          Avg.Sleep.Awake.Dur.min = Avg..Sleep.Awake.Dur...h.min.) %>%
  left_join(df_codelist, by = c("Code" = "Garmin.ID")) %>% # merge with castor ID
  select(-Code) %>%
  relocate(Castor.ID, .after = 2) %>%
  mutate(measurement = case_when(
    Date >= as.Date("2025-03-28") & Date <= as.Date("2025-04-06") ~ "T0",
    Date >= as.Date("2025-04-07") & Date <= as.Date("2025-04-09") ~ "T1",
    Date == as.Date("2025-04-10")                                 ~ "x - Oefening stilgelegd",
    Date == as.Date("2025-04-11")                                 ~ "x - Transitiedag",
    Date >= as.Date("2025-04-12") & Date <= as.Date("2025-04-15") ~ "T2",
    Date >= as.Date("2025-04-16") & Date <= as.Date("2025-04-24") ~ "T3",
    TRUE ~ NA_character_)) %>%  # add measurement column based on date
  relocate(measurement, .after = 3)

write.csv(df_garmin,"E_ResearchData/2_ResearchData/1. Verwerkte data/Garmin/Garmin_data_clean2.csv", row.names = FALSE)

# ------------------------------------------------------------------------------------------------ #
#                             Measurement completion rate
# ------------------------------------------------------------------------------------------------ #

cols_of_interest <- c("HRV.Last.Night.Avg...ms.",
                      "Resting.HR..bpm.",
                      "Avg.Sleep.min",
                      "Nap.Duration.min",
                      "HRV.Last.Night.High..ms.",
                      "Min..Daily.HR..bpm.",
                      "Avg.Sleep.Awake.Dur.min",
                      "Sleep.Score")

for (col in cols_of_interest){
  df <- df_garmin %>%
    select("Castor.ID", "measurement", "Date", cols_of_interest[col])
  
  df$Date <- as.Date(df$Date)
  
  # Dynamically compute non-NA percentage for the current column
  completion_pct_df <- df %>%
    group_by(Date, measurement) %>%
    summarise(
      total = n(),
      ingevuld = sum(!is.na(.data[[cols_of_interest[col]]])),
      completion_pct = ingevuld / total * 100,
      .groups = "drop"
    )
  
  # Average per measurement
  avg_completion_df <- completion_pct_df %>%
    group_by(measurement) %>%
    summarise(avg_completion = mean(completion_pct), .groups = "drop")
  
  # Plot
  final_plot <- ggplot(completion_pct_df, aes(x = Date, y = completion_pct, fill = measurement)) +
    geom_col(position = "dodge") +
    geom_hline(data = avg_completion_df,
               aes(yintercept = avg_completion, color = measurement),
               linetype = "dashed", linewidth = 0.6) +
    labs(
      title = paste("Completion percentage for", cols_of_interest[col]),
      subtitle = "Dashed line = average per measurement point",
      x = "Date",
      y = "Completion (%)",
      fill = "Measurement",
      color = "Average"
    ) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    theme_minimal()
  
  fileName <- paste("CompletionRate_", col, ".svg", sep = "")
  ggsave(paste("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/Garmin/", fileName, sep = ""), plot = final_plot,
         width = 260, height = 200, units = "mm")
}




# ------------------------------------------------------------------------------------------------ #
#                            Assessment of scores over time
# ------------------------------------------------------------------------------------------------ #


# Repeated measures ANOVA/Friedman
# --------------------------------------- #

for (col in cols_of_interest) {
  print(col)
  
  # Filter complete cases per variable
  df_var <- df_garmin %>%
    select(Castor.ID, measurement, all_of(col)) %>%
    filter(!is.na(.data[[col]]))
  
  # Set factor levels for measurement
  df_var$measurement <- factor(df_var$measurement,
                               levels = c("T0","T1", "T2", "T3", "x - Oefening stilgelegd", "x - Transitiedag",))
  
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
        filter(n_distinct(measurement) == 5) %>%
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
        ggsave(paste0("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/Garmin/", fileName),
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
          filter(n_distinct(measurement) == 5) %>%
          ungroup()
        
        pwc <- rstatix::pairwise_t_test(
          data = df_avg_posthoc,
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
        
        fileName <- paste0("boxplot_RepeatedMeasuresAnova_", col, ".svg")
        ggsave(paste0("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/Garmin/", fileName),
               plot = final_plot, width = 260, height = 200, units = "mm")
    }
    
  } else {
    print(paste("Not enough data for:", col))
  }
}





# Repeat with T2 and T3 combined
# --------------------------------------- #

df_garmin_3measures <- df_garmin %>%
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
  df_var <- df_garmin_3measures %>%
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
      ggsave(paste0("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/Garmin/", fileName),
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
      ggsave(paste0("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/Garmin/", fileName),
             plot = final_plot, width = 260, height = 200, units = "mm")
    }
    
  } else {
    print(paste("Not enough data for:", col))
  }
}




