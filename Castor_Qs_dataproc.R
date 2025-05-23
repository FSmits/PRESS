##----------------------------------------------------------------------------------------------- ##
#                                      RELIABILITY ASSESSMENT
##----------------------------------------------------------------------------------------------- ##
#
# Description:  R-script to processe questionnaire data, compute their reliabilies and analyse differences between timepoints.
#             - Study name: PRESS (dossiernr.: 25U-0078)
#             - Measure/instrument: Likert-scale based questionnaires: State and Trait Anxiety Inventory short form (STAI-6, 'Spanning'), Pre-sleep arousal scale (PSAS, 'In slaap vallen', Patient-Reported Outcomes Measurement Information System – fatigue subscale (PROMIS, 'Vermoeidheid')
#             - Data type: Self-report, online (M1, M3, M4) or paper-pencil (M2)
#             - Design: within-subjects longitudinal, 4 timepoints (M1: 27th 03/'25, M2: 10th and 11th 04/'25, M3: 15th 04/'25,, M4: 24th 04/'25)
#
# Date:         May 2025
# R.version:    4.4.0 (2024-04-24)
# Rstudio:      2025.05.0+496
#
## ---------------------------------------------------------------------------------------------- ##

# clear environment
rm(list=ls())

# ------------------------------------------------------------------------------------------------ #
# ----  Settings & Dependencies ----
# ------------------------------------------------------------------------------------------------ #

# ---- Required packages ----
# read files
library(readxl)
library(writexl)
# data manipulation
library(tidyverse)
library(dplyr)
# statistics
library(ltm)
library(psych)
library(rstatix)
# visualizations
library(ggplot2)
library(ggpubr)

# ----- Functions ------



# ----- working directory ------
# external volume, change value if necessary 
#setwd(paste("~/Networkshares/", heronderzoek, "/Groep Geuze/25U-0078_PRESS/", sep = ""))
setwd(paste("/Volumes/heronderzoek-6/Groep Geuze/25U-0078_PRESS/", sep = ""))

# ~/Networkshares/ of ~/Volumes/



# ------------------------------------------------------------------------------------------------ #
# ----  Data Processing ----
# ------------------------------------------------------------------------------------------------ #

# Data In slaap vallen (PSAS) 
# ---------------------------- #
df_in_slaap_vallen <- data.frame(read_delim("E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_In_slaap_vallen_export_20250514.csv", delim = ";", show_col_types = FALSE))
df_in_slaap_vallen <- df_in_slaap_vallen %>%
  filter(Survey.Progress == 100) %>%
  dplyr::select("Castor.Participant.ID", "Survey.Package.Name", starts_with("PSAS_")) %>%
  mutate(sum = rowSums(dplyr::select(., starts_with("PSAS_")), na.rm = TRUE))

# Data Spanning (STAI-6) 
# ---------------------- #
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

# Data Vermoeidheid (PROMIS) 
# -------------------------- #
df_vermoeidheid <- data.frame(read_delim("E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_Vermoeidheid_export_20250514.csv", delim = ";", show_col_types = FALSE))
df_vermoeidheid <- df_vermoeidheid %>%
  filter(Survey.Progress == 100) %>%
  dplyr::select("Castor.Participant.ID", "Survey.Package.Name", starts_with("PROMIS_")) %>%
  mutate(sum = rowSums(dplyr::select(., starts_with("PROMIS_")), na.rm = TRUE))



# ------------------------------------------------------------------------------------------------ #
# ----  Data Analysis ----
# ------------------------------------------------------------------------------------------------ #

# -------  Cronbach's alpha ---------
# ----------------------------------- #
# In slaap vallen
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



# ----- Test-hertest reliability ------- 
# -------------------------------------- #
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


# ---- Repeated measures ANOVA M1, M2, M3, M4 ----
# ------------------------------------------------ #
df %>%
  group_by(Survey.Package.Name) %>%
  get_summary_stats(sum, type = "mean_sd")

# Boxplot
bxp <- ggboxplot(df, x = "Survey.Package.Name", y = "sum", add = "point")
bxp

# Shapiro
df %>%
  group_by(Survey.Package.Name) %>%
  shapiro_test(sum)

# Anova
res.aov <- anova_test(data = df, dv = sum, wid = Castor.Participant.ID, within = Survey.Package.Name)
get_anova_table(res.aov)

# Pairwise comparisons
pwc <- df_volledige_metingen %>%
  pairwise_t_test(sum ~ Survey.Package.Name, paired = TRUE,
                  p.adjust.method = "bonferroni")
pwc

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Survey.Package.Name")
final_plot <- bxp + 
  stat_pvalue_manual(pwc, step.increase = 0.07) +
  labs(title = paste("Vragenlijst:", vragenlijst),
       subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc)) +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Save figure
fileName <- paste("boxplot_repeatedMeasuresAnova_", vragenlijst, ".svg", sep = "")
ggsave(paste("F_DataAnalysis/1f_DataAnalysisScripts/", fileName, sep = ""), plot = final_plot,
       width = 260, height = 200, units = "mm")



# ---- Repeated measures ANOVA (herhaal met gemiddelde van M2 en M3) ----
# ----------------------------------------------------------------------- #
# Stap 1: filter alleen meting 2 en 3
df_23 <- df %>%
  filter(Survey.Package.Name %in% c("Meting 2 vragenlijsten", "Meting 3 vragenlijsten"))

# Stap 2: bereken gemiddelde per ID
df_23_avg <- df_23 %>%
  group_by(Castor.Participant.ID) %>%
  filter(n() == 2) %>%  # alleen deelnemers met beide metingen
  summarise(Survey.Package.Name = "Meting 2-3 gemiddelde",
            sum = mean(sum, na.rm = TRUE),
            .groups = "drop")

# Stap 3: bind deze rijen aan originele dataframe
df_3_metingen <- bind_rows(df, df_23_avg)

# Bewaar alleen meting 1, 2+3 en 4
df_3_metingen <- df_3_metingen %>%
  filter(Survey.Package.Name %in% c("Meting 1 vragenlijsten", 
                                    "Meting 2-3 gemiddelde", 
                                    "Meting 4 vragenlijsten"))
df_3_metingen$Survey.Package.Name <- factor(df_3_metingen$Survey.Package.Name,
                                            levels = c("Meting 1 vragenlijsten", "Meting 2-3 gemiddelde", "Meting 4 vragenlijsten"))

# # ALTERNATIEF voor PLOT: verwijder M2 (behoud alleen M2 van tijdens-oefening)
# df_3_metingen <- df[-which(df$Survey.Package.Name=="Meting 3 vragenlijsten"),] 

sumstats <- df_3_metingen %>%
  group_by(Survey.Package.Name) %>%
  get_summary_stats(sum, type = "mean_sd")
sumstats


# Boxplot
bxp <- ggboxplot(df_3_metingen, x = "Survey.Package.Name", y = "sum", add = "point")
bxp

# Shapiro
df_3_metingen %>%
  group_by(Survey.Package.Name) %>%
  shapiro_test(sum)


# Anova
res.aov <- anova_test(data = df_3_metingen, dv = sum, wid = Castor.Participant.ID, within = Survey.Package.Name)
get_anova_table(res.aov)


# Pairwise comparisons
df_3_metingen_complete <- df_3_metingen %>%
  group_by(Castor.Participant.ID) %>%
  filter(n_distinct(Survey.Package.Name) == 3) %>%
  ungroup()

pwc <- df_3_metingen_complete %>%
  pairwise_t_test(sum ~ Survey.Package.Name, paired = TRUE,
                  p.adjust.method = "bonferroni")
pwc


# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Survey.Package.Name")
final_plot <- bxp + 
  stat_pvalue_manual(pwc, step.increase = 0.07) +
  labs(title = paste("Vragenlijst:", vragenlijst),
       subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc)) +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Save figure
fileName <- paste("boxplot_repeatedMeasuresAnova_", vragenlijst, "_3momenten.svg", sep = "")
ggsave(paste("F_DataAnalysis/1f_DataAnalysisScripts/", fileName, sep = ""), plot = final_plot,
       width = 260, height = 200, units = "mm")

# Visualisation: simple bar plot
showbarplot <- ggplot(sumstats) +
  geom_bar( aes(x=Survey.Package.Name, y=mean), stat="identity", fill="#136497", alpha=0.8) +
  geom_point(aes(x=Survey.Package.Name, y=mean), size=3) +
  geom_errorbar( aes(x=Survey.Package.Name, ymin=mean-sd, ymax=mean+sd), width=0.2, colour="black", alpha=0.9, size=0.3) +
  ggtitle("Vragenlijst: Vermoeidheid") +
  xlab("Meetmoment") +
  ylab("Gemiddelde somscore") +
  theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_line(colour = "grey"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 




# ---- Repeated measures ANOVA (herhaal met verschilscores) ----
# -------------------------------------------------------------- #
results_list <- list()
results_list_w <- list() # for wilcox results


# verschil meting 1 en 2
df_12 <- df %>%
  filter(Survey.Package.Name %in% c("Meting 1 vragenlijsten", "Meting 2 vragenlijsten"))

df_12_complete <- df_12 %>%
  group_by(Castor.Participant.ID) %>%
  filter(n_distinct(Survey.Package.Name) == 2) %>%
  ungroup()

df_wide <- df_12_complete %>%
  select(Castor.Participant.ID, Survey.Package.Name, sum) %>%
  pivot_wider(names_from = Survey.Package.Name, values_from = sum,
              names_prefix = "meting_")  # kolommen: meting_Meting_1, meting_Meting_2

df_wide <- df_wide %>%
  rename(meting1 = 'meting_Meting 1 vragenlijsten',
         meting2 = 'meting_Meting 2 vragenlijsten')

df_wide <- df_wide %>%
  mutate(verschil = meting2 - meting1)

ttest_diff_12 <- t.test(df_wide$verschil, mu = 0, paired = FALSE)
results_list[["Verschil 1-2 vs 0"]] <- tidy(ttest_diff_12)

wilcox_diff_12 <- wilcox.test(df_wide$verschil, mu = 0, paired = FALSE, exact = FALSE)
results_list_w[["Verschil 1-2 vs 0"]] <- tidy(wilcox_diff_12)

ttest_12 <- t.test(df_wide$meting1, df_wide$meting2, paired = TRUE)
results_list[["Meting 1 vs 2"]] <- tidy(ttest_12)

wilcox_12 <- wilcox.test(df_wide$meting1, df_wide$meting2, paired = TRUE, exact = FALSE)
results_list_w[["Meting 1 vs 2"]] <- tidy(wilcox_12)


# Verschil meting 3 en 4
df_34 <- df %>%
  filter(Survey.Package.Name %in% c("Meting 3 vragenlijsten", "Meting 4 vragenlijsten"))

df_34_complete <- df_34 %>%
  group_by(Castor.Participant.ID) %>%
  filter(n_distinct(Survey.Package.Name) == 2) %>%
  ungroup()

df_wide <- df_34_complete %>%
  select(Castor.Participant.ID, Survey.Package.Name, sum) %>%
  pivot_wider(names_from = Survey.Package.Name, values_from = sum,
              names_prefix = "meting_")  # kolommen: meting_Meting_3, meting_Meting_4

df_wide <- df_wide %>%
  rename(meting3 = 'meting_Meting 3 vragenlijsten',
         meting4 = 'meting_Meting 4 vragenlijsten')

df_wide <- df_wide %>%
  mutate(verschil = meting4 - meting3)

ttest_diff_34 <- t.test(df_wide$verschil, mu = 0, paired = FALSE)
results_list[["Verschil 3-4 vs 0"]] <- tidy(ttest_diff_34)

wilcox_diff_34 <- wilcox.test(df_wide$verschil, mu = 0, paired = FALSE, exact = FALSE)
results_list_w[["Verschil 3-4 vs 0"]] <- tidy(wilcox_diff_34)

ttest_34 <- t.test(df_wide$meting3, df_wide$meting4, paired = TRUE)
results_list[["Meting 3 vs 4"]] <- tidy(ttest_34)

wilcox_34 <- wilcox.test(df_wide$meting3, df_wide$meting4, paired = TRUE, exact = FALSE)
results_list_w[["Meting 3 vs 4"]] <- tidy(wilcox_34)



# Verschil meting 1 en 2/3
df_123 <- df_3_metingen %>%
  filter(Survey.Package.Name %in% c("Meting 1 vragenlijsten", "Meting 2-3 gemiddelde"))

df_123_complete <- df_123 %>%
  group_by(Castor.Participant.ID) %>%
  filter(n_distinct(Survey.Package.Name) == 2) %>%
  ungroup()

df_wide <- df_123_complete %>%
  select(Castor.Participant.ID, Survey.Package.Name, sum) %>%
  pivot_wider(names_from = Survey.Package.Name, values_from = sum,
              names_prefix = "meting_")  # kolommen: meting_Meting_1, meting_Meting_2-3

df_wide <- df_wide %>%
  rename(meting1 = 'meting_Meting 1 vragenlijsten',
         meting23 = 'meting_Meting 2-3 gemiddelde')

df_wide <- df_wide %>%
  mutate(verschil = meting23 - meting1)

ttest_diff_123 <- t.test(df_wide$verschil, mu = 0, paired = FALSE)
results_list[["Verschil 1-2/3 vs 0"]] <- tidy(ttest_diff_123)

wilcox_diff_123 <- wilcox.test(df_wide$verschil, mu = 0, paired = FALSE, exact = FALSE)
results_list_w[["Verschil 1-2/3 vs 0"]] <- tidy(wilcox_diff_123)

ttest_123 <- t.test(df_wide$meting1, df_wide$meting23, paired = TRUE)
results_list[["Meting 1 vs 2-3"]] <- tidy(ttest_123)

wilcox_123 <- wilcox.test(df_wide$meting1, df_wide$meting23, paired = TRUE, exact = FALSE)
results_list_w[["Meting 1 vs 2-3"]] <- tidy(wilcox_123)



# Verschil 2/3 en 4
df_234 <- df_3_metingen %>%
  filter(Survey.Package.Name %in% c("Meting 2-3 gemiddelde", "Meting 4 vragenlijsten"))

df_234_complete <- df_234 %>%
  group_by(Castor.Participant.ID) %>%
  filter(n_distinct(Survey.Package.Name) == 2) %>%
  ungroup()

df_wide <- df_234_complete %>%
  select(Castor.Participant.ID, Survey.Package.Name, sum) %>%
  pivot_wider(names_from = Survey.Package.Name, values_from = sum,
              names_prefix = "meting_")  # kolommen: meting_Meting_2-3, meting_Meting_4

df_wide <- df_wide %>%
  rename(meting23 = 'meting_Meting 2-3 gemiddelde',
         meting4 = 'meting_Meting 4 vragenlijsten')

df_wide <- df_wide %>%
  mutate(verschil = meting23 - meting4)

ttest_diff_234 <- t.test(df_wide$verschil, mu = 0, paired = FALSE)
results_list[["Verschil 2/3-4 vs 0"]] <- tidy(ttest_diff_234)

wilcox_diff_234 <- wilcox.test(df_wide$verschil, mu = 0, paired = FALSE, exact = FALSE)
results_list_w[["Verschil 2/3-4 vs 0"]] <- tidy(wilcox_diff_234)

ttest_234 <- t.test(df_wide$meting23, df_wide$meting4, paired = TRUE)
results_list[["Meting 2-3 vs 4"]] <- tidy(ttest_234)

wilcox_234 <- wilcox.test(df_wide$meting23, df_wide$meting4, paired = TRUE, exact = FALSE)
results_list_w[["Meting 2-3 vs 4"]] <- tidy(wilcox_234)



# ------------------------------------------------------------------------------------------------ #
#   Save outcomes
# ------------------------------------------------------------------------------------------------ #

# Resultaten in dataframe
t_test_results <- bind_rows(results_list, .id = "Vergelijking")
wilcox_test_results <- bind_rows(results_list_w, .id = "Vergelijking")

# Voeg Bonferroni-correctie toe
t_test_results <- t_test_results %>%
  mutate(p.value.bonferroni = p.adjust(p.value, method = "bonferroni"))
wilcox_test_results <- wilcox_test_results %>%
  mutate(p.value.bonferroni = p.adjust(p.value, method = "bonferroni"))

# Print en save overzicht
print(t_test_results)
write_xlsx(t_test_results, paste("F_DataAnalysis/1f_DataAnalysisScripts/t_test_results_verschilscores_", vragenlijst, ".xlsx", sep = ""))
print(wilcox_test_results)
write_xlsx(wilcox_test_results, paste("F_DataAnalysis/1f_DataAnalysisScripts/wilcox_test_results_verschilscores_", vragenlijst, ".xlsx", sep = ""))
