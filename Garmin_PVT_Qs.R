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
setwd("/Volumes/Onderzoek/Groep Geuze/25U-0078_PRESS/")


# ~/Networkshares/ of ~/Volumes/





# ------------------------------------------------------------------------------------------------ #
#                                   Data Collection and preparation
# ------------------------------------------------------------------------------------------------ #

# Read demography
# ----------------------------------- #
df_CRF <- read_delim("E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_export_20250514_copy.csv", 
                      delim = ";", quote = "\"", show_col_types = FALSE)
# recode -99 to NA
df_CRF <- df_CRF %>%
  mutate(
    across(where(is.numeric), ~ na_if(., -99)),
    across(where(is.character), ~ na_if(., "-99"))
  )

df_demo <- subset(df_CRF, select = c("Participant Id", "rank"))
colnames(df_demo) <- c("Castor.ID", "rank")
         

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
  mutate(
    measurement = case_when(
      Survey.Package.Name == "Meting 1 vragenlijsten" ~ "T0",
      TRUE ~ NA_character_),
    Castor.ID = Castor.Participant.ID) %>%
  #replace 999 and 99 with NA
  mutate(
    across(c("sleeponset_24h_1","timesawake_24h_1","T0_durationawakening_24h"), 
           ~ na_if(na_if(., 999), 99) ) )

# Repeat for the 4-day diary but only take the past 24h (last night) because responses of previous nights are mostly missing or guesses
df_slaapdagboek4 <- read_csv("E_ResearchData/2_ResearchData/1. Verwerkte data/Castor/PRESS_Slaapdagboek_4_nachten_export_20250514_clean.csv",show_col_types = FALSE)
df_slaapdagboek4$sleeponset_24h <- as.numeric(df_slaapdagboek4$sleeponset_24h)

df_slaapdagboek4 <- df_slaapdagboek4 %>%
  mutate(
    measurement = case_when(
      Survey.Package.Name == "Meting 2 vragenlijsten" ~ "T1",
      Survey.Package.Name == "Meting 3 vragenlijsten" ~ "T2",
      Survey.Package.Name == "Meting 4 vragenlijsten" ~ "T3",
      TRUE ~ NA_character_),
    Castor.ID = Castor.Participant.ID) %>%
  #replace 999 and 99 with NA
  mutate(
    across(c("sleeponset_24h","timesawake_24h","T123_duration_awakening_24h"), 
           ~ na_if(na_if(., 999), 99) ) )

# Keep only rows where at least one column has a non-NA and non-zero value
df_slaapdagboek4 <- df_slaapdagboek4 %>%
  filter(
    if_any(
      c(sleeponset_24h, timesawake_24h, T123_duration_awakening_24h),
      ~ !is.na(.)
    )
  )

# Change Survey Completed On dates where paper-and-pencil was used (T1/M2): use the date_24h
df_slaapdagboek1 <- df_slaapdagboek1 %>%
  mutate(
    Survey.Completed.On = as_datetime( dmy_hms(Survey.Completed.On)  ) )

df_slaapdagboek4 <- df_slaapdagboek4 %>%
  mutate(
    Survey.Completed.On = as_datetime( ifelse( month(dmy_hms(Survey.Completed.On)) == 5,
                                               dmy_hms(paste(date_24h,"14:14:14",sep=" ") ),
                                               dmy_hms(Survey.Completed.On) ) ) )
    

# Calculate total sleep time (TST)
tst1 <- df_slaapdagboek1 %>%
  mutate(
    # Convert sleep_date to Date (strip time part if it exists)
    wake_date = as_date(Survey.Completed.On),
    # Parse onset/wake as times, recode to 24h time where necessary
    sleep_onset_time =  if_else(hour(hms(sleeptime_24h_1)) > 6 & hour(hms(sleeptime_24h_1)) < 12, 
                                hms(sleeptime_24h_1) + hours(12), 
                                if_else(hour(hms(sleeptime_24h_1)) == 12, 
                                        hms(sleeptime_24h_1) - hours(12), 
                                        hms(sleeptime_24h_1))),
    diff_awake      = hms(awakeonset_24h_1) - hms(lastawakening_24h_1),
    wake_onset_time = if_else( abs(hour(diff_awake)) > 0, hms(awakeonset_24h_1), hms(lastawakening_24h_1) ),
    # Adjust date for after-midnight sleep onset or before-midnight awakening (e.g., 00:30 belongs to the next day)
    sleep_onset_date = if_else(hour(sleep_onset_time) < 12, wake_date, wake_date-1),
    wake_onset_date = if_else(hour(wake_onset_time) < 18, wake_date, wake_date-1),
    # Build full datetimes
    sleep_onset_dt = sleep_onset_date + seconds(period_to_seconds(sleep_onset_time)),
    wake_onset_dt  = wake_onset_date + seconds(period_to_seconds(wake_onset_time)),
    # Calculate total sleep time in minutes
    total_sleep_minutes = as.numeric(difftime(wake_onset_dt, sleep_onset_dt, units = "mins")) - 
      replace_na(sleeponset_24h_1,0) - replace_na(T0_durationawakening_24h,0) #when missing (NA), count 0 minutes before sleep onset and 0 minutes awake during
  ) %>%
  group_by(Castor.Participant.ID, Survey.Completed.On) %>%
  summarise(
    sleep_onset_dt = sleep_onset_dt,
    wake_onset_dt = wake_onset_dt,
    total_sleep_minutes = sum(total_sleep_minutes, na.rm = TRUE),
    total_sleep_hours = total_sleep_minutes / 60,
    .groups = "drop")

tst4 <- df_slaapdagboek4 %>%
  mutate(
    # Convert sleep_date to Date (strip time part if it exists) 
    wake_date = as_date(Survey.Completed.On),
    # Parse onset/wake as times, recode to 24h time where necessary
    sleep_onset_time =  if_else(hour(hms(sleeptime_24h)) > 6 & hour(hms(sleeptime_24h)) < 12, 
                                hms(sleeptime_24h) + hours(12), 
                                if_else(hour(hms(sleeptime_24h)) == 12, 
                                        hms(sleeptime_24h) - hours(12), 
                                        hms(sleeptime_24h))),
    diff_awake      = hms(awakeonset_24h) - hms(lastawakening_24h),
    wake_onset_time = if_else( abs(hour(diff_awake)) > 0, hms(awakeonset_24h), hms(lastawakening_24h) ),
    # Adjust date for after-midnight sleep onset or before-midnight awakening (e.g., 00:30 belongs to the next day)
    sleep_onset_date = if_else(hour(sleep_onset_time) < 12, wake_date, wake_date-1),
    wake_onset_date = if_else(hour(wake_onset_time) < 18, wake_date, wake_date-1),
    # Build full datetimes
    sleep_onset_dt = sleep_onset_date + seconds(period_to_seconds(sleep_onset_time)),
    wake_onset_dt  = wake_onset_date + seconds(period_to_seconds(wake_onset_time)),
    # Calculate total sleep time in minutes
    total_sleep_minutes = as.numeric(difftime(wake_onset_dt, sleep_onset_dt, units = "mins")) - 
      replace_na(sleeponset_24h, 0) - replace_na(T123_duration_awakening_24h, 0) #when missing (NA), count 0 minutes before sleep onset and 0 minutes awake during
  ) %>%
  group_by(Castor.Participant.ID, Survey.Completed.On) %>%
  summarise(
    wake_date = wake_date,
    sleep_onset_dt = sleep_onset_dt,
    wake_onset_dt = wake_onset_dt,
    total_sleep_minutes = sum(total_sleep_minutes, na.rm = TRUE),
    total_sleep_hours = total_sleep_minutes / 60,
    .groups = "drop")

#merge with original diary dataframes and merge diary 1 (T0) with diary4 (T1-T3)
df_slaapdagboek1 <- df_slaapdagboek1 %>%
  full_join(tst1, by = c("Castor.Participant.ID", "Survey.Completed.On"))
df_slaapdagboek4 <- df_slaapdagboek4 %>%
  full_join(tst4, by = c("Castor.Participant.ID", "Survey.Completed.On"))

# select subset of variables of interest
df_sd1_subset <- df_slaapdagboek1 %>% 
  select(Castor.Participant.ID, Survey.Package.Name, Survey.Completed.On, 
         date_24h_1, bedtime_24h_1, sleeptime_24h_1, sleeponset_24h_1, timesawake_24h_1, 
         T0_durationawakening_24h, lastawakening_24h_1, awakeonset_24h_1, sleeprating_24h_1, 
         countcaffeine_24h_1, lasttimecaffeine_24h_1, cafeine_night_1, countcaffeine_nighht_1, 
         measurement, Castor.ID, sleep_onset_dt, wake_onset_dt, total_sleep_minutes, 
         total_sleep_hours)
#rename to same names as in slaapdagboek4
colnames(df_sd1_subset) <- c("Castor.Participant.ID", "Survey.Package.Name", "Survey.Completed.On", 
                             "date_24h", "bedtime_24h", "sleeptime_24h", "sleeponset_24h", "timesawake_24h", 
                             "T123_duration_awakening_24h", "lastawakening_24h", "awakeonset_24h", "sleeprating_24h", 
                             "countcaffeine_24h", "lasttimecaffeine_24h", "cafeine_night", "countcaffeine_nighht", 
                             "measurement", "Castor.ID", "sleep_onset_dt", "wake_onset_dt", "total_sleep_minutes", 
                             "total_sleep_hours")
#repeat for slaapdagboek4
df_sd4_subset <- df_slaapdagboek4 %>% 
  select(Castor.Participant.ID, Survey.Package.Name, Survey.Completed.On, 
         date_24h, bedtime_24h, sleeptime_24h, sleeponset_24h, timesawake_24h, 
         T123_duration_awakening_24h, lastawakening_24h, awakeonset_24h, sleeprating_24h, 
         countcaffeine_24h, lasttimecaffeine_24h, cafeine_night, countcaffeine_nighht, 
         measurement, Castor.ID, sleep_onset_dt, wake_onset_dt, total_sleep_minutes, 
         total_sleep_hours)

# merge by binding
df_slaapdagboek_all <- rbind(df_sd1_subset,df_sd4_subset)

# get average by measurement
df_slaapdagboek_avg <- df_slaapdagboek_all %>%
  group_by(Castor.ID, measurement) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")



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
  select(Castor.ID, measurement, rt_mean, rt_sd) %>%
  filter(measurement != "T1") #skip T1 (M2) due to missing data for almost all participants

df_in_slaap_vallen_sum <- convert_measurement(df_in_slaap_vallen) %>%
  rename(sum_slaap = sum)

df_spanning_sum <- convert_measurement(df_spanning) %>%
  rename(sum_spanning = sum)

df_vermoeidheid_sum <- convert_measurement(df_vermoeidheid) %>%
  rename(sum_vermoeidheid = sum)

df_slaapdagboek <- df_slaapdagboek_all %>%
  select(Castor.ID, measurement, total_sleep_minutes, T123_duration_awakening_24h)

df_vas <- df_vas %>%
  select("Castor.ID", "measurement", "alertness", "irritability", "skittish")
  
df_garmin$measurement[df_garmin$measurement=="x - Oefening stilgelegd"] <- "T1"
df_garmin_avg$measurement[df_garmin_avg$measurement=="x - Oefening stilgelegd"] <- "T1"

# Merge all Qs with df_garmin_avg and slaapdagboek_avg
df_merged <- df_garmin_avg %>%
  full_join(df_demo, by = "Castor.ID") %>%
  full_join(df_pvt, by = c("Castor.ID", "measurement")) %>%
  full_join(df_in_slaap_vallen_sum, by = c("Castor.ID", "measurement")) %>%
  full_join(df_spanning_sum, by = c("Castor.ID", "measurement")) %>%
  full_join(df_vermoeidheid_sum, by = c("Castor.ID", "measurement")) %>%
  full_join(df_vas, by = c("Castor.ID", "measurement")) %>%
  full_join(df_slaapdagboek, by = c("Castor.ID", "measurement"))  

# Merge sleep diary with df_garmin
#df_slaapdagboek_all$Survey.Completed.On <- as.Date(dmy_hms(df_slaapdagboek_all$Survey.Completed.On))
df_slaapdagboek_all$Date <- as_date(ymd_hms(df_slaapdagboek_all$Survey.Completed.On))
df_garmin$Date <- as_date(df_garmin$Date )
df_sleep <- df_garmin %>%
  full_join(df_slaapdagboek_all, by = c("Castor.ID", "measurement","Date"))  

# Find difference in total sleep time between Garmin and Sleep diary
df_sleep$Avg.Sleep_Nap.min <- df_sleep$Avg.Sleep.min + 
  if_else(!is.na(df_sleep$Nap.Duration.min), df_sleep$Nap.Duration.min, 0)

garmin_diary_diff <- df_sleep %>%
  rowwise() %>%
  mutate(diff = if_all(c(Avg.Sleep_Nap.min, total_sleep_minutes), ~ !is.na(.)) * (Avg.Sleep_Nap.min - total_sleep_minutes)) %>%
  ungroup()
 
# plot the difference
garmin_diary_diff <- garmin_diary_diff %>%
  mutate(outcome_cat = case_when(
    is.na(diff)   ~ "missing",
    diff < -30    ~ "Garmin TST lager dan slaapdagboek TST (>30 min. verschil)",
    diff > 30     ~ "Garmin TST hoger dan slaapdagboek TST (>30 min. verschil)",
    diff > -29 & diff < 29 ~ "Garmin TST ~ slaapdagboek TST (<30 min. verschil)"))

plotdiff <- subset(garmin_diary_diff, Date == "2025-03-27" | Date == "2025-04-10" | Date == "2025-04-15" | Date == "2025-04-24" )

# Calculate percentages per Date + outcome_cat
plotdiff_perc <- plotdiff %>%
  group_by(measurement, outcome_cat) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(measurement) %>%
  mutate(perc = count / sum(count))

ggplot(plotdiff_perc, aes(x = measurement, y = perc, fill = outcome_cat)) +
  geom_col() +
  geom_text(aes(label = scales::percent(perc, accuracy = 1)), 
            position = position_stack(vjust = 0.5), color = "white", size = 3) +
  labs(x = "Meetmoment", y = "Percentage", fill = "Verschil in totale slaaptijd (TST)") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("darkgreen", "darkorange", "purple", "lightgrey")) +
  theme_minimal()

ggplot(plotdiff[!is.na(plotdiff$diff), ], aes(x = measurement, y = diff, color = outcome_cat)) +
  geom_jitter(size = 3, width = 0.1, alpha = 0.7) +
  labs(x = "Meetmoment", y = "Verschil in minuten (Garmin TST - slaapdagboek TST)") +
  scale_color_manual(values = c("darkgreen", "darkorange", "purple", "lightgrey")) +
  theme_minimal()



# ------------------------------------------------------------------------------------------------ #
#                               Add T2/T3 and delta's
# ------------------------------------------------------------------------------------------------ #

# Combine T1 and T2 into one measurement
# -------------------------------------- #
df_t1_t2 <- df_merged %>%
  filter(measurement %in% c("T1", "T2")) %>%
  group_by(Castor.ID) %>%
  summarise(across(
    .cols = where(is.numeric),
    .fns  = ~ {
      m <- mean(., na.rm = TRUE)
      if (is.nan(m)) NA_real_ else m
    }
  ), .groups = "drop") %>%
  mutate(measurement = "T1/2") %>%
  select(all_of(colnames(df_merged)))

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
    values_from = c(HRV.Last.Night.Avg...ms., Avg.Sleep.min, Avg.Sleep.Awake.Dur.min, 
                    Sleep.Score, rt_mean, rt_sd, sum_slaap, sum_spanning, sum_vermoeidheid, 
                    alertness, irritability, skittish, total_sleep_minutes, T123_duration_awakening_24h),
    names_sep = "_"
  ) %>%
  mutate(across(
    ends_with("_T1/2"),
    ~ ifelse(is.na(.x) | is.na(get(sub("_T1/2$", "_T0", cur_column()))), NA_real_, .x - get(sub("_T1/2$", "_T0", cur_column()))),
    .names = "{sub('_T1/2$', '', .col)}"
  )) %>%
  select(Castor.ID, all_of(c("HRV.Last.Night.Avg...ms.", "Avg.Sleep.min", 
                             "Avg.Sleep.Awake.Dur.min", "Sleep.Score", 
                             "rt_mean", "rt_sd", "sum_slaap", 
                             "sum_spanning", "sum_vermoeidheid", 
                             "alertness", "irritability", "skittish", 
                             "total_sleep_minutes","T123_duration_awakening_24h"))) %>%
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
vars_to_plot <- c("HRV.Last.Night.Avg...ms.","Avg.Sleep.min", "sum_vermoeidheid", 
                  "sum_slaap", "sum_spanning", "alertness", "irritability", "skittish", 
                  "rt_mean", "rt_sd","total_sleep_minutes", "T123_duration_awakening_24h")

df_plot <- df_merged %>%
  filter(measurement %in% c("T0", "T1/2", "T3")) %>%
  mutate(
    measurement = factor(measurement, levels = c("T0", "T1/2", "T3")),
    sleep_cat = case_when(
      total_sleep_minutes < 300  ~ "Slaap: minder dan 5h",
      Avg.Sleep.min < 300        ~ "Slaap: minder dan 5h",
      total_sleep_minutes >= 300 ~ "Slaap: meer dan 5h",
      Avg.Sleep.min >= 300       ~ "Slaap: meer dan 5h",
      TRUE ~ NA_character_
    ),
    rank_cat = case_when(
      rank > 1  ~ "(onder)officier",
      rank == 1 ~ "manschap",
      is.na(rank) ~ "manschap",
      TRUE ~ NA_character_
    )
  )
varname <-   "alertness"  "irritability", "skittish", 
"rt_mean", "rt_sd","total_sleep_minutes", "T123_duration_awakening_24h"
  
for (varname in vars_to_plot) {
  # # Plot met facets per deelnemer
  # p_facet <- ggplot(df_plot, aes(x = measurement, y = .data[[varname]], group = Castor.ID)) +
  #   geom_line(alpha = 0.7) +
  #   geom_point(alpha = 0.7) +
  #   facet_wrap(~ Castor.ID, scales = "free_y") +
  #   labs(
  #     title = paste("Verloop van", varname, "per deelnemer (facets)"),
  #     x = "Measurement",
  #     y = varname
  #   ) +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 
  # fileName <- paste0("verloop_", varname, ".svg")
  # ggsave(paste0("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/", fileName),
  #        plot = p_facet, width = 260, height = 200, units = "mm")
  
  # Plot met alle deelnemers in 1 figuur zonder facets
  p_all_sleep <- ggplot(df_plot, aes(x = measurement, y = .data[[varname]], group = Castor.ID, color = factor(sleep_cat))) +
    geom_line(alpha = 0.7, color="black", linewidth=0.2) +
    geom_point(alpha = 0.7, size = 3) +
    labs(
      title = "Color split: Slaapduur",
      x = "Meetmoment",
      y = varname,
      color = "Slaap"
    ) +
    scale_color_manual(values = c("darkgreen", "darkorange", "lightgrey")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p_all_rank <- ggplot(df_plot, aes(x = measurement, y = .data[[varname]], group = Castor.ID, color = factor(rank_cat))) +
    geom_line(alpha = 0.7, color="black", linewidth=0.2) +
    geom_point(alpha = 0.7, size = 3) +
    labs(
      title = "Color split: Rang (categorie)",
      x = "Meetmoment",
      y = varname,
      color = "Rang_categorie"
    ) +
    scale_color_manual(values = c( "brown", "darkblue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  fileName <- paste0("verloop_alleDeelnemers", varname, ".svg")
  ggsave(paste0("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/", fileName),
         plot = p_all_sleep, width = 260, height = 200, units = "mm")
  ggsave(paste0("F_DataAnalysis/1f_DataAnalysisScripts/Outcomes Xandra/", fileName),
         plot = p_all_rank, width = 260, height = 200, units = "mm")
  
}


ggplot(df_plot, aes(x = measurement, y = Avg.Sleep.min, group = Castor.ID, color = factor(sleep_cat))) +
  geom_jitter()

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












