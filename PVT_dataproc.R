# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
#
# Description: R-script to read and process PVT data
#             - Study name: PRESS (dossiernr.: 25U-0078)
#             - Measure/instrument: smartphone-based Psychommotor Vigilance Task (PVT)
#             - Data type: behavioral performance (reaction times and accuracy)
#             - Design: within-subjects longitudinal, 4 timepoints (M1: 27th 03/'25, M2: 10th and 11th 04/'25, M3: 15th 04/'25,, M4: 24th 04/'25)
#
# Date:         May 2025
# R.version:    4.4.0 (2024-04-24) -- "Puppy Cup"
# Rstudio:      2025.05.0+496
#
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

# clear environment
rm(list=ls())


# ------------------------------------------------------------------------------------------------ #
# ----    Settings & Dependencies ----
# ------------------------------------------------------------------------------------------------ #

# ---- Required packages ----
require("reshape2")
require("ggplot2")
require("ggpubr")
require("ggstatsplot") #install.packages("ggstatsplot")
require("tidyverse")
require("dplyr")
require("plyr")
require("rstatix")

# ----- Functions ------
pvt_func <- function(filepath){
  # load the raw data file
  raw <- read.csv(filepath)
  
  # replace NAN by NA
  raw[raw=="NAN"] <- NA
  
  # extract participant number (PVT-ID) and date and time
  PVTid    <- raw$participant[1]
  datetime  <- raw$date[1]
  
  # false responses: count cells with a value in a 'dont_respond' column
  falseresponse_trls <- which(!is.na(raw$dont_respond.x))
  falseresp          <- length(falseresponse_trls)
  # lapses (misses/no response): count cells with a realistic difference (>60 ms) between 'Target_started' and 'Target_stopped' but without a value in both the 'dont_respond' and 'response' columns
  lapses_trls <- which( (raw$Target.stopped-raw$Target.started)>.060 & is.na(raw$dont_respond.x) & is.na(raw$response.x))
  lapses      <- length(lapses_trls)
  # calculate mean and standard deviation of the reaction time
  rt_mean <- mean( as.numeric(raw$RTms), na.rm=TRUE)
  rt_sd   <- sd( as.numeric(raw$RTms), na.rm=TRUE)
  
  # return the computed values:
  return(c(PVTid, datetime, falseresp, lapses, rt_mean, rt_sd))
}


# ---- Working directory and Subject ID ----
# fill in the path to the PRESS research folder
path2dir <- '/Volumes/heronderzoek-6/Groep Geuze/25U-0078_PRESS/E_ResearchData/2_ResearchData/'

# ---- Load subject ID (not required for this script) ----
# NOTE: Subject IDs differ per instrument, even though it's the same study (privacy and data sharing reasons). 
# link different subject IDs per participant with SubjectID_koppelbestand_Castor-PVT-Garmin.csv in the PRESS research folder (/heronderzoek/Groep Geuze/25U-0078_PRESS/E_ResearchData/2_ResearchData/)
# read subject IDs voor PVT-data:
subID     <- read.csv( paste(path2dir,'SubjectID_koppelbestand_Castor-PVT-Garmin.csv', sep="") )


# ------------------------------------------------------------------------------------------------ #
# ----  Data Processing ----
# ------------------------------------------------------------------------------------------------ #

# ---- Process raw PVT data ----
# create a list of all PVT raw data files (.csv)
path2file <- paste(path2dir,'0. Ruwe data (niet in werken)/PVT/data_PVT_T0-T3_export2/csv files/', sep="")
data_list <- list.files(path=path2file, pattern=NULL, all.files=FALSE,full.names=TRUE)

# create a dataframe to save outcomes
data <- data.frame()
# loop through data files to extract outcomes
for(i in 1:length(data_list)){
  # pass the data file into the function
  filepath <- data_list[i]
  
  tryCatch({
    # save the function outcomes in 'out'
    out <- pvt_func(filepath)},
    error=function(e) -99999)
  
  tryCatch({
    # write the outcomes to dataframe
    data[i,1:6] <- out},
    error=function(e) -88888)
  
  #remove filename and outcomes before next iteration
  rm(filepath)
  rm(out)
}

# Clean the dataframe
# add column names
colnames(data) <- c("PVTid", "datetime", "falseresp", "lapses", "rt_mean", "rt_sd")
# remove empty rows
data <- data[!is.na(data$rt_mean), ]
# remove practice runs (for practice runs we entered participant IDs with more than 3 characters)
data <- data[-which(str_length(data$PVTid)>3), ]

# fill in the measurement time point based on the date
data$timepoint <- NA
data$timepoint[substr(data$datetime,1,10)=="2025-03-27"] <- "M1"
data$timepoint[substr(data$datetime,1,10)=="2025-04-10"] <- "M2"
data$timepoint[substr(data$datetime,1,10)=="2025-04-11"] <- "M2"
data$timepoint[substr(data$datetime,1,10)=="2025-04-12"] <- "M2"
data$timepoint[substr(data$datetime,1,10)=="2025-04-15"] <- "M3"
data$timepoint[substr(data$datetime,1,10)=="2025-04-24"] <- "M4"


# ------------------------------------------------------------------------------------------------ #
# ----  Data Analysis ----
# ------------------------------------------------------------------------------------------------ #

# ---- Repeated measures ANOVA M1, M3, M4 ----
# --------------------------------------------- #
# Note: Laat M2 uit analyses vanwege groot aantal missende data (data van slechts 3 participanten beschikbaar)

# create data frame for analyses with only data from timepoint M1, M3 and M4 
df <- data[which(data$timepoint=="M1"|data$timepoint=="M3"|data$timepoint=="M4"),]
# add fake PVT-ID for the dataset with missing ID
df$PVTid[is.na(df$PVTid)] <- 9999

df %>%
  get_summary_stats(rt_mean, type = "mean_sd")

# Boxplot
bxp <- ggboxplot(df, x = "timepoint", y = "rt_mean", fill="timepoint") +
  scale_colour_manual(values=c("#b8baba","#136497", "#c1d9e8")) + scale_fill_manual(values=c("#b8baba","#136497", "#c1d9e8"))
bxp 

# Histogram
ggplot(df, aes(x=rt_mean, color=timepoint)) +
  geom_histogram(fill="white")

# Shapiro
df %>%
  group_by(timepoint) %>%
  shapiro_test(rt_mean)

# Anova
res.aov <- anova_test(data = df, dv = rt_mean, wid = PVTid, within = timepoint)
get_anova_table(res.aov)

# Pairwise comparisons
pwc <- df %>%
  pairwise_t_test(rt_mean ~ timepoint, paired = FALSE, #Note: paired should be TRUE, but because of missing data the number of participants with full datasets is too small to make the pairwise comparisons
                  p.adjust.method = "bonferroni")
pwc

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "timepoint")
stats_plot <- bxp + 
  stat_pvalue_manual(pwc, step.increase = 0.07) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(title = "M1: Voor - M3: Oefening - M4: Na",
       subtitle = get_test_label(res.aov, detailed = TRUE),
       caption = get_pwc_label(pwc)) +
  xlab("Meetmoment") +
  ylab("Gemiddelde reactietijd") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Visualisation: violin plots with individual points and stats
plt <- ggbetweenstats(
  data = df,
  x = timepoint,
  y = rt_mean)

# Visualisation: simple bar plot
datp <- ddply(df, c("timepoint"), summarise,
                 N    = sum(!is.na(rt_mean)),
                 mean = mean(rt_mean, na.rm=TRUE),
                 sd   = sd(rt_mean, na.rm=TRUE),
                 se   = sd / sqrt(N) )
showbarplot <- ggplot(datp) +
  geom_bar( aes(x=timepoint, y=mean), stat="identity", fill="#136497", alpha=0.8) +
  geom_point(aes(x=timepoint, y=mean), size=3) +
  geom_errorbar( aes(x=timepoint, ymin=mean-sd, ymax=mean+sd), width=0.2, colour="black", alpha=0.9, size=0.3) +
  ggtitle("Reactietijden Alertheidstaak") +
  xlab("Meetmoment") +
  ylab("Gemiddelde reactietijd") +
  theme(panel.background = element_rect(fill="white",colour="white"), panel.border = element_blank(), panel.grid.major = element_line(colour = "grey"), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
  



# ------------------------------------------------------------------------------------------------ #
# ----  Save outcomes ----
# ------------------------------------------------------------------------------------------------ #

# Write dataframe to .csv file
write.csv(data, 
          file = paste(path2dir,'1. Verwerkte data/PVT/', 'PVT_readouts.csv',sep=""),
          row.names = FALSE)
# Re-load this data if needed:
# data <- read.csv(paste(path2dir,'1. Verwerkte data/PVT/', 'PVT_readouts.csv',sep=""))
