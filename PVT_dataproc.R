# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
#
# Description: R-script to read and process PVT data
#             - Study name: PRESS (dossiernr.: 25U-0078)
#             - Measure/instrument: smartphone-based Psychommotor Vigilance Task (PVT)
#             - Data type: behavioral performance (reaction times and accuracy)
#             - Design: within-subjects longitudinal, 4 timepoints (M1: 27th 03/'25, M2: 10th-11th 04/'25, M3: 15th 04/'25,, M4: 24th 04/'25)
#
# Date:         May 2025
# R.version:    4.4.0 (2024-04-24) -- "Puppy Cup"
# Rstudio:      2025.05.0+496
#
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

# load required packages
require("reshape2")
require("ggplot2")
require("tidyverse")
require("dplyr")

# clear environment
rm(list=ls())

# fill in the path to the PRESS research folder
path2dir <- '/Volumes/heronderzoek-6/Groep Geuze/25U-0078_PRESS/E_ResearchData/2_ResearchData/'

# Load subject ID
# NOTE: Subject IDs differ per instrument, even though it's the same study (privacy and data sharing reasons). 
# link different subject IDs per participant with SubjectID_koppelbestand_Castor-PVT-Garmin.csv in the PRESS research folder (/heronderzoek/Groep Geuze/25U-0078_PRESS/E_ResearchData/2_ResearchData/)
# read subject IDs voor PVT-data:
subID     <- read.csv( paste(path2dir,'SubjectID_koppelbestand_Castor-PVT-Garmin.csv', sep="") )


# ---- Process raw data ----
# create a list of all PVT raw data files (.csv)
path2file <- paste(path2dir,'0. Ruwe data (niet in werken)/PVT/data_PVT_T0-T3_export2/csv files/', sep="")
data_list <- list.files(path=path2file, pattern=NULL, all.files=FALSE,full.names=TRUE)

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


# ---- Clean the dataframe ----
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


# ---- Save data frame ----
write.csv(data, 
          file = paste(path2dir,'1. Verwerkte data/PVT/', 'PVT_readouts.csv',sep=""),
          row.names = FALSE)

