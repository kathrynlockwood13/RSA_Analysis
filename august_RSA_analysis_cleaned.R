# Clear the Global Environment
rm(list = ls())

# Load necessary packages
pacman::p_load(ggplot2, tidyverse, stringr, googledrive, 
               psych, ez, reshape2, Rmisc, ggsignif, car, emmeans,
               lme4, lmerTest, sjPlot)

# Set directory paths
directory <- ("/data/Kathryn/Projects/Priority/derivatives/rsa_dfs/scripts")
setwd(directory)

# Define regions
regions <- c("ffa", "hipp", "hipp_ant", "hipp_post", "loc", "mpfc", "ppa")

# Loop through each region
for (region in regions) {
  # Create a vector of subjects from 2 to 30 with leading zeros
  subjects <- sprintf("%02d", 2:30)
  
  # Create an empty list to store the data frames for this region
  region_list <- list()
  
  # Loop through subjects from sub-02 to sub-30
  for (subject_id in 2:30) {
    # Format subject number with leading zeros
    subject_formatted <- sprintf("%02d", subject_id)
    
    # Generate the path for each subject
    subject_path <- file.path("/data/Kathryn/Projects/Priority/derivatives/rsa_dfs", region, paste0("sub-", subject_formatted, "/", region, "_corr"))
    
    # Use list.files to find the specific file in each subject's folder
    subject_files <- list.files(path = subject_path,
                                pattern = paste0("*corr_matrix_", region, "_TR_Run-01_day1_by_Trial-Run-01.csv"),
                                full.names = TRUE)
    
    # Read the CSV files into data frames and add them to the region_list
    for (file_path in subject_files) {
      df_region <- read.csv(file_path)
      region_list[[length(region_list) + 1]] <- df_region
    }
  }
  # Assign the region_list to a variable named after the region
  assign(paste0(region, "_Run1day1byTrialR1_list"), region_list)
}
