# Clear the Global Environment
rm(list = ls())

# Load necessary packages
pacman::p_load(ggplot2, tidyverse, stringr, googledrive, 
               psych, ez, reshape2, Rmisc, ggsignif, car, emmeans,
               lme4, lmerTest, sjPlot, broom,ez)

# Set directory paths
directory <- ("/data/Kathryn/Projects/Priority/derivatives/rsa_dfs/scripts")
setwd(directory)

##############################################################################
#START OF LOADING IN FILES 
# Define regions and file types
regions <- c("ffa", "hipp", "hipp_ant", "hipp_post", "loc", "mpfc", "peri", "ppa", "rsp", "vtc")
file_types <- c("TR_Run-01_day1_by_Trial-Run-01.csv", "TR_Run-01_day1_by_Trial-Run-02.csv",
                "TR_Run-01_day2_by_Trial-Run-01.csv", "TR_Run-01_day2_by_Trial-Run-02.csv",
                "TR_Run-02_day2_by_Trial-Run-01.csv", "TR_Run-03_day2_by_Trial-Run-02.csv")

# Loop through each region
for (region in regions) {
  # Loop through each file type
  for (file_type in file_types) {
    # Create a vector of subjects from 2 to 30 with leading zeros
    subjects <- sprintf("%02d", 2:30)
    
    # Create an empty list to store the data frames for this region and file type
    region_file_list <- list()
    
    # Loop through subjects from sub-02 to sub-30
    for (subject_id in 2:30) {
      # Format subject number with leading zeros
      subject_formatted <- sprintf("%02d", subject_id)
      
      # Generate the path for each subject
      subject_path <- file.path("/data/Kathryn/Projects/Priority/derivatives/rsa_dfs", region, paste0("sub-", subject_formatted, "/", region, "_corr"))
      
      # Create the full file path by replacing "ppa" with the current region and file type
      full_file_path <- gsub("ppa", region, file_type)
      
      # Use list.files to find the specific file in each subject's folder
      subject_files <- list.files(path = subject_path,
                                  pattern = full_file_path,
                                  full.names = TRUE)
      
      # Read the CSV files into data frames and add them to the region_file_list
      for (file_path in subject_files) {
        df_region_file <- read.csv(file_path)
        region_file_list[[length(region_file_list) + 1]] <- df_region_file
      }
    }
    
    # Assign the region_file_list to a variable named after the region and file type
    assign(paste0(region, "_", str_remove(file_type, ".csv"), "_list"), region_file_list)
  }
}

##############################################################################
# Generate a list of all the list names
list_names <- ls(pattern = "_list")
list_names

##getting singular dataframes from the lists 
region <- "ffa"
file_type <- "TR_Run-01_day1_by_Trial-Run-01.csv"
list_name <- paste0(region, "_", str_remove(file_type, ".csv"), "_list")

data_frame_list <- get(list_name)

first_data_frame <- data_frame_list[[1]]
View(first_data_frame)
##END OF LOADING IN FILES 
################################################################################
##START OF FISHER Z TRANSFORMING FILES USING DAY 1 REST 

##getting fisher z percentiles and applying them to dfs 
# Fisher Z transformation
fisher_z_transform <- function(r) {
  0.5 * log((1 + r) / (1 - r))
}

# Initialize a data frame to store percentiles
percentile_results_df1 <- data.frame(Region = character(0), 
                                     Participant = character(0),
                                     Day = integer(0),
                                     Run = integer(0),
                                     Percentile_95 = numeric(0),
                                     Percentile_99 = numeric(0))

regions <- c("ffa", "hipp", "hipp_ant", "hipp_post", "loc", "mpfc", "peri", "ppa", "rsp", "vtc")
day <- 1
run <- 1

for (region in regions) {
  data_framesR1 <- get(paste0(region, "_TR_Run-01_day", day, "_by_Trial-Run-01_list"))
  
  for (i in seq_along(data_framesR1)) {
    df_R1 <- data_framesR1[[i]][,-1]  # Exclude the first column
    df_R1 <- apply(df_R1, 2, fisher_z_transform)  # Apply the Fisher Z transform
    
    values_list_R1 <- unlist(df_R1)
    
    percentile_95_R1 <- quantile(values_list_R1, probs = 0.95)
    percentile_99_R1 <- quantile(values_list_R1, probs = 0.99)
    
    percentile_results_df1 <- rbind(percentile_results_df1,
                                    data.frame(Region = region,
                                               Participant = sprintf("%02d", i + 1),
                                               Day = day,
                                               Run = run,
                                               Percentile_95 = percentile_95_R1,
                                               Percentile_99 = percentile_99_R1))
  }
}

View(percentile_results_df1)
###############################################################################
#thresholding 

# Fisher's z-transform function
fisher_z_transform <- function(r) {
  0.5 * log((1 + r) / (1 - r))
}

# Define a function to perform the thresholding
threshold_matrix_for_region_and_subject <- function(region, subject) {
  subject_str <- sprintf("%02d", subject)  # Convert subject number to string with leading zeros
  
  # Base path for the current region and subject
  base_path <- paste0("/data/Kathryn/Projects/Priority/derivatives/rsa_dfs/", region, "/sub-", subject_str, "/", region, "_corr")
  
  # Generate the input file name with base path
  infile <- file.path(base_path, paste0("corr_matrix_", region, "_TR_Run-02_day2_by_Trial-Run-01.csv"))
  
  # Check if the file exists
  if (!file.exists(infile)) {
    warning(paste("File", infile, "does not exist!"))
    return(NULL)
  }
  
  # Read the matrix
  matrix_df <- read.csv(infile, row.names = 1)
  
  # Apply Fisher's z-transform to the matrix
  matrix_df <- apply(matrix_df, c(1, 2), fisher_z_transform)
  
  # Find the corresponding percentile for this region and subject
  percentile <- subset(percentile_results_df1, Region == region & Participant == as.character(subject_str))$Percentile_95
  
  # If no threshold is found, warn and exit the function
  if (is.null(percentile) || length(percentile) == 0) {
    warning(paste("No threshold available for", region, "and subject", subject_str, ". Skipping..."))
    return(NULL)
  }
  
  # Threshold the matrix
  matrix_df[matrix_df < percentile] <- NA
  
  # Generate the output file name with base path
  outfile <- file.path(base_path, paste0(substr(basename(infile), 1, nchar(basename(infile))-4), "_95perc.csv"))
  
  # Save the thresholded matrix
  write.csv(matrix_df, outfile, row.names = TRUE)
}

# Given regions and subjects
regions <- c("ffa", "hipp", "hipp_ant", "hipp_post", "loc", "mpfc", "peri", "ppa", "rsp", "vtc")
subjects <- 2:30  # Numbers from 2 to 30

# Iterate through regions and subjects
for (region in regions) {
  for (subject in subjects) {
    threshold_matrix_for_region_and_subject(region, subject)
  }
}
#END OF DAY1 REST THRESHOLDING FILES 


##control shift C to comment out multiple lines in R 
###testing singular run through of thresholding 
# # Read the CSV file
# df <- read_csv("/data/Kathryn/Projects/Priority/derivatives/rsa_dfs/ffa/sub-02/ffa_corr/corr_matrix_ffa_TR_Run-02_day2_by_Trial-Run-01.csv")
# 
# 
# # Fisher Z transformation for vectors
# fisher_z_transform_vector <- function(r) {
#   0.5 * log((1 + r) / (1 - r))
# }
# 
# # Apply the function to df starting from column 2
# df[, 2:ncol(df)] <- lapply(df[, 2:ncol(df)], fisher_z_transform_vector)
# 
# 
# # Threshold the values in the dataframe
# threshold = 0.11003933
# df[df < threshold] = NA
# 
# # View the dataframe
# View(df)


###############################################
#START OF GETTING MEAN VALUE FOR ENCODING BY ENCODING MATRICES 
# Initial Data
regions <- c("ffa", "hipp", "hipp_ant", "hipp_post", "loc", "mpfc", "peri", "ppa", "rsp", "vtc")
runs <- c("Run-01", "Run-02")

# Fisher Z Transformation Function
fisher_z <- function(x) {
  if (!is.na(x)) {
    return(0.5 * log((1 + x) / (1 - x)))
  } else {
    return(NA)
  }
}

# Compute Mean excluding NA values
compute_mean_of_non_na_values <- function(df, start_col, end_col, start_row, end_row) {
  non_na_values <- c()
  for (col in start_col:end_col) {
    for (row in start_row:end_row) {
      value <- df[row, col]
      if (!is.na(value)) {
        non_na_values <- c(non_na_values, value)
      }
    }
  }
  return(mean(non_na_values, na.rm = TRUE))
}

# Results Data Frame
results_df <- data.frame(subject=character(),
                         region=character(),
                         category=character(),
                         run=character(),
                         value=numeric(),
                         stringsAsFactors=FALSE)

# Loop Through: Subjects -> Regions -> Runs
for (subject in 2:30) {
  for (region in regions) {
    for (run in runs) {
      # File Path
      file_path <- paste0("/data/Kathryn/Projects/Priority/derivatives/rsa_dfs/", region, "/sub-", sprintf("%02d", subject), "/", region, "_corr/corr_matrix_", region, "_Trial_by_Trial_", run, ".csv")
      
      # Debug: Print File Path
      cat("Processing file:", file_path, "\n")
      
      # Read CSV
      testdf <- read.csv(file_path, stringsAsFactors = FALSE)
      
      # Debug: Inspect Loaded Data
      cat("Head of dataframe after loading:\n")
      print(head(testdf))
      
      # Replace 1s with NA
      testdf[testdf == 1] <- NA
      
      # Apply Fisher Z Transformation
      for (i in 2:ncol(testdf)) {
        for (j in 2:nrow(testdf)) {
          testdf[j, i] <- fisher_z(testdf[j, i])
        }
      }
      
      # Compute Mean for Different Categories
      new_val <- compute_mean_of_non_na_values(testdf, start_col = 2, end_col = 19, start_row = 2, end_row = 19)
      old3x_val <- compute_mean_of_non_na_values(testdf, start_col = 20, end_col = 37, start_row = 20, end_row = 37)
      old1x_val <- compute_mean_of_non_na_values(testdf, start_col = 38, end_col = 55, start_row = 38, end_row = 55)
      
      # Append to Results Data Frame
      new_row <- data.frame(subject=paste0("sub-", sprintf("%02d", subject)),
                            region=region,
                            category=c("new", "old3x", "old1x"),
                            run=run,
                            value=c(new_val, old3x_val, old1x_val),
                            stringsAsFactors=FALSE)
      results_df <- rbind(results_df, new_row)
      
      # Debug: Print Updated Results DF
      cat("Current results:\n")
      print(head(results_df))
    }
  }
}

 #Final Results
 view(results_df)
 ##END OF GETTING THOSE MEAN VALUES 
 
 
###START OF ANALYZING AND VISUALIZING THE MEAN VALUES 
 # Use aggregate to calculate the mean of Run-01 and Run-02 for each subject, region, and category
 averaged_results <- aggregate(value ~ subject + region + category, data = results_df, FUN = mean)
 
 # Rename the 'value' column to 'mean_value' for clarity
 colnames(averaged_results)[which(names(averaged_results) == "value")] <- "mean_value"
 # View the averaged results
 view(averaged_results)
 
 ##hipp 
 hipp_avg_data <- averaged_results[averaged_results$region == "hipp",]
 # Run the ezANOVA
 anova_results <- ezANOVA(
   data = hipp_avg_data,               # specify the data frame
   dv = .(mean_value),                   # specify the dependent variable
   wid = .(subject),                # specify the subject identifier
   within = .(category),            # specify the within-subjects factor
   #detailed = TRUE                  # set to TRUE to get more detailed output
 )
 # View the results
 print(anova_results)
 
# Example using pairwise t-tests with correction
pairwise.t.test(hipp_avg_data$mean_value, 
                hipp_avg_data$category, 
                paired = TRUE, 
                p.adjust.method = "fdr")


regions <- c("ffa", "hipp", "hipp_ant", "hipp_post", "loc", "mpfc", "peri", "ppa", "rsp", "vtc")

# Create an empty list to store results
results_list <- list()

# Loop over each region
for(region in regions) {
  
  region_data <- averaged_results[averaged_results$region == region,]
  
  # Run the ezANOVA
  anova_results <- ezANOVA(
    data = region_data,               
    dv = .(mean_value),               
    wid = .(subject),                
    within = .(category)         
  )
  
  # Run pairwise t-tests with correction
  pairwise_results <- pairwise.t.test(region_data$mean_value, 
                                      region_data$category, 
                                      paired = TRUE, 
                                      p.adjust.method = "fdr")
  
  # Store the results in the list
  results_list[[region]] <- list(
    ANOVA = anova_results$ANOVA,
    Sphericity_Test = anova_results$`Mauchly's Test for Sphericity`,
    Sphericity_Corrections = anova_results$`Sphericity Corrections`,
    Pairwise_Tests = pairwise_results
  )
}

# View the results for a specific region, e.g., "hipp"
print(results_list[["hipp"]])
print(results_list[["hipp_ant"]])
print(results_list[["hipp_post"]])
print(results_list[["ffa"]])
print(results_list[["ppa"]])
print(results_list[["loc"]])
print(results_list[["mpfc"]])
print(results_list[["peri"]])
print(results_list[["rsp"]])
print(results_list[["vtc"]])

#EXTRACTING P VALUES FROM PAIRED T TEST
# Initialize an empty list to store the data frames
p_value_dfs <- list()

# Loop over each region in results_list
for(region in names(results_list)) {
  # Extract the p-value matrix from the Pairwise_Tests component
  p_value_matrix <- results_list[[region]]$Pairwise_Tests$p.value
  
  # Convert the p-value matrix to a data frame
  p_value_data <- as.data.frame(p_value_matrix)
  
  # Add a column for the region
  p_value_data$Region <- region
  
  # Store the data frame in the list
  p_value_dfs[[region]] <- p_value_data
}

# Combine all the data frames into a single data frame
p_value_df <- bind_rows(p_value_dfs)

# View the compiled p-value data frame
view(p_value_df)


# Function to extract ANOVA summary
extract_anova_summary <- function(anova_result) {
  data.frame(
    Effect = anova_result$Effect,
    DFn = anova_result$DFn,
    DFd = anova_result$DFd,
    F = anova_result$F,
    p = anova_result$p,
    p_less_05 = anova_result$p<.05,
    ges = anova_result$ges
  )
}

#EXTRACTING ANOVA SUMMARIES FOR EACH REGION 
# Extract ANOVA summary for each region
anova_summaries <- lapply(results_list, function(x) extract_anova_summary(x$ANOVA))

# Combine ANOVA summaries into one data frame
anova_summary_df <- do.call(rbind, anova_summaries)

# Print ANOVA summary table
print(anova_summary_df)

f_values <- sapply(results_list, function(x) x$ANOVA$F[x$ANOVA$Effect == "category"])
p_values <- sapply(results_list, function(x) x$ANOVA$p[x$ANOVA$Effect == "category"])

# Create a dataframe for plotting
anova_df <- data.frame(Region = names(f_values), F_Value = unlist(f_values), P_Value = unlist(p_values))

# Add significance column based on p-value threshold
alpha <- 0.05
anova_df$Significant <- anova_df$P_Value < alpha

view(anova_df)


# write.csv(results_df, "results_df_means.csv", row.names = FALSE)


################################################################################
# ##plotting singular portion of results df 
# # Plotting the data
# ggplot(data=results_df, aes(x=run, y=value, fill=category)) +
#   geom_bar(stat="identity", position="dodge") +
#   labs(title="Results for sub-02", x="Run", y="Value") +
#   theme_minimal()
# 
# ##ffa
# # Subset the data for 'ffa' region
# ffa_data <- averaged_results[averaged_results$region == "ffa",]
# view(ffa_data)
# # Create the plot
# ggplot(data=ffa_data, aes(x=category, y=mean_value)) +
#   geom_bar(stat="summary", fun="mean", position="dodge") +
#   labs(title="hipp", x="Category", y="Voxel Value") +
#   theme_minimal() 

# Define a function to create a plot for a given region
create_region_plot <- function(region) {
  # Subset the data for the specified region
  region_data <- averaged_results %>%
    filter(region == !!region)
  
  # Create a plot with different colored bars for each category
  p <- ggplot(data=region_data, aes(x=category, y=mean_value, fill=category)) +
    geom_bar(stat="summary", fun="mean", position="dodge") +
    labs(title=paste0(region, " Encoding Similarity Avg"), x="Category", y="Voxel Value") +
    theme_minimal() +
    scale_fill_brewer(palette="Set3")  # Use a color palette that provides distinct colors
  
  return(p)
}

# List of region names
regions <- c("ffa", "hipp", "hipp_ant", "hipp_post", "loc", "mpfc", "peri", "ppa", "rsp", "vtc")

# Use purrr::map to apply the create_region_plot function to each region
plots <- map(regions, create_region_plot)
plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]
plots[[5]]
plots[[6]]
plots[[7]]
plots[[8]]
plots[[9]]
plots[[10]]

##END OF VISUALIZATION AND ANALYZING 



####START OF TESTING TO SEE IF THERE IS A DIFFERENCE BETWEEN DAY 1 AND DAY 2 BASELINE REST 
# Initialize a data frame to store t-test results at 99% and 95% 
ttest_results_df2 <- data.frame(Region = character(0), 
                                Category = character(0), 
                                Percentile = numeric(0),
                                P_Value = numeric(0), 
                                Mean_Difference = numeric(0), 
                                Significance = character(0))

regions <- c("ffa", "hipp", "hipp_ant", "hipp_post", "loc", "mpfc", "peri", "ppa", "rsp", "vtc")
categories <- list(new = 2:19, old3x = 20:37, old1x = 38:55)
percentiles <- c(0.95, 0.99)

for (region in regions) {
  for (category_name in names(categories)) {
    for (percentile in percentiles) {
      category_range <- categories[[category_name]]
      result_df <- data.frame(Participant = sprintf("%02d", 2:30))
      
      for (day in 1:2) {
        all_values_list_R1 <- list()
        all_values_list_R2 <- list()
        
        data_framesR1 <- get(paste0(region, "_TR_Run-01_day", day, "_by_Trial-Run-01_list"), envir = .GlobalEnv)
        data_framesR2 <- get(paste0(region, "_TR_Run-01_day", day, "_by_Trial-Run-02_list"), envir = .GlobalEnv)
        
        if (day == 2) {
          data_framesR1 <- get(paste0(region, "_TR_Run-01_day", day, "_by_Trial-Run-02_list"), envir = .GlobalEnv)
          data_framesR2 <- get(paste0(region, "_TR_Run-01_day", day, "_by_Trial-Run-01_list"), envir = .GlobalEnv)
        }
        
        if(length(data_framesR1) != length(data_framesR2)) {
          stop("Mismatched lengths: data_framesR1 and data_framesR2")
        }
        
        for (i in seq_along(data_framesR1)) {
          df_R1 <- data_framesR1[[i]]
          df_R2 <- data_framesR2[[i]]
          
          values_list_R1 <- unlist(df_R1[, category_range])
          values_list_R2 <- unlist(df_R2[, category_range])
          
          all_values_list_R1[[i]] <- sort(values_list_R1)
          all_values_list_R2[[i]] <- sort(values_list_R2)
        }
        
        percentile_list_R1 <- lapply(all_values_list_R1, function(values) {
          quantile(values, probs = percentile)
        })
        
        percentile_list_R2 <- lapply(all_values_list_R2, function(values) {
          quantile(values, probs = percentile)
        })
        
        run_day_R1 <- paste0("Run1day", day, "byTrialR1")
        run_day_R2 <- paste0("Run1day", day, "byTrialR2")
        mean_day <- paste0("MeanRun1day", day)
        
        result_df[, run_day_R1] <- unlist(percentile_list_R1)
        result_df[, run_day_R2] <- unlist(percentile_list_R2)
        result_df[, mean_day] <- rowMeans(data.frame(R1 = unlist(percentile_list_R1), R2 = unlist(percentile_list_R2)))
      }
      
      ttest_result <- t.test(result_df$MeanRun1day1, result_df$MeanRun1day2, paired = TRUE)
      
      significance <- ifelse(ttest_result$p.value < 0.001, "***",
                             ifelse(ttest_result$p.value < 0.01, "**",
                                    ifelse(ttest_result$p.value < 0.05, "*", "NS")))
      
      ttest_results_df2 <- rbind(ttest_results_df2, 
                                 data.frame(Region = region, 
                                            Category = category_name, 
                                            Percentile = percentile,
                                            P_Value = ttest_result$p.value, 
                                            Mean_Difference = ttest_result$estimate, 
                                            Significance = significance))
    }
  }
}

# View the results
print(ttest_results_df2)




# Initialize a data frame to store percentiles
percentile_results_df <- data.frame(Region = character(0), 
                                    Participant = character(0),
                                    Day = integer(0),
                                    Run = integer(0),
                                    Percentile_95 = numeric(0),
                                    Percentile_99 = numeric(0))
####END OF TESTING BASELINE 1 AND 2 

########THRESHOLD DF ANALYSIS START##########################
# Define the path patterns for both encoding runs with a placeholder for region
path_pattern_run01 <- "/data/Kathryn/Projects/Priority/derivatives/rsa_dfs/%s/sub-%s/ffa_corr/corr_matrix_ffa_TR_Run-02_day2_by_Trial-Run-01_95perc.csv"
path_pattern_run02 <- "/data/Kathryn/Projects/Priority/derivatives/rsa_dfs/%s/sub-%s/ffa_corr/corr_matrix_ffa_TR_Run-03_day2_by_Trial-Run-02_95perc.csv"

# Define the regions to loop through
regions <- c("ffa", "hipp", "hipp_ant", "hipp_post", "loc", "mpfc", "peri", "ppa", "rsp", "vtc")

# Initialize an empty dataframe to store results with additional columns for encoding 1, encoding 2, and region
results_df_TR_groupmax <- tibble(subject = character(),
                                 TR = character(),
                                 col_num = character(),
                                 encoding_1 = numeric(),
                                 encoding_2 = numeric(),
                                 region = character())

# Function to safely get the maximum value
safe_max <- function(values) {
  max_val <- max(values, na.rm = TRUE)
  if (is.infinite(max_val)) return(NA)
  return(max_val)
}

# Loop through each subject and each region
for (region in regions) {
  for (subject_num in 2:30) {
    subject_str <- sprintf("%02d", subject_num)
    file_path_run01 <- sprintf(path_pattern_run01, region, subject_str)
    file_path_run02 <- sprintf(path_pattern_run02, region, subject_str)
    
    # Check if both files exist
    if (file.exists(file_path_run01) & file.exists(file_path_run02)) {
      df_run01 <- read.csv(file_path_run01)
      df_run02 <- read.csv(file_path_run02)
      
      # Loop through each row
      for (i in 1:nrow(df_run01)) {
        TR <- as.character(df_run01[i, 1])
        
        # Extract maximum value for each column group in both encoding runs
        max_val_group1_run01 <- safe_max(df_run01[i, 2:19])
        max_val_group2_run01 <- safe_max(df_run01[i, 20:37])
        max_val_group3_run01 <- safe_max(df_run01[i, 38:55])
        
        max_val_group1_run02 <- safe_max(df_run02[i, 2:19])
        max_val_group2_run02 <- safe_max(df_run02[i, 20:37])
        max_val_group3_run02 <- safe_max(df_run02[i, 38:55])
        
        # Append results to results_df with separate columns for each encoding run and the current region
        results_df_TR_groupmax <- rbind(results_df_TR_groupmax, 
                                        tibble(subject = subject_str,
                                               TR = TR,
                                               col_num = c("new", "old3x", "old1x"),
                                               encoding_1 = c(max_val_group1_run01, max_val_group2_run01, max_val_group3_run01),
                                               encoding_2 = c(max_val_group1_run02, max_val_group2_run02, max_val_group3_run02),
                                               region = region))
      }
    }
  }
}

# View the results
view(results_df_TR_groupmax)


# Add binary indicator columns
results_df_TR_groupmax <- results_df_TR_groupmax %>%
  mutate(encoding_1_bi = ifelse(is.na(encoding_1), 0, 1),
         encoding_2_bi = ifelse(is.na(encoding_2), 0, 1))

# Add the encoding_both column
results_df_TR_groupmax <- results_df_TR_groupmax %>%
  mutate(encoding_both = ifelse(encoding_1_bi == 1 & encoding_2_bi == 1, 1, 0))

# View the updated dataframe
view(results_df_TR_groupmax)

# Create an empty dataframe to store the counts for each subject and category
final_counts_df <- data.frame(subject = character(),
                              new = integer(),
                              old3x = integer(),
                              old1x = integer(),
                              stringsAsFactors = FALSE)

# Unique subjects and categories
unique_subjects <- unique(results_df_TR_groupmax$subject)
unique_col_nums <- unique(results_df_TR_groupmax$col_num)

# Loop over each subject
for (subj in unique_subjects) {
  # Initialize counts for each category
  counts <- setNames(numeric(length(unique_col_nums)), unique_col_nums)
  
  # Loop over each category
  for (col_num in unique_col_nums) {
    # Calculate the sum of encoding_both == 1 for each subject-category combination
    counts[col_num] <- sum(results_df_TR_groupmax$encoding_both[results_df_TR_groupmax$subject == subj & results_df_TR_groupmax$col_num == col_num])
  }
  
  # Combine the subject and its counts into a dataframe
  subj_df <- data.frame(subject = subj, t(counts), stringsAsFactors = FALSE)
  
  # Bind this subject's counts to the final dataframe
  final_counts_df <- rbind(final_counts_df, subj_df)
}

# Replace NA with 0 if there are any NAs in the final dataframe
final_counts_df[is.na(final_counts_df)] <- 0

# Print the final counts dataframe
view(final_counts_df)

final_counts_long <- pivot_longer(final_counts_df, 
                                  cols = c(new, old3x, old1x), 
                                  names_to = "category", 
                                  values_to = "count")

ggplot(final_counts_long, aes(x = count, y = subject, color = category)) +
  geom_point(position = position_dodge(width = 0.5)) +
  theme_minimal() +
  labs(x = "Count", y = "Subject")



# Determine the category with the highest TR value for each subject and TR
highest_category <- results_df_TR_groupmax %>%
  group_by(subject, TR) %>%
  filter(TR_value == max(TR_value, na.rm = TRUE)) %>%
  select(subject, TR, col_num, TR_value) %>%
  distinct()

view(highest_category)


# Corrected path pattern
path_pattern <- "/data/Kathryn/Projects/Priority/derivatives/rsa_dfs/ffa/sub-%s/ffa_corr/corr_matrix_ffa_TR_%s_day2_by_%s_95perc.csv"

# Function to safely get the maximum value
safe_max <- function(values) {
  max_val <- max(values, na.rm = TRUE)
  if (is.infinite(max_val)) return(NA)
  return(max_val)
}

# Temporary storage for values from both runs
temp_df <- tibble(subject = character(),
                  TR = character(),
                  col_num = character(),
                  TR_value = numeric(),
                  Run = character())

# Loop through each subject
for (subject_num in 2:30) {
  subject_str <- sprintf("%02d", subject_num)
  
  # Define mapping of runs to trial runs
  run_trial_mapping <- list("Run-02" = "Trial-Run-01", "Run-03" = "Trial-Run-02")
  
  # Loop through each run
  for (run_str in names(run_trial_mapping)) {
    trial_run_str <- run_trial_mapping[[run_str]]
    file_path <- sprintf(path_pattern, subject_str, run_str, trial_run_str)
    
    # Check if the file exists
    if (file.exists(file_path)) {
      df <- read.csv(file_path)
      
      # Loop through each row
      for (i in 1:nrow(df)) {
        TR <- as.character(df[i, 1])
        
        # Extract maximum value for each column group
        max_val_group1 <- safe_max(df[i, 2:19])
        max_val_group2 <- safe_max(df[i, 20:37])
        max_val_group3 <- safe_max(df[i, 38:55])
        
        # Append results to temp_df
        temp_df <- rbind(temp_df, 
                         tibble(subject = subject_str,
                                TR = TR,
                                col_num = c("new", "old3x", "old1x"),
                                TR_value = c(max_val_group1, max_val_group2, max_val_group3),
                                Run = run_str))
      }
    }
    
  }
}

# Retain only the highest TR value across both runs for each TR and category
df_encoding_avg_highest_TR <- temp_df %>%
  group_by(subject, TR, col_num) %>%
  summarize(TR_value = max(TR_value, na.rm = TRUE),
            Run = ifelse(TR_value == max(temp_df$TR_value[temp_df$Run == "Run-02" & 
                                                            temp_df$subject == first(subject) & 
                                                            temp_df$TR == first(TR) & 
                                                            temp_df$col_num == first(col_num)], na.rm = TRUE),
                         "Run-02", "Run-03")) %>%
  ungroup()

# View the results
view(df_encoding_avg_highest_TR)


##HIPPP 

# Define the path patterns for both encoding runs
path_pattern_run01 <- "/data/Kathryn/Projects/Priority/derivatives/rsa_dfs/hipp/sub-%s/hipp_corr/corr_matrix_hipp_TR_Run-02_day2_by_Trial-Run-01_95perc.csv"
path_pattern_run02 <- "/data/Kathryn/Projects/Priority/derivatives/rsa_dfs/hipp/sub-%s/hipp_corr/corr_matrix_hipp_TR_Run-03_day2_by_Trial-Run-02_95perc.csv"

# Initialize an empty dataframe to store results with additional columns for encoding 1 and encoding 2
results_df_TR_groupmax_hipp <- tibble(subject = character(),
                                 TR = character(),
                                 col_num = character(),
                                 encoding_1 = numeric(),
                                 encoding_2 = numeric())

# Function to safely get the maximum value
safe_max <- function(values) {
  max_val <- max(values, na.rm = TRUE)
  if (is.infinite(max_val)) return(NA)
  return(max_val)
}

# Loop through each subject
for (subject_num in 2:30) {
  subject_str <- sprintf("%02d", subject_num)
  file_path_run01 <- sprintf(path_pattern_run01, subject_str)
  file_path_run02 <- sprintf(path_pattern_run02, subject_str)
  
  # Check if both files exist
  if (file.exists(file_path_run01) & file.exists(file_path_run02)) {
    df_run01 <- read.csv(file_path_run01)
    df_run02 <- read.csv(file_path_run02)
    
    # Loop through each row
    for (i in 1:nrow(df_run01)) {
      TR <- as.character(df_run01[i, 1])
      
      # Extract maximum value for each column group in both encoding runs
      max_val_group1_run01 <- safe_max(df_run01[i, 2:19])
      max_val_group2_run01 <- safe_max(df_run01[i, 20:37])
      max_val_group3_run01 <- safe_max(df_run01[i, 38:55])
      
      max_val_group1_run02 <- safe_max(df_run02[i, 2:19])
      max_val_group2_run02 <- safe_max(df_run02[i, 20:37])
      max_val_group3_run02 <- safe_max(df_run02[i, 38:55])
      
      # Append results to results_df with separate columns for each encoding run
      results_df_TR_groupmax_hipp <- rbind(results_df_TR_groupmax_hipp, 
                                      tibble(subject = subject_str,
                                             TR = TR,
                                             col_num = c("new", "old3x", "old1x"),
                                             encoding_1 = c(max_val_group1_run01, max_val_group2_run01, max_val_group3_run01),
                                             encoding_2 = c(max_val_group1_run02, max_val_group2_run02, max_val_group3_run02)))
    }
  }
}
 
# View the results
view(results_df_TR_groupmax_hipp)



# Add binary indicator columns
results_df_TR_groupmax_hipp <- results_df_TR_groupmax_hipp %>%
  mutate(encoding_1_bi = ifelse(is.na(encoding_1), 0, 1),
         encoding_2_bi = ifelse(is.na(encoding_2), 0, 1))

# Add the encoding_both column
results_df_TR_groupmax_hipp <- results_df_TR_groupmax_hipp %>%
  mutate(encoding_both = ifelse(encoding_1_bi == 1 & encoding_2_bi == 1, 1, 0))


view(results_df_TR_groupmax_hipp)
# Create an empty dataframe to store the counts for each subject and category
final_counts_df_hipp <- data.frame(subject = character(),
                              new = integer(),
                              old3x = integer(),
                              old1x = integer(),
                              stringsAsFactors = FALSE)

# Unique subjects and categories
unique_subjects <- unique(results_df_TR_groupmax_hipp$subject)
unique_col_nums <- unique(results_df_TR_groupmax_hipp$col_num)
# Loop over each subject
for (subj in unique_subjects) {
  # Initialize counts for each category
  counts <- setNames(numeric(length(unique_col_nums)), unique_col_nums)
  
  # Loop over each category
  for (col_num in unique_col_nums) {
    # Calculate the sum of encoding_both == 1 for each subject-category combination
    counts[col_num] <- sum(results_df_TR_groupmax_hipp$encoding_both[results_df_TR_groupmax_hipp$subject == subj & results_df_TR_groupmax_hipp$col_num == col_num])
  }
  
  # Combine the subject and its counts into a dataframe
  subj_df <- data.frame(subject = subj, t(counts), stringsAsFactors = FALSE)
  
  # Bind this subject's counts to the final dataframe
  final_counts_df_hipp <- rbind(final_counts_df_hipp, subj_df)
}

# Replace NA with 0 if there are any NAs in the final dataframe
final_counts_df_hipp[is.na(final_counts_df_hipp)] <- 0

# Print the final counts dataframe
view(final_counts_df_hipp)
final_counts_long_hipp <- pivot_longer(final_counts_df_hipp, 
                                  cols = c(new, old3x, old1x), 
                                  names_to = "category", 
                                  values_to = "count")

ggplot(final_counts_long_hipp, aes(x = subject, y = count, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Subject", y = "Count", fill = "Category") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) # Rotate x labels for readability


# Assuming final_counts_long is already created and has 'subject', 'category', and 'count' columns
ggplot(final_counts_long_hipp, aes(x = subject, y = count, fill = category)) +
  geom_bar(stat = "identity") + # Use dodge to place bars side by side
  scale_fill_brewer(palette = "Set1") + # Optional: Use a color palette for the fill
  theme_minimal() +
  labs(x = "Subject", y = "Count", fill = "Category") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x labels for readability


# Plot the total count for each category
ggplot(category_totals, aes(x = category, y = total_count, fill = category)) +
  geom_bar(stat = "identity") + # Use bars to represent the total count
  theme_minimal() +
  labs(x = "Category", y = "Total Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 