# load packages 
pacman::p_load("meta", "metafor", "readxl", "openxlsx", "tidyverse", "kableExtra", "robumeta", "clubSandwich") 

# set working directory
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence")

# settings
settings.meta(CIbracket = "(") 
settings.meta(CIseparator = "-") 

# columns
leftcols_recent <- c("study", "study_num", "effect_num", "exposure_definition_short", "exposure_time_frame", "perpetrator", "country")
leftlabs_recent <- c("Study", "Study number", "Effect number", "Exposure definition", "Exposure time frame", "Perpetrator", "Country")
leftcols_lifetime <- c("study", "study_num", "effect_num", "exposure_definition_short", "perpetrator", "country")
leftlabs_lifetime <- c("Study", "Study number", "Effect number", "Exposure definition", "Perpetrator", "Country")
rightcols <- c("effect", "ci")
rightlabs = c("Estimate", "95% CI")

# lists for loops and functions
subgroup_columns <- c("ldc_bin", "lmic_bin", "hiv_decrim", "gbv_law", "pre_2016", "recruitment", "perpetrator", "who_region", "rob_score")
dataframes <- list(fsw_data_pv_recent, fsw_data_sv_recent, fsw_data_psv_recent, fsw_data_pv_ever, fsw_data_sv_ever, fsw_data_psv_ever)  
dataframe_names <- c("fsw_data_pv_recent", "fsw_data_sv_recent", "fsw_data_psv_recent", "fsw_data_pv_ever", "fsw_data_sv_ever", "fsw_data_psv_ever")
analyses <- c("unadj", "adj", "best")
exposures <- c("recent", "ever")


#### separate variable into different dataframes for each level

# Create a list to store the dataframes for all subgroup columns
all_subgroup_dataframes <- list()

# Loop through each column in subgroup_columns
for (column in subgroup_columns) {
  # Get unique levels of the current column
  unique_levels <- unique(fsw_data_pv_recent[[column]])
  
  # Debug: Print the current column and its unique levels
  message(paste("Processing column:", column))
  print(unique_levels)
  
  # Loop through each level of the current column
  for (level in unique_levels) {
    # Filter the dataframe for the current level
    filtered_df <- fsw_data_pv_recent %>% filter(.data[[column]] == level)
    
    # Assign the filtered dataframe to the list with a meaningful name
    dataframe_name <- paste0("fsw_data_pv_recent_", column, "_", level)
    all_subgroup_dataframes[[dataframe_name]] <- filtered_df
  }
}

# assign each dataframe to the global environment
for (name in names(all_subgroup_dataframes)) {
  assign(name, all_subgroup_dataframes[[name]])
}

# Print the names of the created dataframes
message("Created dataframes: ", paste(names(all_subgroup_dataframes), collapse = ", "))

## create empty dataframe with subgroups as rows

# Initialize an empty dataframe to store subgroup results
fsw_data_pv_recent_subgroup <- data.frame(
  subgroup_level = character(),
  subgroup = character(),
  stringsAsFactors = FALSE
)

# Loop through each dataframe in the global environment to populate the subgroup dataframe
for (name in names(all_subgroup_dataframes)) {
  # Remove the prefix "fsw_data_pv_recent_" from the subgroup name
  subgroup_name <- sub("fsw_data_pv_recent_", "", name)
  
  # Extract the subgroup (e.g., "ldc_bin" from "ldc_bin_no")
  subgroup <- sub("^(.*?)_.*$", "\\1", subgroup_name)
  
  # Extract only the level (e.g., "no" from "ldc_bin_no")
  subgroup_level <- sub(".*_(.*)$", "\\1", subgroup_name)
  
  # Append the cleaned subgroup name and level to the results dataframe
  fsw_data_pv_recent_subgroup <- rbind(
    fsw_data_pv_recent_subgroup,
    data.frame(subgroup_level = subgroup_level, subgroup = subgroup, stringsAsFactors = FALSE)
  )
}

# Print the final subgroup dataframe
View(fsw_data_pv_recent_subgroup)









## run a model over each subgroup dataframe

# Loop through each dataframe in the global environment to fit models
for (name in names(all_subgroup_dataframes)) {
  # Access the dataframe directly from the global environment
  filtered_df <- get(name)
  
  # Debug: Print the name of the dataframe being processed
  message(paste("Processing dataframe:", name))
  
  # Skip if the dataframe is empty
  if (nrow(filtered_df) == 0) {
    message(paste("Skipping dataframe:", name, "- it is empty."))
    next
  }
  
  # Check for missing columns and attempt to create them
  if (!"study_num" %in% colnames(filtered_df) || !"effect_num" %in% colnames(filtered_df)) {
    message(paste("Creating study_num and effect_num for dataframe:", name))
    filtered_df <- create_study_effect_nums(filtered_df)
  }
  
  # Check if effect size columns are missing
  if (!"effect_best_var_ln" %in% colnames(filtered_df) || !"effect_best_ln" %in% colnames(filtered_df)) {
    message(paste("Skipping dataframe:", name, "- missing effect size columns."))
    next
  }
  
  # Create a covariance matrix assuming constant sampling correlation
  rho <- 0.6
  tryCatch({
    V_mat <- impute_covariance_matrix(filtered_df$effect_best_var_ln,
                                      cluster = filtered_df$study_num,
                                      r = rho,
                                      smooth_vi = TRUE)
    
    # Fit a multilevel random effects model using `rma.mv` from metafor
    result <- rma.mv(yi = filtered_df$effect_best_ln,  # Specify the effect size column
                     V = V_mat, 
                     random = ~ 1 | study_num / effect_num, 
                     data = filtered_df, 
                     sparse = TRUE)
    
    # Print the results
    print(result)
    print(exp(coef(result)))
  }, error = function(e) {
    message(paste("Error processing dataframe:", name))
    message(e)
  })
}





























# Initialize an empty dataframe to store subgroup results
fsw_data_pv_recent_subgroup <- data.frame(
  subgroup_level = character(),
  subgroup = character(),
  stringsAsFactors = FALSE
)

# Loop through each dataframe in the global environment
for (name in names(all_subgroup_dataframes)) {
  # Access the dataframe directly from the global environment
  filtered_df <- get(name)
  
  # Debug: Print the name of the dataframe being processed
  message(paste("Processing dataframe:", name))
  
  # Skip if the dataframe is empty
  if (nrow(filtered_df) == 0) {
    message(paste("Skipping dataframe:", name, "- it is empty."))
    next
  }
  
  # Check for missing columns and attempt to create them
  if (!"study_num" %in% colnames(filtered_df) || !"effect_num" %in% colnames(filtered_df)) {
    message(paste("Creating study_num and effect_num for dataframe:", name))
    filtered_df <- create_study_effect_nums(filtered_df)
  }
  
  # Check if effect size columns are missing
  if (!"effect_best_var_ln" %in% colnames(filtered_df) || !"effect_best_ln" %in% colnames(filtered_df)) {
    message(paste("Skipping dataframe:", name, "- missing effect size columns."))
    next
  }
  
  # Create a covariance matrix assuming constant sampling correlation
  rho <- 0.6
  tryCatch({
    V_mat <- impute_covariance_matrix(filtered_df$effect_best_var_ln,
                                      cluster = filtered_df$study_num,
                                      r = rho,
                                      smooth_vi = TRUE)
    
    # Fit a multilevel random effects model using `rma.mv` from metafor
    result <- rma.mv(yi = filtered_df$effect_best_ln,  # Specify the effect size column
                     V = V_mat, 
                     random = ~ 1 | study_num / effect_num, 
                     data = filtered_df, 
                     sparse = TRUE)
    
    # Print the results
    print(result)
    print(exp(coef(result)))
    
    # Remove the prefix "fsw_data_pv_recent_" from the subgroup name
    subgroup_name <- sub("fsw_data_pv_recent_", "", name)
    
    # Extract the subgroup (e.g., "who_region" from "who_region_African Region")
    subgroup <- sub("^(.*?)_.*$", "\\1", subgroup_name)
    
    # Extract only the level (e.g., "African Region" from "who_region_African Region")
    subgroup_level <- sub("^.*?_(.*)$", "\\1", subgroup_name)
    
    # Append the subgroup name and subgroup to the results dataframe
    fsw_data_pv_recent_subgroup <- rbind(
      fsw_data_pv_recent_subgroup,
      data.frame(subgroup_level = subgroup_level, subgroup = subgroup, stringsAsFactors = FALSE)
    )
  }, error = function(e) {
    message(paste("Error processing dataframe:", name))
    message(e)
  })
}

# Print the final results dataframe
View(fsw_data_pv_recent_subgroup)