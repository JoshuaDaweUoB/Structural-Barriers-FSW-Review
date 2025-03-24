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

# Optionally, assign each dataframe to the global environment
for (name in names(all_subgroup_dataframes)) {
  assign(name, all_subgroup_dataframes[[name]])
}

# Print the names of the created dataframes
message("Created dataframes: ", paste(names(all_subgroup_dataframes), collapse = ", "))

filtered_df <- fsw_data_pv_recent_ldc_bin_no
filtered_df <- create_study_effect_nums(filtered_df)

  # Create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df$effect_best_var_ln,
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)

result <- rma.mv(yi = filtered_df$effect_best_ln,  # Specify the effect size column
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num, 
                 data = filtered_df, 
                 sparse = TRUE)

# Print the results
print(result)
print(exp(coef(result)))
  