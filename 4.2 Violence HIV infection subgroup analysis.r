# load packages 
pacman::p_load("meta", "metafor", "readxl", "openxlsx", "tidyverse", "kableExtra", "robumeta", "clubSandwich") 

# set working directory
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence")

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

## create empty dataframe with subgroup strata as rows

# Create the fsw_data_pv_recent_subgroup dataframe
fsw_data_pv_recent_subgroup <- data.frame(
  subgroup_level = character(),
  subgroup = character(),
  stringsAsFactors = FALSE
)

# Populate the fsw_data_pv_recent_subgroup dataframe
for (name in names(all_subgroup_dataframes)) {
  # Remove the prefix "fsw_data_pv_recent_" from the subgroup name
  subgroup_name <- sub("fsw_data_pv_recent_", "", name)
  
  # Extract the subgroup (e.g., "ldc_bin" from "ldc_bin_no")
  subgroup <- sub("^(.*?)_.*$", "\\1", subgroup_name)
  
  # Extract only the level (e.g., "no" from "ldc_bin_no")
  subgroup_level <- sub("^[^_]*_", "", subgroup_name)
  
  # Append the cleaned subgroup name and level to the results dataframe
  fsw_data_pv_recent_subgroup <- rbind(
    fsw_data_pv_recent_subgroup,
    data.frame(subgroup_level = subgroup_level, subgroup = subgroup, stringsAsFactors = FALSE)
  )
}

# Add columns to store model results
fsw_data_pv_recent_subgroup$model_coef <- NA
fsw_data_pv_recent_subgroup$model_ci_lower <- NA
fsw_data_pv_recent_subgroup$model_ci_upper <- NA
fsw_data_pv_recent_subgroup$studies<- NA
fsw_data_pv_recent_subgroup$estimates <- NA

## run a model over each subgroup dataframe that conducts the model for each strata of subgroups
## and stores the results in the fsw_data_pv_recent_subgroup dataframe

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
    filtered_df <- filtered_df %>% filter(use == "yes")
  }
  
  # Check if effect size columns are missing
  if (!"effect_best_var_ln" %in% colnames(filtered_df) || !"effect_best_ln" %in% colnames(filtered_df)) {
    message(paste("Skipping dataframe:", name, "- missing effect size columns."))
    next
  }
  
  # Create a covariance matrix assuming constant sampling correlation
  rho <- 0.6
  tryCatch({
    if (nrow(filtered_df) > 1) {
      # Fit a multilevel random effects model using `rma.mv` from metafor
      V_mat <- impute_covariance_matrix(filtered_df$effect_best_var_ln,
                                        cluster = filtered_df$study_num,
                                        r = rho,
                                        smooth_vi = TRUE)
      
      result <- rma.mv(yi = filtered_df$effect_best_ln,  # Specify the effect size column
                       V = V_mat, 
                       random = ~ 1 | study_num / effect_num, 
                       data = filtered_df, 
                       sparse = TRUE)
      
      # Extract model results
      coef <- exp(coef(result))  # Exponentiated coefficient
      ci_lower <- exp(result$ci.lb)  # Lower bound of confidence interval
      ci_upper <- exp(result$ci.ub)  # Upper bound of confidence interval
    } else {
      # Handle cases where k == 1
      coef <- exp(filtered_df$effect_best_ln)
      ci_lower <- exp(filtered_df$effect_best_ln - 1.96 * sqrt(filtered_df$effect_best_var_ln))
      ci_upper <- exp(filtered_df$effect_best_ln + 1.96 * sqrt(filtered_df$effect_best_var_ln))
    }
    
    # Extract subgroup and level
    subgroup_name <- sub("fsw_data_pv_recent_", "", name)
    subgroup <- sub("^(.*?)_.*$", "\\1", subgroup_name)
    subgroup_level <- sub("^[^_]*_", "", subgroup_name)
    
    # Debug: Print extracted values
    message(paste("Extracted subgroup:", subgroup, "and level:", subgroup_level))
    
    # Find the corresponding row in fsw_data_pv_recent_subgroup
    row_index <- which(fsw_data_pv_recent_subgroup$subgroup_level == subgroup_level & 
                       fsw_data_pv_recent_subgroup$subgroup == subgroup)
    
    # Debug: Check if row_index is found
    if (length(row_index) == 0) {
      message(paste("No matching row found for subgroup:", subgroup, "and level:", subgroup_level))
      next
    }
    
    # Count the number of unique studies and estimates
    num_studies <- length(unique(filtered_df$study_num))
    num_estimates <- nrow(filtered_df)
    
    # Store the results in the dataframe
    fsw_data_pv_recent_subgroup$model_coef[row_index] <- coef
    fsw_data_pv_recent_subgroup$model_ci_lower[row_index] <- ci_lower
    fsw_data_pv_recent_subgroup$model_ci_upper[row_index] <- ci_upper
    fsw_data_pv_recent_subgroup$studies[row_index] <- num_studies
    fsw_data_pv_recent_subgroup$estimates[row_index] <- num_estimates
    
    # Print the results
    if (nrow(filtered_df) > 1) {
      print(result)
    } else {
      message(paste("Single study result for subgroup:", subgroup, "and level:", subgroup_level))
      message(paste("Effect size:", coef, "CI:", ci_lower, "-", ci_upper))
    }
  }, error = function(e) {
    message(paste("Error processing dataframe:", name))
    message(e)
  })
}

# Drop rows with missing values from fsw_data_pv_recent_subgroup
fsw_data_pv_recent_subgroup <- fsw_data_pv_recent_subgroup %>%
  drop_na()

# create log of OR and 95% CIs
fsw_data_pv_recent_subgroup <- fsw_data_pv_recent_subgroup %>%
  mutate(
    log_model_coef = log(model_coef),
    log_model_ci_lower = log(model_ci_lower),
    log_model_ci_upper = log(model_ci_upper)
  )

## renaming rows

# Rename values in the subgroup_level column
fsw_data_pv_recent_subgroup <- fsw_data_pv_recent_subgroup %>%
  mutate(subgroup_level = case_when(
    subgroup_level == "bin_no" ~ "No",
    subgroup_level == "bin_yes" ~ "Yes",
    subgroup_level == "decrim_partial" ~ "Partial",
    subgroup_level == "decrim_yes" ~ "Yes",
    subgroup_level == "decrim_no" ~ "No",
    subgroup_level == "law_yes" ~ "Yes",
    subgroup_level == "law_no" ~ "No",
    subgroup_level == "2016_FALSE" ~ "No",
    subgroup_level == "2016_TRUE" ~ "Yes",
    subgroup_level == "region_African Region" ~ "African Region",
    subgroup_level == "region_Region of the Americas" ~ "Region of the Americas",
    subgroup_level == "region_South-East Asia Region" ~ "South-East Asia Region",
    subgroup_level == "region_European Region" ~ "European Region",
    subgroup_level == "region_Eastern Mediterranean Region" ~ "Eastern Mediterranean Region",
    subgroup_level == "region_Western Pacific Region" ~ "Western Pacific Region",
    subgroup_level == "score_Good" ~ "Good",
    subgroup_level == "score_Very good" ~ "Very good",
    subgroup_level == "score_Satisfactory" ~ "Satisfactory",
    subgroup_level == "score_Unsatisfactory" ~ "Unsatisfactory",
    TRUE ~ subgroup_level  # Keep other values unchanged
  ))

fsw_data_pv_recent_subgroup <- fsw_data_pv_recent_subgroup %>%
  mutate(subgroup = case_when(
    subgroup == "ldc" ~ "Least developed country",
    subgroup == "lmic" ~ "Lower-middle income country",
    subgroup == "hiv" ~ "HIV decriminalisation",
    subgroup == "gbv" ~ "Gender-based violence law",
    subgroup == "pre" ~ "Published before 2016",
    subgroup == "recruitment" ~ "Recruitment",
    subgroup == "perpetrator" ~ "Perpetrator",
    subgroup == "who" ~ "WHO region",
    subgroup == "rob" ~ "Risk of bias score",
    TRUE ~ subgroup  # Keep other values unchanged
  ))

# View the updated dataframe
View(fsw_data_pv_recent_subgroup)

# Perform meta-analysis
fsw_data_pv_recent_subgroup_forest <- metagen(
  TE = log_model_coef,
  lower = log_model_ci_lower,
  upper = log_model_ci_upper,
  data = fsw_data_pv_recent_subgroup,
  sm = "OR",
  method.tau = "DL",
  comb.fixed = FALSE,
  comb.random = FALSE, 
  backtransf = TRUE,
  byvar = subgroup, 
  text.random = "Overall"
)

# Print summary
print(summary(fsw_data_pv_recent_subgroup_forest))

# Save forest plot
png(filename = "Plots/overall plots/recent_pv_subgroup.png", width = 30, height = 30, units = "cm", res = 600)
forest(
  fsw_data_pv_recent_subgroup_forest,
  sortvar = subgroup,
  xlim = c(0.2, 4),             
  leftcols = c("subgroup_level", "studies", "estimates"), 
  leftlabs = c("Category", "Nb studies", "Nb estimates"),
  pooled.totals = FALSE,
  xintercept = 1,
  addrow.overall = TRUE,
  test.subgroup = FALSE,
  overall.hetstat = FALSE,
  overall = FALSE,
  labeltext = TRUE,
  col.subgroup = "black",
  print.subgroup.name = FALSE
)
dev.off(


## ever exposed to physical violence
#### separate variable into different dataframes for each level

# Create a list to store the dataframes for all subgroup columns
all_subgroup_dataframes <- list()

# Loop through each column in subgroup_columns
for (column in subgroup_columns) {
  # Get unique levels of the current column
  unique_levels <- unique(fsw_data_pv_ever[[column]])
  
  # Debug: Print the current column and its unique levels
  message(paste("Processing column:", column))
  print(unique_levels)
  
  # Loop through each level of the current column
  for (level in unique_levels) {
    # Filter the dataframe for the current level
    filtered_df <- fsw_data_pv_ever %>% filter(.data[[column]] == level)
    
    # Assign the filtered dataframe to the list with a meaningful name
    dataframe_name <- paste0("fsw_data_pv_ever_", column, "_", level)
    all_subgroup_dataframes[[dataframe_name]] <- filtered_df
  }
}

# assign each dataframe to the global environment
for (name in names(all_subgroup_dataframes)) {
  assign(name, all_subgroup_dataframes[[name]])
}

# Print the names of the created dataframes
message("Created dataframes: ", paste(names(all_subgroup_dataframes), collapse = ", "))

## create empty dataframe with subgroup strata as rows

# Create the fsw_data_pv_ever_subgroup dataframe
fsw_data_pv_ever_subgroup <- data.frame(
  subgroup_level = character(),
  subgroup = character(),
  stringsAsFactors = FALSE
)

# Populate the fsw_data_pv_ever_subgroup dataframe
for (name in names(all_subgroup_dataframes)) {
  # Remove the prefix "fsw_data_pv_ever_" from the subgroup name
  subgroup_name <- sub("fsw_data_pv_ever_", "", name)
  
  # Extract the subgroup (e.g., "ldc_bin" from "ldc_bin_no")
  subgroup <- sub("^(.*?)_.*$", "\\1", subgroup_name)
  
  # Extract only the level (e.g., "no" from "ldc_bin_no")
  subgroup_level <- sub("^[^_]*_", "", subgroup_name)
  
  # Append the cleaned subgroup name and level to the results dataframe
  fsw_data_pv_ever_subgroup <- rbind(
    fsw_data_pv_ever_subgroup,
    data.frame(subgroup_level = subgroup_level, subgroup = subgroup, stringsAsFactors = FALSE)
  )
}

# Add columns to store model results
fsw_data_pv_ever_subgroup$model_coef <- NA
fsw_data_pv_ever_subgroup$model_ci_lower <- NA
fsw_data_pv_ever_subgroup$model_ci_upper <- NA
fsw_data_pv_ever_subgroup$studies<- NA
fsw_data_pv_ever_subgroup$estimates <- NA

## run a model over each subgroup dataframe that conducts the model for each strata of subgroups
## and stores the results in the fsw_data_pv_ever_subgroup dataframe

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
    filtered_df <- filtered_df %>% filter(use == "yes")
  }
  
  # Check if effect size columns are missing
  if (!"effect_best_var_ln" %in% colnames(filtered_df) || !"effect_best_ln" %in% colnames(filtered_df)) {
    message(paste("Skipping dataframe:", name, "- missing effect size columns."))
    next
  }
  
  # Create a covariance matrix assuming constant sampling correlation
  rho <- 0.6
  tryCatch({
    if (nrow(filtered_df) > 1) {
      # Fit a multilevel random effects model using `rma.mv` from metafor
      V_mat <- impute_covariance_matrix(filtered_df$effect_best_var_ln,
                                        cluster = filtered_df$study_num,
                                        r = rho,
                                        smooth_vi = TRUE)
      
      result <- rma.mv(yi = filtered_df$effect_best_ln,  # Specify the effect size column
                       V = V_mat, 
                       random = ~ 1 | study_num / effect_num, 
                       data = filtered_df, 
                       sparse = TRUE)
      
      # Extract model results
      coef <- exp(coef(result))  # Exponentiated coefficient
      ci_lower <- exp(result$ci.lb)  # Lower bound of confidence interval
      ci_upper <- exp(result$ci.ub)  # Upper bound of confidence interval
    } else {
      # Handle cases where k == 1
      coef <- exp(filtered_df$effect_best_ln)
      ci_lower <- exp(filtered_df$effect_best_ln - 1.96 * sqrt(filtered_df$effect_best_var_ln))
      ci_upper <- exp(filtered_df$effect_best_ln + 1.96 * sqrt(filtered_df$effect_best_var_ln))
    }
    
    # Extract subgroup and level
    subgroup_name <- sub("fsw_data_pv_ever_", "", name)
    subgroup <- sub("^(.*?)_.*$", "\\1", subgroup_name)
    subgroup_level <- sub("^[^_]*_", "", subgroup_name)
    
    # Debug: Print extracted values
    message(paste("Extracted subgroup:", subgroup, "and level:", subgroup_level))
    
    # Find the corresponding row in fsw_data_pv_ever_subgroup
    row_index <- which(fsw_data_pv_ever_subgroup$subgroup_level == subgroup_level & 
                       fsw_data_pv_ever_subgroup$subgroup == subgroup)
    
    # Debug: Check if row_index is found
    if (length(row_index) == 0) {
      message(paste("No matching row found for subgroup:", subgroup, "and level:", subgroup_level))
      next
    }
    
    # Count the number of unique studies and estimates
    num_studies <- length(unique(filtered_df$study_num))
    num_estimates <- nrow(filtered_df)
    
    # Store the results in the dataframe
    fsw_data_pv_ever_subgroup$model_coef[row_index] <- coef
    fsw_data_pv_ever_subgroup$model_ci_lower[row_index] <- ci_lower
    fsw_data_pv_ever_subgroup$model_ci_upper[row_index] <- ci_upper
    fsw_data_pv_ever_subgroup$studies[row_index] <- num_studies
    fsw_data_pv_ever_subgroup$estimates[row_index] <- num_estimates
    
    # Print the results
    if (nrow(filtered_df) > 1) {
      print(result)
    } else {
      message(paste("Single study result for subgroup:", subgroup, "and level:", subgroup_level))
      message(paste("Effect size:", coef, "CI:", ci_lower, "-", ci_upper))
    }
  }, error = function(e) {
    message(paste("Error processing dataframe:", name))
    message(e)
  })
}

# Drop rows with missing values from fsw_data_pv_ever_subgroup
fsw_data_pv_ever_subgroup <- fsw_data_pv_ever_subgroup %>%
  drop_na()

# create log of OR and 95% CIs
fsw_data_pv_ever_subgroup <- fsw_data_pv_ever_subgroup %>%
  mutate(
    log_model_coef = log(model_coef),
    log_model_ci_lower = log(model_ci_lower),
    log_model_ci_upper = log(model_ci_upper)
  )

## renaming rows

# Rename values in the subgroup_level column
fsw_data_pv_ever_subgroup <- fsw_data_pv_ever_subgroup %>%
  mutate(subgroup_level = case_when(
    subgroup_level == "bin_no" ~ "No",
    subgroup_level == "bin_yes" ~ "Yes",
    subgroup_level == "decrim_partial" ~ "Partial",
    subgroup_level == "decrim_yes" ~ "Yes",
    subgroup_level == "decrim_no" ~ "No",
    subgroup_level == "law_yes" ~ "Yes",
    subgroup_level == "law_no" ~ "No",
    subgroup_level == "2016_FALSE" ~ "No",
    subgroup_level == "2016_TRUE" ~ "Yes",
    subgroup_level == "region_African Region" ~ "African Region",
    subgroup_level == "region_Region of the Americas" ~ "Region of the Americas",
    subgroup_level == "region_South-East Asia Region" ~ "South-East Asia Region",
    subgroup_level == "region_European Region" ~ "European Region",
    subgroup_level == "region_Eastern Mediterranean Region" ~ "Eastern Mediterranean Region",
    subgroup_level == "region_Western Pacific Region" ~ "Western Pacific Region",
    subgroup_level == "score_Good" ~ "Good",
    subgroup_level == "score_Very good" ~ "Very good",
    subgroup_level == "score_Satisfactory" ~ "Satisfactory",
    subgroup_level == "score_Unsatisfactory" ~ "Unsatisfactory",
    TRUE ~ subgroup_level  # Keep other values unchanged
  ))

fsw_data_pv_ever_subgroup <- fsw_data_pv_ever_subgroup %>%
  mutate(subgroup = case_when(
    subgroup == "ldc" ~ "Least developed country",
    subgroup == "lmic" ~ "Lower-middle income country",
    subgroup == "hiv" ~ "HIV decriminalisation",
    subgroup == "gbv" ~ "Gender-based violence law",
    subgroup == "pre" ~ "Published before 2016",
    subgroup == "recruitment" ~ "Recruitment",
    subgroup == "perpetrator" ~ "Perpetrator",
    subgroup == "who" ~ "WHO region",
    subgroup == "rob" ~ "Risk of bias score",
    TRUE ~ subgroup  # Keep other values unchanged
  ))

# View the updated dataframe
View(fsw_data_pv_ever_subgroup)

# Perform meta-analysis
fsw_data_pv_ever_subgroup_forest <- metagen(
  TE = log_model_coef,
  lower = log_model_ci_lower,
  upper = log_model_ci_upper,
  data = fsw_data_pv_ever_subgroup,
  sm = "OR",
  method.tau = "DL",
  comb.fixed = FALSE,
  comb.random = FALSE, 
  backtransf = TRUE,
  byvar = subgroup, 
  text.random = "Overall"
)

# Print summary
print(summary(fsw_data_pv_ever_subgroup_forest))

# Save forest plot
png(filename = "Plots/overall plots/ever_pv_subgroup.png", width = 30, height = 30, units = "cm", res = 600)
forest(
  fsw_data_pv_ever_subgroup_forest,
  sortvar = subgroup,
  xlim = c(0.2, 4),             
  leftcols = c("subgroup_level", "studies", "estimates"), 
  leftlabs = c("Category", "Nb studies", "Nb estimates"),
  pooled.totals = FALSE,
  xintercept = 1,
  addrow.overall = TRUE,
  test.subgroup = FALSE,
  overall.hetstat = FALSE,
  overall = FALSE,
  labeltext = TRUE,
  col.subgroup = "black",
  print.subgroup.name = FALSE
)
dev.off(
