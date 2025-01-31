# Subgroup random effects model with constant sampling correlation working model

# covariance matrix assuming constant sampling correlation within subgroups
V_subgroup <- impute_covariance_matrix(fsw_data_pv_recent$effect_best_var_ln, 
                                       cluster = fsw_data_pv_recent$study_num, 
                                       r = rho,
                                       smooth_vi = TRUE,
                                       subgroup = fsw_data_pv_recent$ldc_bin)

# random effects working model in metafor
subgroup_model <- rma.mv(effect_best_ln ~ 0 + ldc_bin + pre_2016 + recruitment + perpetrator + who_region,
                         V = V_subgroup, 
                         random = list(~ ldc_bin | study_num), struct = "DIAG",
                         data = fsw_data_pv_recent, sparse = TRUE)

subgroup_model 

# RVE standard errors
CI_subgroup <- conf_int(subgroup_model, vcov = "CR2")
CI_subgroup

# Robust F-test
Wald_subgroup <- Wald_test(subgroup_model, 
                           constraints = constrain_equal(1:5), 
                           vcov = "CR2")
Wald_subgroup

# Exponentiate the model estimates and CIs to convert them back to the original scale
exp_estimates <- exp(subgroup_model$b)
exp_CI_lower <- exp(subgroup_model$ci.lb)
exp_CI_upper <- exp(subgroup_model$ci.ub)

# Combine the results into a data frame for easier interpretation
results <- data.frame(
  Estimate = exp_estimates,
  CI_Lower = exp_CI_lower,
  CI_Upper = exp_CI_upper
)

print(results)

# Ensure that the categorical variables are factors
fsw_data_pv_recent$ldc_bin <- factor(fsw_data_pv_recent$ldc_bin)
fsw_data_pv_recent$pre_2016 <- factor(fsw_data_pv_recent$pre_2016)
fsw_data_pv_recent$recruitment <- factor(fsw_data_pv_recent$recruitment)
fsw_data_pv_recent$perpetrator <- factor(fsw_data_pv_recent$perpetrator)
fsw_data_pv_recent$who_region <- factor(fsw_data_pv_recent$who_region)

# Check the levels of each factor variable
ldc_bin_levels <- levels(fsw_data_pv_recent$ldc_bin)
pre_2016_levels <- levels(fsw_data_pv_recent$pre_2016)
recruitment_levels <- levels(fsw_data_pv_recent$recruitment)
perpetrator_levels <- levels(fsw_data_pv_recent$perpetrator)
who_region_levels <- levels(fsw_data_pv_recent$who_region)

# Display the reference category for each factor variable
cat("Reference category for ldc_bin:", ldc_bin_levels[1], "\n")
cat("Reference category for pre_2016:", pre_2016_levels[1], "\n")
cat("Reference category for recruitment:", recruitment_levels[1], "\n")
cat("Reference category for perpetrator:", perpetrator_levels[1], "\n")
cat("Reference category for who_region:", who_region_levels[1], "\n")

# Define rho (correlation coefficient)
rho <- 0.6

# List of subgroups and dataframes to loop over
subgroup_vars <- c("ldc_bin", "pre_2016", "recruitment", "perpetrator", "who_region")
dataframe_names <- c("fsw_data_pv_recent", "fsw_data_sv_recent", "fsw_data_psv_recent")

# Function to perform the analysis and save results
perform_analysis <- function(df, dataframe_name, subgroup_var) {
  # Create a covariance matrix assuming constant sampling correlation within subgroups
  V_subgroup <- impute_covariance_matrix(df$effect_best_var_ln, 
                                         cluster = df$study_num, 
                                         r = rho,
                                         smooth_vi = TRUE,
                                         subgroup = df[[subgroup_var]])
  
  # Fit random effects working model in metafor with different optimizer and increased iterations
  subgroup_model <- tryCatch({
    rma.mv(effect_best_ln ~ ldc_bin + pre_2016 + recruitment + perpetrator + who_region,
           V = V_subgroup, 
           random = list(~ get(subgroup_var) | study_num), struct = "DIAG",
           data = df, sparse = TRUE,
           control = list(optimizer = "optim", maxit = 10000))
  }, error = function(e) {
    warning(paste("Model did not converge for", subgroup_var, "in", dataframe_name, ":", e$message))
    return(NULL)
  })
  
  # Check if the model is NULL (did not converge)
  if (is.null(subgroup_model)) {
    return(NULL)
  }
  
  # Exponentiate the model estimates and CIs to convert them back to the original scale
  exp_estimates <- exp(subgroup_model$b)
  exp_CI_lower <- exp(subgroup_model$ci.lb)
  exp_CI_upper <- exp(subgroup_model$ci.ub)
  
  # Combine the results into a data frame for easier interpretation
  results <- data.frame(
    Estimate = exp_estimates,
    CI_Lower = exp_CI_lower,
    CI_Upper = exp_CI_upper
  )
  
  # Save the results to an Excel file with the variable and dataframe names as suffixes
  write.xlsx(results, paste0("results_", subgroup_var, "_", dataframe_name, ".xlsx"))
  
  # Print the results
  print(results)
}

# Loop over each dataframe and each subgroup variable to perform the analysis
for (dataframe_name in dataframe_names) {
  df <- get(dataframe_name)
  for (subgroup_var in subgroup_vars) {
    perform_analysis(df, dataframe_name, subgroup_var)
  }
}


# Define rho (correlation coefficient)
rho <- 0.6

# List of subgroups and dataframes to loop over
subgroup_vars <- c("ldc_bin", "pre_2016", "recruitment", "perpetrator", "who_region")
dataframe_names <- c("fsw_data_pv_recent", "fsw_data_sv_recent", "fsw_data_psv_recent")

# Function to perform the analysis and save results
perform_analysis <- function(df, dataframe_name, subgroup_var) {
  # Create a covariance matrix assuming constant sampling correlation within subgroups
  V_subgroup <- impute_covariance_matrix(df$effect_best_var_ln, 
                                         cluster = df$study_num, 
                                         r = rho,
                                         smooth_vi = TRUE,
                                         subgroup = df[[subgroup_var]])
  
  # Fit random effects working model in metafor with different optimizer and increased iterations
  subgroup_model <- tryCatch({
    rma.mv(effect_best_ln ~ ldc_bin + pre_2016 + recruitment + perpetrator + who_region,
           V = V_subgroup, 
           random = list(~ get(subgroup_var) | study_num), struct = "DIAG",
           data = df, sparse = TRUE,
           control = list(optimizer = "optim", maxit = 10000))
  }, error = function(e) {
    warning(paste("Model did not converge for", subgroup_var, "in", dataframe_name, ":", e$message))
    return(NULL)
  })
  
  # Check if the model is NULL (did not converge)
  if (is.null(subgroup_model)) {
    return(NULL)
  }
  
  # Exponentiate the model estimates and CIs to convert them back to the original scale
  exp_estimates <- exp(subgroup_model$b)
  exp_CI_lower <- exp(subgroup_model$ci.lb)
  exp_CI_upper <- exp(subgroup_model$ci.ub)
  
  # Combine the results into a data frame for easier interpretation
  results <- data.frame(
    Category = rownames(subgroup_model$b),
    Estimate = exp_estimates,
    CI_Lower = exp_CI_lower,
    CI_Upper = exp_CI_upper
  )
  
  # Save the results to an Excel file with the variable and dataframe names as suffixes
  write.xlsx(results, paste0("results_", subgroup_var, "_", dataframe_name, ".xlsx"))
  
  # Print the results
  print(results)
}

# Loop over each dataframe and each subgroup variable to perform the analysis
for (dataframe_name in dataframe_names) {
  df <- get(dataframe_name)
  for (subgroup_var in subgroup_vars) {
    perform_analysis(df, dataframe_name, subgroup_var)
  }
}

