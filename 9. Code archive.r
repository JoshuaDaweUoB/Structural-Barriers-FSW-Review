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

# load packages 
pacman::p_load("meta", "metafor", "readxl", "tidyverse") 

# settings
settings.meta(CIbracket = "(") 
settings.meta(CIseparator = "-") 

# who subgroups
for (var in dfs_subgroup) {
  df <- get(var)
  meta_analysis_best_who(df, var)
}

for (var in dfs_subgroup) {
  df <- get(var)
  meta_analysis_best_who_sv(df, var)
}

# lic subgroups
for (var in dfs_subgroup) {
  df <- get(var)
  meta_analysis_best_ldc(df, var)
}

# recruitment subgroups
for (var in dfs_subgroup) {
  df <- get(var)
  meta_analysis_best_recruit(df, var)
}

# study setting subgroups
for (var in dfs_subgroup) {
  df <- get(var)
  meta_analysis_best_setting(df, var)
}

# perpetrator subgroups
for (var in dfs_subgroup) {
  df <- get(var)
  meta_analysis_best_perp(df, var)
}

# perpetrator subgroups
for (var in dfs_subgroup_perp) {
  df <- get(var)
  meta_analysis_best_perp_all(df, var)
}

# pre 2016 subgroups
for (var in dfs_subgroup) {
  df <- get(var)
  meta_analysis_best_pre2016(df, var)
}

# create sexual violence subgroup forest plots
for (df_name in dfs_subgroups) {
  subgroup_plot_forest_sv(df_name)
}

# create physical violence subgroup forest plots
for (df_name in dfs_subgroups) {
  subgroup_plot_forest_pv(df_name)
}

# create sexual or physical violence subgroup forest plots
for (df_name in dfs_subgroups) {
  subgroup_plot_forest_psv(df_name)
}
# load packages 
pacman::p_load("meta", "metafor", "readxl", "tidyverse") 

# settings
settings.meta(CIbracket = "(") 
settings.meta(CIseparator = "-") 

# ever summary plot
summary_hiv_violence_ever  <- metagen(TE = effect_ln,
                                      lower = lower_ln,
                                      upper = upper_ln,
                                      studlab = name,
                                      data = summary_violence_ever,
                                      sm = "OR",
                                      method.tau = "DL",
                                      comb.fixed = FALSE,
                                      comb.random = FALSE, 
                                      backtransf = TRUE,
                                      byvar = Adjust,
                                      text.random = "Overall")

summary(summary_hiv_violence_ever) 

png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/summary_hiv_violence_ever.png", width = 24, height = 12, units = "cm",
    res = 600)

summary_hiv_violence_ever <- forest.meta(summary_hiv_violence_ever, 
                                         sortvar = name,
                                         xlim=c(0.2, 4),             
                                         leftcols = c("name", "studies"), 
                                         leftlabs = c("Pooled exposure", "Studies"),
                                         rightcols = c("or_95_2", "fuck"), 
                                         rightlabs = c("OR (95% CI)", "I²"), 
                                         pooled.totals = F,
                                         xintercept=1,
                                         addrow.overall = T,
                                         test.subgroup = F,
                                         overall.hetstat = F,
                                         overall = F,
                                         labeltext = TRUE,
                                         col.subgroup = "black",
                                         print.subgroup.name = FALSE) 
dev.off()

## recent summary plot
summary_hiv_violence_rec  <- metagen(TE = effect_ln,
                                     lower = lower_ln,
                                     upper = upper_ln,
                                     studlab = name,
                                     data = summary_violence_rec,
                                     sm = "OR",
                                     method.tau = "DL",
                                     comb.fixed = FALSE,
                                     comb.random = FALSE, 
                                     backtransf = TRUE,
                                     byvar = Adjust,
                                     text.random = "Overall")

summary(summary_hiv_violence_rec) 

png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/summary_hiv_violence_recent.png", width = 24, height = 12, units = "cm",
    res = 600)

summary_hiv_violence_rec <- forest.meta(summary_hiv_violence_rec, 
                                        sortvar = name,
                                        xlim=c(0.2, 4),             
                                        leftcols = c("name", "studies"), 
                                        leftlabs = c("Pooled exposure", "Studies"),
                                        rightcols = c("or_95_2", "fuck"), 
                                        rightlabs = c("OR (95% CI)", "I²"), 
                                        pooled.totals = F,
                                        xintercept=1,
                                        addrow.overall = T,
                                        test.subgroup = F,
                                        overall.hetstat = F,
                                        overall = F,
                                        labeltext = TRUE,
                                        col.subgroup = "black",
                                        print.subgroup.name = FALSE) 
dev.off()

# all physical violence and hiv infection data


# all ART data
summary_hiv_art  <- metagen(TE = effect_best_ln,
                            lower = effect_best_lower_ln,
                            upper = effect_best_upper_ln,
                            studlab = study,
                            data = fsw_data_art,
                            sm = "OR",
                            method.tau = "DL",
                            comb.fixed = FALSE,
                            comb.random = FALSE, 
                            backtransf = TRUE,
                            byvar = exposure_type,
                            text.random = "Overall")

summary(summary_hiv_art) 

png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/hiv_art.png", width = 20, height = 18, units = "cm",
    res = 500)

summary_hiv_art <- forest.meta(summary_hiv_art, 
                               sortvar = study,
                               xlim=c(0.2, 4),             
                               leftcols = c("study", "outcome"), 
                               leftlabs = c("Study", "Outcome"),
                               rightcols = c("exposure_time_frame", "perpetrator"), 
                               rightlabs = c("Time frame", "Perpetrator"), 
                               pooled.totals = F,
                               xintercept=1,
                               addrow.overall = T,
                               test.subgroup = F,
                               overall.hetstat = F,
                               overall = F,
                               labeltext = TRUE,
                               col.subgroup = "black",
                               print.subgroup.name = FALSE) 
dev.off()

# all VS data
summary_hiv_vs  <- metagen(TE = effect_best_ln,
                           lower = effect_best_lower_ln,
                           upper = effect_best_upper_ln,
                           studlab = study,
                           data = fsw_data_vs,
                           sm = "OR",
                           method.tau = "DL",
                           comb.fixed = FALSE,
                           comb.random = FALSE, 
                           backtransf = TRUE,
                           byvar = exposure_type,
                           text.random = "Overall")

summary(summary_hiv_vs) 

png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/hiv_vs.png", width = 22, height = 10, units = "cm",
    res = 500)

summary_hiv_vs <- forest.meta(summary_hiv_vs, 
                              sortvar = study,
                              xlim=c(0.2, 4),             
                              leftcols = c("study", "outcome_definition"), 
                              leftlabs = c("Study", "VS definition"),
                              rightcols = c("exposure_time_frame", "perpetrator"), 
                              rightlabs = c("Time frame", "Perpetrator"), 
                              pooled.totals = F,
                              xintercept=1,
                              addrow.overall = T,
                              test.subgroup = F,
                              overall.hetstat = F,
                              overall = F,
                              labeltext = TRUE,
                              col.subgroup = "black",
                              print.subgroup.name = FALSE) 
dev.off()

# SV and HIV data
summary_hiv_sv  <- metagen(TE = effect_best_ln,
                           lower = effect_best_lower_ln,
                           upper = effect_best_upper_ln,
                           studlab = study,
                           data = fsw_data_sv_hiv,
                           sm = "OR",
                           method.tau = "DL",
                           comb.fixed = FALSE,
                           comb.random = FALSE, 
                           backtransf = TRUE,
                           text.random = "Overall")

summary(summary_hiv_sv) 

png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/hiv_sv.png", width = 25, height = 25, units = "cm",
    res = 500)

summary_hiv_sv <- forest.meta(summary_hiv_sv, 
                              sortvar = study,
                              xlim=c(0.2, 4),             
                              leftcols = c("study"), 
                              leftlabs = c("Study"),
                              rightcols = c("exposure_time_frame", "perpetrator"), 
                              rightlabs = c("Time frame", "Perpetrator"), 
                              pooled.totals = F,
                              xintercept=1,
                              addrow.overall = T,
                              test.subgroup = F,
                              overall.hetstat = F,
                              overall = F,
                              labeltext = TRUE,
                              col.subgroup = "black",
                              print.subgroup.name = FALSE) 
dev.off()










#### multilevel random effects model with constant sampling correlation ####

# constant sampling correlation
rho <- 0.6

# Unadjusted analyses with subgroups for exposure_type
{
  # Filter the dataframe 
  filtered_df <- fsw_data_art %>% 
    filter(exposure_tf_bin == "Recent") %>% 
    filter(unadj_est != "NA") %>% 
    filter(outcome_bin == "ART use")
  
   # Create study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)

  # Get unique levels of the exposure_type variable
  unique_exposure_types <- unique(filtered_df$exposure_type)
  
  # Initialize lists to store results
  metagen_list <- list()
  subgroup_labels <- list()
  study_labels <- list()
  outcome_labels <- list()
  exposure_definitions <- list()
  exposure_time_frames <- list()
  perpetrators <- list()
  countries <- list()
  
  # Loop through each unique exposure type
  for (exposure in unique_exposure_types) {
    cat("Processing exposure type:", exposure, "\n")
    
    # Filter data for the current exposure type
    subgroup_df <- filtered_df %>% filter(exposure_type == exposure)
    
    # Debug: Print the number of rows for the current subgroup
    cat("Number of rows for", exposure, ":", nrow(subgroup_df), "\n")
    
    # Create a covariance matrix assuming constant sampling correlation
    V_mat <- impute_covariance_matrix(subgroup_df$unadj_var_ln,
                                      cluster = subgroup_df$study_num,
                                      r = rho,
                                      smooth_vi = TRUE)
    
    # Fit a multilevel random effects model using `rma.mv` from metafor
    result <- rma.mv(unadj_est_ln, 
                     V = V_mat, 
                     random = ~ 1 | study_num / effect_num,
                     data = subgroup_df,   
                     sparse = TRUE)       
    
    cat("Model fitting complete for exposure type:", exposure, "\n")
    
    result 
    exp(coef(result))
    
    result2 <- metagen(TE = subgroup_df$unadj_est_ln,
                       lower = subgroup_df$un_lower_ln,
                       upper = subgroup_df$un_upper_ln,
                       studlab = subgroup_df$study,
                       data = subgroup_df,
                       sm = "OR",
                       method.tau = "REML",
                       common = FALSE,
                       random = TRUE,
                       backtransf = TRUE)
    
    print(summary(result2))
    
    # Debug: Print the length of the current metagen object
    cat("Length of TE for", exposure, ":", length(result2$TE), "\n")
    
    # Store the metagen object and additional labels
    metagen_list[[length(metagen_list) + 1]] <- result2
    subgroup_labels[[length(subgroup_labels) + 1]] <- rep(exposure, length(result2$TE))
    study_labels[[length(study_labels) + 1]] <- subgroup_df$study
    outcome_labels[[length(outcome_labels) + 1]] <- subgroup_df$outcome_definition_short
    exposure_definitions[[length(exposure_definitions) + 1]] <- subgroup_df$exposure_definition_short
    exposure_time_frames[[length(exposure_time_frames) + 1]] <- subgroup_df$exposure_time_frame
    perpetrators[[length(perpetrators) + 1]] <- subgroup_df$perpetrator
    countries[[length(countries) + 1]] <- subgroup_df$country
  }
  
  # Combine the metagen objects into one
  combined_TE <- unlist(lapply(metagen_list, function(x) x$TE))
  combined_lower <- unlist(lapply(metagen_list, function(x) x$lower))
  combined_upper <- unlist(lapply(metagen_list, function(x) x$upper))
  combined_studlab <- unlist(lapply(metagen_list, function(x) x$studlab))
  combined_subgroup <- unlist(subgroup_labels)
  combined_study <- unlist(study_labels)
  combined_outcome <- unlist(outcome_labels)
  combined_exposure_definition <- unlist(exposure_definitions)
  combined_exposure_time_frame <- unlist(exposure_time_frames)
  combined_perpetrator <- unlist(perpetrators)
  combined_country <- unlist(countries)
  
  # Debug: Print lengths of combined vectors
  cat("Length of combined_TE:", length(combined_TE), "\n")
  cat("Length of combined_studlab:", length(combined_studlab), "\n")
  cat("Length of combined_subgroup:", length(combined_subgroup), "\n")
  cat("Length of combined_study:", length(combined_study), "\n")
  cat("Length of combined_outcome:", length(combined_outcome), "\n")
  cat("Length of combined_exposure_definition:", length(combined_exposure_definition), "\n")
  cat("Length of combined_exposure_time_frame:", length(combined_exposure_time_frame), "\n")
  cat("Length of combined_perpetrator:", length(combined_perpetrator), "\n")
  cat("Length of combined_country:", length(combined_country), "\n")
  
  # Ensure lengths match before creating the combined metagen object
  if (length(combined_TE) == length(combined_studlab) && length(combined_TE) == length(combined_subgroup) && length(combined_TE) == length(combined_study) && length(combined_TE) == length(combined_outcome) && length(combined_TE) == length(combined_exposure_definition) && length(combined_TE) == length(combined_exposure_time_frame) && length(combined_TE) == length(combined_perpetrator) && length(combined_TE) == length(combined_country)) {
    combined_data <- data.frame(
      study = combined_study,
      subgroup = combined_subgroup,
      outcome_definition_short = combined_outcome,
      exposure_definition_short = combined_exposure_definition,
      exposure_time_frame = combined_exposure_time_frame,
      perpetrator = combined_perpetrator,
      country = combined_country
    )
    
    combined_metagen <- metagen(
      TE = combined_TE,
      lower = combined_lower,
      upper = combined_upper,
      studlab = combined_data$study,
      subgroup = combined_data$subgroup,
      sm = "OR",
      method.tau = "REML",
      common = FALSE,
      random = TRUE,
      backtransf = TRUE
    )
    
    # Sanitize the filename
    filename <- "Plots/recent_unadj_art_combined.png"
    
    cat("Saving combined plot to:", filename, "\n")
    
    png(filename = filename, width = 50, height = 18, units = "cm", res = 600)
    
    forest(combined_metagen,
           xlim = c(0.2, 4),             
           leftcols = c("combined_studlab", "combined_exposure_definition", "combined_exposure_time_frame", "combined_perpetrator", "combined_country"), 
           leftlabs = c("Study", "Exposure definition", "Exposure time frame", "Perpetrator", "Country"),
           rightcols = c("combined_outcome"),
           rightlabs = c("Outcome"),
           pooled.totals = TRUE,
           xintercept = 1,
           addrow.overall = TRUE,
           overall.hetstat = TRUE,
           overall = TRUE,
           labeltext = TRUE,
           col.subgroup = "black")
    
    dev.off()
    
    cat("Combined plot saved.\n")
  } else {
    cat("Error: Lengths of combined vectors do not match.\n")
  }
}











































































































































# Unadjusted analyses with subgroups for exposure_type
{
  # Filter the dataframe for the desired conditions
  filtered_df <- fsw_data_art %>% 
    filter(exposure_tf_bin == "Recent") %>% 
    filter(unadj_est != "NA") %>% 
    filter(outcome_bin == "ART use")
  
  # Get unique levels of the exposure_type variable
  unique_exposure_types <- unique(filtered_df$exposure_type)
  
  # Loop through each unique exposure type
  for (exposure in unique_exposure_types) {
    # Filter data for the current exposure type
    subgroup_df <- filtered_df %>% filter(exposure_type == exposure)
    
    # Create a covariance matrix assuming constant sampling correlation
    V_mat <- impute_covariance_matrix(subgroup_df$unadj_var_ln,
                                      cluster = subgroup_df$study_num,
                                      r = rho,
                                      smooth_vi = TRUE)
    
    # Fit a multilevel random effects model using `rma.mv` from metafor
    result <- rma.mv(unadj_est_ln, 
                     V = V_mat, 
                     random = ~ 1 | study_num / effect_num,
                     data = subgroup_df,   
                     sparse = TRUE)       
    
    result 
    exp(coef(result))
    
    result2 <- metagen(TE = unadj_est_ln,
                       lower = un_lower_ln,
                       upper = un_upper_ln,
                       studlab = study,
                       data = subgroup_df,
                       sm = "OR",
                       method.tau = "REML",
                       common = FALSE,
                       random = TRUE, 
                       backtransf = TRUE,
                       text.random = "Overall")
    
    print(summary(result2))
    
    result2$TE.random <- result$b
    result2$lower.random <- result$ci.lb
    result2$upper.random <- result$ci.ub
    
    # Sanitize the filename
    sanitized_exposure <- gsub("[^a-zA-Z0-9]", "_", exposure)
    filename <- paste0("Plots/recent_unadj_art_", sanitized_exposure, ".png")
    
    png(filename = filename, width = 50, height = 18, units = "cm", res = 600)
    
    forest(result2,
           sortvar = study,
           xlim = c(0.2, 4),             
           leftcols = leftcols_recent, 
           leftlabs = leftlabs_recent,
           rightcols = rightcols,
           rightlabs = rightlabs,
           pooled.totals = TRUE,
           xintercept = 1,
           addrow.overall = TRUE,
           overall.hetstat = TRUE,
           overall = TRUE,
           labeltext = TRUE,
           col.subgroup = "black")
    
    dev.off()
  }
}

## unadjusted analyses ##

{ 
  
  filtered_df <- fsw_data_art %>% filter(exposure_tf_bin == "Recent")
  filtered_df <- filtered_df %>% filter(unadj_est != "NA")
  filtered_df <- filtered_df %>% filter(outcome_bin == "ART use")

  # create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df$unadj_var_ln,
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # fit a multilevel random effects model using `rma.mv` from metafor
  result <- rma.mv(unadj_est_ln, 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE)       
  
  result 
  exp(coef(result))
  
  result2 <- metagen(TE = unadj_est_ln,
                     lower = un_lower_ln,
                     upper = un_upper_ln,
                     studlab = study,
                     data = filtered_df,
                     sm = "OR",
                     method.tau = "REML",
                     common = FALSE,
                     random = TRUE, 
                     backtransf = TRUE,
                     text.random = "Overall")
  
  print(summary(result2))
  
  result2$TE.random <- result$b
  result2$lower.random <- result$ci.lb
  result2$upper.random <- result$ci.ub
  
  filename <- paste0("Plots/recent_unadj_art.png")
  png(filename = filename, width = 50, height = 18, units = "cm", res = 600)
  
  forest(result2,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols_recent, 
         leftlabs = leftlabs_recent,
         rightcols = rightcols,
         rightlabs = rightlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  dev.off()  
  
}  

## adjusted analyses ##

{ 
  filtered_df <- fsw_data_art %>% filter(exposure_tf_bin == "Recent")
  filtered_df <- filtered_df %>% filter(adj_est != "NA")
  filtered_df <- filtered_df %>% filter(outcome_bin == "ART use")

  # create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df$adj_var_ln,
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # fit a multilevel random effects model using `rma.mv` from metafor
  result <- rma.mv(adj_est_ln, 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE)       
  
  result 
  exp(coef(result))
  
  
  result2 <- metagen(TE = adj_est_ln,
                     lower = adj_lower_ln,
                     upper = adj_upper_ln,
                     studlab = study,
                     data = filtered_df,
                     sm = "OR",
                     method.tau = "REML",
                     common = FALSE,
                     random = TRUE, 
                     backtransf = TRUE,
                     text.random = "Overall")
  
  print(summary(result2))
  
  result2$TE.random <- result$b
  result2$lower.random <- result$ci.lb
  result2$upper.random <- result$ci.ub
  
  filename <- paste0("Plots/recent_adj_art.png")
  png(filename = filename, width = 50, height = 16, units = "cm", res = 600)
  
  forest(result2,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols_recent, 
         leftlabs = leftlabs_recent,
         rightcols = rightcols,
         rightlabs = rightlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  dev.off()  
  
}  

## best analyses ## 

{ 
  filtered_df <- fsw_data_art %>% filter(exposure_tf_bin == "Recent")
  filtered_df <- filtered_df %>% filter(effect_best != "NA")
  filtered_df <- filtered_df %>% filter(outcome_bin == "ART use")

  # create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df$effect_best_var_ln,
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # fit a multilevel random effects model using `rma.mv` from metafor
  result <- rma.mv(effect_best_ln, 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE)       
  
  result 
  exp(coef(result))
  
  
  result2 <- metagen(TE = effect_best_ln,
                     lower = effect_best_lower_ln,
                     upper = effect_best_upper_ln,
                     studlab = study,
                     data = filtered_df,
                     sm = "OR",
                     method.tau = "REML",
                     common = FALSE,
                     random = TRUE, 
                     backtransf = TRUE,
                     text.random = "Overall")
  
  print(summary(result2))
  
  result2$TE.random <- result$b
  result2$lower.random <- result$ci.lb
  result2$upper.random <- result$ci.ub
  
  filename <- paste0("Plots/recent_best_art.png")
  png(filename = filename, width = 50, height = 18, units = "cm", res = 600)
  
  forest(result2,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols_recent, 
         leftlabs = leftlabs_recent,
         rightcols = rightcols,
         rightlabs = rightlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  dev.off()  
  
}  

### lifetime violence (all types) ###

## unadjusted analyses ##

{ 
  
  filtered_df <- fsw_data_art %>% filter(exposure_tf_bin == "Ever")
  filtered_df <- filtered_df %>% filter(unadj_est != "NA")
  filtered_df <- filtered_df %>% filter(outcome_bin == "ART use")

  # create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df$unadj_var_ln,
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # fit a multilevel random effects model using `rma.mv` from metafor
  result <- rma.mv(unadj_est_ln, 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE)       
  
  result 
  exp(coef(result))
  
  
  result2 <- metagen(TE = unadj_est_ln,
                     lower = un_lower_ln,
                     upper = un_upper_ln,
                     studlab = study,
                     data = filtered_df,
                     sm = "OR",
                     method.tau = "REML",
                     common = FALSE,
                     random = TRUE, 
                     backtransf = TRUE,
                     text.random = "Overall")
  
  print(summary(result2))
  
  result2$TE.random <- result$b
  result2$lower.random <- result$ci.lb
  result2$upper.random <- result$ci.ub
  
  filename <- paste0("Plots/ever_unadj_art.png")
  png(filename = filename, width = 38, height = 14, units = "cm", res = 600)
  
  forest(result2,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols_recent, 
         leftlabs = leftlabs_recent,
         rightcols = rightcols,
         rightlabs = rightlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  dev.off()  
  
}  

## adjusted analyses ##

{ 
  filtered_df <- fsw_data_art %>% filter(exposure_tf_bin == "Ever")
  filtered_df <- filtered_df %>% filter(adj_est != "NA")
  filtered_df <- filtered_df %>% filter(outcome_bin == "ART use")

  # create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df$adj_var_ln,
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # fit a multilevel random effects model using `rma.mv` from metafor
  result <- rma.mv(adj_est_ln, 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE)       
  
  result 
  exp(coef(result))
  
  
  result2 <- metagen(TE = adj_est_ln,
                     lower = adj_lower_ln,
                     upper = adj_upper_ln,
                     studlab = study,
                     data = filtered_df,
                     sm = "OR",
                     method.tau = "REML",
                     common = FALSE,
                     random = TRUE, 
                     backtransf = TRUE,
                     text.random = "Overall")
  
  print(summary(result2))
  
  result2$TE.random <- result$b
  result2$lower.random <- result$ci.lb
  result2$upper.random <- result$ci.ub
  
  filename <- paste0("Plots/ever_adj_art.png")
  png(filename = filename, width = 38, height = 14, units = "cm", res = 600)
  
  forest(result2,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols_recent, 
         leftlabs = leftlabs_recent,
         rightcols = rightcols,
         rightlabs = rightlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  dev.off()  
  
}  

## best analyses ## 

{ 
  filtered_df <- fsw_data_art %>% filter(exposure_tf_bin == "Ever")
  filtered_df <- filtered_df %>% filter(effect_best != "NA")
  filtered_df <- filtered_df %>% filter(outcome_bin == "ART use")
 
  # create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df$effect_best_var_ln,
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # fit a multilevel random effects model using `rma.mv` from metafor
  result <- rma.mv(effect_best_ln, 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE)       
  
  result 
  exp(coef(result))
  
  
  result2 <- metagen(TE = effect_best_ln,
                     lower = effect_best_lower_ln,
                     upper = effect_best_upper_ln,
                     studlab = study,
                     data = filtered_df,
                     sm = "OR",
                     method.tau = "REML",
                     common = FALSE,
                     random = TRUE, 
                     backtransf = TRUE,
                     text.random = "Overall")
  
  print(summary(result2))
  
  result2$TE.random <- result$b
  result2$lower.random <- result$ci.lb
  result2$upper.random <- result$ci.ub
  
  filename <- paste0("Plots/ever_best_art.png")
  png(filename = filename, width = 38, height = 14, units = "cm", res = 600)
  
  forest(result2,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols_recent, 
         leftlabs = leftlabs_recent,
         rightcols = rightcols,
         rightlabs = rightlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  dev.off()  
  
}  

#### overall forest plot ####

## load dataframe

summary_violence_art <- read_excel("Violence estimates.xlsx", "ART") 

## forest plot

summary_hiv_violence_art  <- metagen(TE = effect_ln,
                                      lower = lower_ln,
                                      upper = upper_ln,
                                      studlab = name,
                                      data = summary_violence_art,
                                      sm = "OR",
                                      method.tau = "DL",
                                      comb.fixed = FALSE,
                                      comb.random = FALSE, 
                                      backtransf = TRUE,
                                      byvar = name,
                                      text.random = "Overall")

summary(summary_hiv_violence_art) 

filename <- paste0("Plots/overall plots/violence_art_overall.png")
png(filename = filename, width = 25, height = 14, units = "cm", res = 600)

summary_hiv_violence_art <- forest(summary_hiv_violence_art, 
                                    sortvar = name,
                                    xlim=c(0.2, 4),             
                                    leftcols = c("adjust", "studies", "estimates"), 
                                    leftlabs = c("Model type", "Studies", "Estimates"),
                                    rightcols = c("or_95", "i2"), 
                                    rightlabs = c("OR (95% CI)", "I²"), 
                                    pooled.totals = F,
                                    xintercept=1,
                                    addrow.overall = T,
                                    test.subgroup = F,
                                    overall.hetstat = F,
                                    overall = F,
                                    labeltext = TRUE,
                                    col.subgroup = "black",
                                    print.subgroup.name = FALSE) 
dev.off()
