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
