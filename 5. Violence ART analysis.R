# load packages 
pacman::p_load("dplyr", "meta", "metafor", "readxl", "tidyverse", "kableExtra", "robumeta", "clubSandwich") 

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

## ART dataframes
fsw_data_art_uptake <- fsw_data_art %>% filter(outcome_bin == "ART use")
fsw_data_art_adherence <- fsw_data_art %>% filter(outcome_bin == "ART adherence")

#### Multilevel random effects model with constant sampling correlation ####

# Constant sampling correlation
rho <- 0.6

# Define the different analyses
analyses <- c("unadj", "adj", "best")

# Define the corresponding variable names for each analysis
var_names <- list(
  unadj = list(est = "unadj_est_ln", var = "unadj_var_ln", lower = "un_lower_ln", upper = "un_upper_ln"),
  adj = list(est = "adj_est_ln", var = "adj_var_ln", lower = "adj_lower_ln", upper = "adj_upper_ln"),
  best = list(est = "effect_best_ln", var = "effect_best_var_ln", lower = "effect_best_lower_ln", upper = "effect_best_upper_ln")
)

## recent violence data 

# Define the corresponding plot filenames for each analysis
plot_filenames <- list(
  unadj = "Plots/recent_unadj_art_uptake.png",
  adj = "Plots/recent_adj_art_uptake.png",
  best = "Plots/recent_best_art_uptake.png"
)

# Function to perform the analysis and create forest plots
perform_analysis <- function(df, analysis) {
 
  # Filter the dataframe
  filtered_df <- df %>% filter(exposure_tf_bin == "Recent", use == "yes", !is.na(.[[var_names[[analysis]]$est]]))
  
  # Create study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)

  # Create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # Fit a multilevel random effects model using `rma.mv` from metafor
  result <- rma.mv(filtered_df[[var_names[[analysis]]$est]], 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE)       
  
  print(result)
  print(exp(coef(result)))
  
  result2 <- metagen(TE = filtered_df[[var_names[[analysis]]$est]],
                     lower = filtered_df[[var_names[[analysis]]$lower]],
                     upper = filtered_df[[var_names[[analysis]]$upper]],
                     studlab = filtered_df$study,
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
  
  filename <- plot_filenames[[analysis]]
  png(filename = filename, width = 45, height = 14, units = "cm", res = 600)
  
  forest(result2,
         sortvar = filtered_df$study,
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

# Loop over each analysis to perform the analysis and create forest plots
for (analysis in analyses) {
  perform_analysis(fsw_data_art_uptake, analysis)
}

# Define the corresponding plot filenames for each analysis
plot_filenames <- list(
  unadj = "Plots/recent_unadj_art_adherence.png",
  adj = "Plots/recent_adj_art_adherence.png",
  best = "Plots/recent_best_art_adherence.png"
)

# Function to perform the analysis and create forest plots
perform_analysis <- function(df, analysis) {
 
  # Filter the dataframe
  filtered_df <- df %>% filter(exposure_tf_bin == "Recent", use == "yes", !is.na(.[[var_names[[analysis]]$est]]))
  
  # Create study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)

  # Create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # Fit a multilevel random effects model using `rma.mv` from metafor
  result <- rma.mv(filtered_df[[var_names[[analysis]]$est]], 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE)       
  
  print(result)
  print(exp(coef(result)))
  
  result2 <- metagen(TE = filtered_df[[var_names[[analysis]]$est]],
                     lower = filtered_df[[var_names[[analysis]]$lower]],
                     upper = filtered_df[[var_names[[analysis]]$upper]],
                     studlab = filtered_df$study,
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
  
  filename <- plot_filenames[[analysis]]
  png(filename = filename, width = 45, height = 14, units = "cm", res = 600)
  
  forest(result2,
         sortvar = filtered_df$study,
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

# Loop over each analysis to perform the analysis and create forest plots
for (analysis in analyses) {
  perform_analysis(fsw_data_art_adherence, analysis)
}

## lifetime violence data

# Define the corresponding plot filenames for each analysis
plot_filenames <- list(
  unadj = "Plots/ever_unadj_art_uptake.png",
  adj = "Plots/ever_adj_art_uptake.png",
  best = "Plots/ever_best_art_uptake.png"
)

# Function to perform the analysis and create forest plots
perform_analysis <- function(df, analysis) {
  # Filter the dataframe
  filtered_df <- df %>% filter(exposure_tf_bin == "Ever", use == "yes", !is.na(.[[var_names[[analysis]]$est]]))
  
  # Create study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)

  # Create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # Fit a multilevel random effects model using `rma.mv` from metafor
  result <- rma.mv(filtered_df[[var_names[[analysis]]$est]], 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE,
                   control = list(optimizer = "optim", method = "BFGS"))
  
  print(result)
  print(exp(coef(result)))
  
  result2 <- metagen(TE = filtered_df[[var_names[[analysis]]$est]],
                     lower = filtered_df[[var_names[[analysis]]$lower]],
                     upper = filtered_df[[var_names[[analysis]]$upper]],
                     studlab = filtered_df$study,
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
  
  filename <- plot_filenames[[analysis]]
  png(filename = filename, width = 45, height = 14, units = "cm", res = 600)
  
  forest(result2,
         sortvar = filtered_df$study,
         xlim = c(0.2, 4),
         leftcols = leftcols_lifetime,
         leftlabs = leftlabs_lifetime,
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

# Loop over each analysis to perform the analysis and create forest plots for lifetime violence
for (analysis in analyses) {
  perform_analysis(fsw_data_art_uptake, analysis)
}

## lifetime

# Define the corresponding plot filenames for each analysis
plot_filenames <- list(
  unadj = "Plots/ever_unadj_art_adherence.png",
  adj = "Plots/ever_adj_art_adherence.png",
  best = "Plots/ever_best_art_adherence.png"
)

# Function to perform the analysis and create forest plots
perform_analysis <- function(df, analysis) {
  # Filter the dataframe
  filtered_df <- df %>% filter(exposure_tf_bin == "Ever", use == "yes", !is.na(.[[var_names[[analysis]]$est]]))
  
  # Create study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)

  # Create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # Fit a multilevel random effects model using `rma.mv` from metafor
  result <- rma.mv(filtered_df[[var_names[[analysis]]$est]], 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE,
                   control = list(optimizer = "optim", method = "BFGS"))
  
  print(result)
  print(exp(coef(result)))
  
  result2 <- metagen(TE = filtered_df[[var_names[[analysis]]$est]],
                     lower = filtered_df[[var_names[[analysis]]$lower]],
                     upper = filtered_df[[var_names[[analysis]]$upper]],
                     studlab = filtered_df$study,
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
  
  filename <- plot_filenames[[analysis]]
  png(filename = filename, width = 45, height = 14, units = "cm", res = 600)
  
  forest(result2,
         sortvar = filtered_df$study,
         xlim = c(0.2, 4),
         leftcols = leftcols_lifetime,
         leftlabs = leftlabs_lifetime,
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

# Loop over each analysis to perform the analysis and create forest plots for lifetime violence
for (analysis in analyses) {
  perform_analysis(fsw_data_art_adherence, analysis)
}

## combined recent and lifetime

# Define the corresponding plot filenames for each analysis
plot_filenames <- list(
  unadj = "Plots/all_unadj_art_uptake.png",
  adj = "Plots/all_adj_art_uptake.png",
  best = "Plots/all_best_art_uptake.png"
)

# Function to perform the analysis and create forest plots without outlier
perform_analysis <- function(df, analysis) {
 
  # Filter the dataframe
  filtered_df <- df %>% filter(use == "yes", !is.na(.[[var_names[[analysis]]$est]]))
  
  # Create study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)

  # Create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # Fit a multilevel random effects model using `rma.mv` from metafor
  result <- rma.mv(filtered_df[[var_names[[analysis]]$est]], 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE,
                   control = list(optimizer = "optim", method = "BFGS"))
  
  print(result)
  print(exp(coef(result)))
  
  result2 <- metagen(TE = filtered_df[[var_names[[analysis]]$est]],
                     lower = filtered_df[[var_names[[analysis]]$lower]],
                     upper = filtered_df[[var_names[[analysis]]$upper]],
                     studlab = filtered_df$study,
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
  
  filename <- plot_filenames[[analysis]]
  png(filename = filename, width = 45, height = 14, units = "cm", res = 600)
  
  forest(result2,
         sortvar = filtered_df$study,
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

# Loop over each analysis to perform the analysis and create forest plots
for (analysis in analyses) {
  perform_analysis(fsw_data_art_uptake, analysis)
}


# Define the corresponding plot filenames for each analysis
plot_filenames <- list(
  unadj = "Plots/all_unadj_art_adherence.png",
  adj = "Plots/all_adj_art_adherence.png",
  best = "Plots/all_best_art_adherence.png"
)

# Function to perform the analysis and create forest plots without outlier
perform_analysis <- function(df, analysis) {
 
  # Filter the dataframe
  filtered_df <- df %>% filter(use == "yes", !is.na(.[[var_names[[analysis]]$est]]))
  
  # Create study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)

  # Create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # Fit a multilevel random effects model using `rma.mv` from metafor
  result <- rma.mv(filtered_df[[var_names[[analysis]]$est]], 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE,
                   control = list(optimizer = "optim", method = "BFGS"))
  
  print(result)
  print(exp(coef(result)))
  
  result2 <- metagen(TE = filtered_df[[var_names[[analysis]]$est]],
                     lower = filtered_df[[var_names[[analysis]]$lower]],
                     upper = filtered_df[[var_names[[analysis]]$upper]],
                     studlab = filtered_df$study,
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
  
  filename <- plot_filenames[[analysis]]
  png(filename = filename, width = 45, height = 14, units = "cm", res = 600)
  
  forest(result2,
         sortvar = filtered_df$study,
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

# Loop over each analysis to perform the analysis and create forest plots
for (analysis in analyses) {
  perform_analysis(fsw_data_art_adherence, analysis)
}



























































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
