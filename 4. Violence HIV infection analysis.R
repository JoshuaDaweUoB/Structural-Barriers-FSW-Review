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
dataframes <- list(fsw_data_pv_recent, fsw_data_sv_recent, fsw_data_psv_recent, fsw_data_pv_ever, fsw_data_sv_ever, fsw_data_psv_ever)  
dataframe_names <- c("fsw_data_pv_recent", "fsw_data_sv_recent", "fsw_data_psv_recent", "fsw_data_pv_ever", "fsw_data_sv_ever", "fsw_data_psv_ever")
analyses <- c("unadj", "adj", "best")
exposures <- c("recent", "ever")

#### multilevel random effects model with constant sampling correlation ####

# constant sampling correlation
rho <- 0.6

# Define the corresponding variable names for each analysis
var_names <- list(
  unadj = list(est = "unadj_est_ln", var = "unadj_var_ln", lower = "un_lower_ln", upper = "un_upper_ln"),
  adj = list(est = "adj_est_ln", var = "adj_var_ln", lower = "adj_lower_ln", upper = "adj_upper_ln"),
  best = list(est = "effect_best_ln", var = "effect_best_var_ln", lower = "effect_best_lower_ln", upper = "effect_best_upper_ln")
)

# Define the corresponding dataframes for each exposure and type of violence
dataframes <- list(
  physical = list(recent = "fsw_data_pv_recent", ever = "fsw_data_pv_ever"),
  sexual = list(recent = "fsw_data_sv_recent", ever = "fsw_data_sv_ever"),
  physical_sexual = list(recent = "fsw_data_psv_recent", ever = "fsw_data_psv_ever")
)

# Define the corresponding plot filenames for each analysis, exposure, and type of violence
plot_filenames <- list(
  physical = list(
    unadj = list(recent = "Plots/pv_recent_unadj_csc.png", ever = "Plots/pv_ever_unadj_csc.png"),
    adj = list(recent = "Plots/pv_recent_adj_csc.png", ever = "Plots/pv_ever_adj_csc.png"),
    best = list(recent = "Plots/pv_recent_best_csc.png", ever = "Plots/pv_ever_best_csc.png")
  ),
  sexual = list(
    unadj = list(recent = "Plots/sv_recent_unadj_csc.png", ever = "Plots/sv_ever_unadj_csc.png"),
    adj = list(recent = "Plots/sv_recent_adj_csc.png", ever = "Plots/sv_ever_adj_csc.png"),
    best = list(recent = "Plots/sv_recent_best_csc.png", ever = "Plots/sv_ever_best_csc.png")
  ),
  physical_sexual = list(
    unadj = list(recent = "Plots/psv_recent_unadj_csc.png", ever = "Plots/psv_ever_unadj_csc.png"),
    adj = list(recent = "Plots/psv_recent_adj_csc.png", ever = "Plots/psv_ever_adj_csc.png"),
    best = list(recent = "Plots/psv_recent_best_csc.png", ever = "Plots/psv_ever_best_csc.png")
  )
)

# Define the corresponding left columns and labels for each exposure
leftcols <- list(
  recent = leftcols_recent,
  ever = leftcols_lifetime
)
leftlabs <- list(
  recent = leftlabs_recent,
  ever = leftlabs_lifetime
)

# Function to perform the analysis and create forest plots
perform_analysis <- function(df, analysis, exposure, violence_type) {
  # Filter the dataframe
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  filtered_df <- filtered_df %>% filter(!is.na(filtered_df[[var_names[[analysis]]$est]]))
  
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
  
  filename <- plot_filenames[[violence_type]][[analysis]][[exposure]]
  png(filename = filename, width = 45, height = 22, units = "cm", res = 600)
  
  forest(result2,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols[[exposure]], 
         leftlabs = leftlabs[[exposure]],
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

# Loop over each type of violence, analysis, and exposure to perform the analysis and create forest plots
for (violence_type in names(dataframes)) {
  for (analysis in analyses) {
    for (exposure in exposures) {
      df <- get(dataframes[[violence_type]][[exposure]])
      perform_analysis(df, analysis, exposure, violence_type)
    }
  }
}

## other violence types

# Function to create forest plots for other violence types
create_other_violence_plots <- function(data, outcome_filter = "HIV prevalence") {
  
  # Loop through each analysis type (unadjusted, adjusted, best)
  for (analysis in analyses) {
    
    # Filter the dataframe
    filtered_df <- data %>% filter(outcome == outcome_filter, exposure_tf_bin == "Recent")
    filtered_df <- filtered_df %>% filter(!is.na(.data[[var_names[[analysis]]$est]]))

    # Skip if no data after filtering
    if (nrow(filtered_df) == 0) {
      message(paste("No data available for", analysis, "analysis"))
      next
    }
    
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
    
    print(paste("Results for", analysis, "analysis:"))
    print(result)
    print(exp(coef(result)))
    
    # Create metagen object for forest plot
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
    
    # Update the metagen object with rma.mv results
    result2$TE.random <- result$b
    result2$lower.random <- result$ci.lb
    result2$upper.random <- result$ci.ub
    
    # Create forest plot filename
    filename <- paste0("Plots/hiv_infection_other_violence_", analysis, ".png")
    
    # Create forest plot
    png(filename = filename, width = 45, height = 14, units = "cm", res = 600)
    
    forest(result2,
           sortvar = filtered_df$study,
           xlim = c(0.2, 4),             
           leftcols = c("studlab", "exposure_definition_short", "exposure_time_frame", "perpetrator", "country"), 
           leftlabs = c("Study", "Exposure definition", "Exposure time frame", "Perpetrator", "Country"),
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
    
    message(paste("Forest plot saved as:", filename))
  }
}

# Call the function to create all three plots
create_other_violence_plots(fsw_data_other)

## any violence forest plots

#### overall forest plots ####

# Function to create forest plots
create_forest_plot <- function(data, effect_col, lower_col, upper_col, studlab_col, byvar_col, filename) {
  # Perform meta-analysis
  summary_hiv_violence <- metagen(TE = data[[effect_col]],
                                  lower = data[[lower_col]],
                                  upper = data[[upper_col]],
                                  studlab = data[[studlab_col]],
                                  data = data,
                                  sm = "OR",
                                  method.tau = "DL",
                                  comb.fixed = FALSE,
                                  comb.random = FALSE, 
                                  backtransf = TRUE,
                                  byvar = data[[byvar_col]],
                                  text.random = "Overall")
  
  # Print summary
  print(summary(summary_hiv_violence))
  
  # Save forest plot
  png(filename = filename, width = 25, height = 14, units = "cm", res = 600)
  forest(summary_hiv_violence, 
         sortvar = data[[studlab_col]],
         xlim = c(0.2, 4),             
         leftcols = c("name", "studies", "estimates"), 
         leftlabs = c("Pooled exposure", "Studies", "Estimates"),
         rightcols = c("or_95_2", "i2"), 
         rightlabs = c("OR (95% CI)", "I²"), 
         pooled.totals = FALSE,
         xintercept = 1,
         addrow.overall = TRUE,
         test.subgroup = FALSE,
         overall.hetstat = FALSE,
         overall = FALSE,
         labeltext = TRUE,
         col.subgroup = "black",
         print.subgroup.name = FALSE)
  dev.off()
}

# Load dataframes
summary_violence_ever <- read_excel("Violence estimates.xlsx", "HIV infection - ever")
summary_violence_recent <- read_excel("Violence estimates.xlsx", "HIV infection - recent")

# Create forest plots
create_forest_plot(summary_violence_ever, "effect_ln_2", "lower_ln_2", "upper_ln_2", "name", "Adjust", "Plots/overall plots/violence_ever_overall.png")
create_forest_plot(summary_violence_recent, "effect_ln_2", "lower_ln_2", "upper_ln_2", "name", "Adjust", "Plots/overall plots/violence_recent_overall.png")

## subgroup analysis 

# Call the function for "recently exposed to physical violence"
process_and_plot(
  data = fsw_data_pv_recent,
  data_name = "fsw_data_pv_recent",
  output_plot_filename = "Plots/subgroups/recent_pv_subgroup.png"
)

# Call the function for "ever exposed to physical violence"
process_and_plot(
  data = fsw_data_pv_ever,
  data_name = "fsw_data_pv_ever",
  output_plot_filename = "Plots/subgroups/ever_pv_subgroup.png"
)

# Call the function for "recently exposed to sexual violence"
process_and_plot(
  data = fsw_data_sv_recent,
  data_name = "fsw_data_sv_recent",
  output_plot_filename = "Plots/subgroups/recent_sv_subgroup.png"
)

# Call the function for "ever exposed to sexual violence"
process_and_plot(
  data = fsw_data_sv_ever,
  data_name = "fsw_data_sv_ever",
  output_plot_filename = "Plots/subgroups/ever_sv_subgroup.png"
)

# Call the function for "recently exposed to physical and/or sexual violence"
process_and_plot(
  data = fsw_data_psv_recent,
  data_name = "fsw_data_psv_recent",
  output_plot_filename = "Plots/subgroups/recent_psv_subgroup.png"
)

# Call the function for "ever exposed to physical and/or sexual violence"
process_and_plot(
  data = fsw_data_psv_ever,
  data_name = "fsw_data_psv_ever",
  output_plot_filename = "Plots/subgroups/ever_psv_subgroup.png"
)

## sensitivity analysis

# constant sampling correlation alternatives
rho_1 <- 0.4
rho_2 <- 0.8

perform_analysis_rho1 <- function(df, analysis, exposure, violence_type) {
  # Filter the dataframe
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  filtered_df <- filtered_df %>% filter(!is.na(filtered_df[[var_names[[analysis]]$est]]))
  
  # Create study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)
  
  # Create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho_1,
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
  
  # Ensure the folder exists
  dir.create("Plots/sensitivity_analysis", recursive = TRUE, showWarnings = FALSE)
  
  # Append the suffix to the filename
  filename <- paste0("Plots/sensitivity_analysis/", basename(plot_filenames[[violence_type]][[analysis]][[exposure]]), "_rho1.png")
  
  png(filename = filename, width = 45, height = 22, units = "cm", res = 600)
  
  forest(result2,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols[[exposure]], 
         leftlabs = leftlabs[[exposure]],
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

# Loop over each type of violence, analysis, and exposure to perform the analysis and create forest plots
for (violence_type in names(dataframes)) {
  for (analysis in analyses) {
    for (exposure in exposures) {
      df <- get(dataframes[[violence_type]][[exposure]])
      perform_analysis_rho1(df, analysis, exposure, violence_type)
    }
  }
}

perform_analysis_rho2 <- function(df, analysis, exposure, violence_type) {
  # Filter the dataframe
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  filtered_df <- filtered_df %>% filter(!is.na(filtered_df[[var_names[[analysis]]$est]]))
  
  # Create study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)
  
  # Create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho_2,
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
  
  # Ensure the folder exists
  dir.create("Plots/sensitivity_analysis", recursive = TRUE, showWarnings = FALSE)
  
  # Append the suffix to the filename
  filename <- paste0("Plots/sensitivity_analysis/", basename(plot_filenames[[violence_type]][[analysis]][[exposure]]), "_rho2.png")
  
  png(filename = filename, width = 45, height = 22, units = "cm", res = 600)
  
  forest(result2,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols[[exposure]], 
         leftlabs = leftlabs[[exposure]],
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

# Loop over each type of violence, analysis, and exposure to perform the analysis and create forest plots
for (violence_type in names(dataframes)) {
  for (analysis in analyses) {
    for (exposure in exposures) {
      df <- get(dataframes[[violence_type]][[exposure]])
      perform_analysis_rho2(df, analysis, exposure, violence_type)
    }
  }
}
