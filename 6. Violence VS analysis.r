# load packages 
pacman::p_load("meta", "metafor", "readxl", "openxlsx", "tidyverse", "kableExtra", "robumeta", "clubSandwich") 

# set working directory
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence")

# settings
settings.meta(CIbracket = "(") 
settings.meta(CIseparator = "-") 

# columns
leftcols_recent <- c("study", "study_num", "effect_num", "exposure_definition_short", "exposure_time_frame", "perpetrator", "country", "outcome_definition")
leftlabs_recent <- c("Study", "Study number", "Effect number", "Exposure definition", "Exposure time frame", "Perpetrator", "Country", "VS cutoff")
leftcols_lifetime <- c("study", "study_num", "effect_num", "exposure_definition_short", "perpetrator", "country", "outcome_definition")
leftlabs_lifetime <- c("Study", "Study number", "Effect number", "Exposure definition", "Perpetrator", "Country", "VS cutoff")
rightcols <- c("effect", "ci")
rightlabs = c("Estimate", "95% CI")

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
  unadj = "Plots/recent_unadj_vs.png",
  adj = "Plots/recent_adj_vs.png",
  best = "Plots/recent_best_vs.png"
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
  perform_analysis(fsw_data_vs, analysis)
}

## lifetime violence data

# Define the corresponding plot filenames for each analysis
plot_filenames <- list(
  unadj = "Plots/ever_unadj_vs.png",
  adj = "Plots/ever_adj_vs.png",
  best = "Plots/ever_best_vs.png"
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
  perform_analysis(fsw_data_vs, analysis)
}


## combined recent and lifetime

# Define the corresponding plot filenames for each analysis
plot_filenames <- list(
  unadj = "Plots/all_unadj_vs.png",
  adj = "Plots/all_adj_vs.png",
  best = "Plots/all_best_vs.png"
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
  perform_analysis(fsw_data_vs, analysis)
}

## no outlier

# Define the corresponding plot filenames for each analysis
plot_filenames <- list(
  unadj = "Plots/all_unadj_vs_no_outlier.png",
  adj = "Plots/all_adj_vs_no_outlier.png",
  best = "Plots/all_best_vs_no_outlier.png"
)

# Function to perform the analysis and create forest plots
perform_analysis <- function(df, analysis) {
 
  # Filter the dataframe
  filtered_df <- df %>% filter(use == "yes", author != "Wilson", !is.na(.[[var_names[[analysis]]$est]]))
  
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
                   control = list(optimizer = "nlminb", start = list(tau2 = 0.1)))
  
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
  perform_analysis(fsw_data_vs, analysis)
}
