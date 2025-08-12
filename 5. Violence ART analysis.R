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

# constant sampling correlation
rho <- 0.6

# models
analyses <- c("best", "unadj", "adj")

# variable names for each analysis
var_names <- list(
  unadj = list(est = "unadj_est_ln", var = "unadj_var_ln", lower = "un_lower_ln", upper = "un_upper_ln"),
  adj = list(est = "adj_est_ln", var = "adj_var_ln", lower = "adj_lower_ln", upper = "adj_upper_ln"),
  best = list(est = "effect_best_ln", var = "effect_best_var_ln", lower = "effect_best_lower_ln", upper = "effect_best_upper_ln")
)

# function for recent violence and art use
perform_analysis_recent_uptake <- function(df, analysis) {
 
  # filter
  filtered_df <- df %>% filter(exposure_tf_bin == "Recent", !is.na(.[[var_names[[analysis]]$est]]))

  # study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)

  # covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # multilevel random effects model using `rma.mv` from metafor
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
    
  # filename
  filename <- paste0("Plots/art use/art uptake/recent_", analysis, "_art_uptake.png")
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

# loop to perform the analysis
for (analysis in analyses) {
  perform_analysis_recent_uptake(fsw_data_art_uptake, analysis)
}

# function for lifetime violence and art uptake
perform_analysis_ever_uptake <- function(df, analysis) {
 
  # filter 
  filtered_df <- df %>% filter(exposure_tf_bin == "Ever", !is.na(.[[var_names[[analysis]]$est]]))
  
  # study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)

  # covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # multilevel random effects model using `rma.mv` from metafor
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
  
  # filename
  filename <- paste0("Plots/art use/art uptake/ever_", analysis, "_art_uptake.png")
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

# loop to perform the analysis
for (analysis in analyses) {
  perform_analysis(fsw_data_art_uptake, analysis)
}

## recent violence and art adherence

# function for recent violence and art adherence
perform_analysis_recent_adherence <- function(df, analysis) {
 
  # filter
  filtered_df <- df %>% filter(exposure_tf_bin == "Recent", !is.na(.[[var_names[[analysis]]$est]]))

  # study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)

  # covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # multilevel random effects model using `rma.mv` from metafor
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
    
  # filename
  filename <- paste0("Plots/art use/art adherence/recent_", analysis, "_art_adherence.png")
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

# loop to perform the analysis
for (analysis in analyses) {
  perform_analysis_recent_adherence(fsw_data_art_adherence, analysis)
}

## subgroup analyses

# recent violence and art adherence
filtered_df <- fsw_data_art_adherence %>%
  filter(exposure_tf_bin == "Recent")

# function for "recent expoure to violence for ART adherence"
process_and_plot(
  data = filtered_df,
  data_name = "filtered_df",
  output_plot_filename = "Plots/subgroups/recent_art_adherence_subgroup.png"
)

# recent violence and art use
filtered_df <- fsw_data_art_uptake %>%
  filter(exposure_tf_bin == "Recent")

# function for "recent expoure to violence for ART use"
process_and_plot(
  data = filtered_df,
  data_name = "filtered_df",
  output_plot_filename = "Plots/subgroups/recent_art_use_subgroup.png"
)

# lifetime violence and art adherence
filtered_df <- fsw_data_art_adherence %>%
  filter(exposure_tf_bin == "Ever")

# function for "ever expoure to violence for ART use"
process_and_plot(
  data = filtered_df,
  data_name = "filtered_df",
  output_plot_filename = "Plots/subgroups/ever_art_adherence_subgroup.png"
)

# recent violence and art use
filtered_df <- fsw_data_art_uptake %>%
  filter(exposure_tf_bin == "Ever")

# function for "ever expoure to violence for ART use"
process_and_plot(
  data = filtered_df,
  data_name = "filtered_df",
  output_plot_filename = "Plots/subgroups/ever_art_use_subgroup.png"
)

## sensitivity analysis

# run for recent and ever violence, art adherence and rho = 0.4
for (exposure in c("Recent", "Ever")) {
  for (analysis in analyses) {
    perform_all_violence_analysis_rho1(fsw_data_art_adherence, analysis, exposure)
  }
}

# run for recent and ever violence, art adherence and rho = 0.8
for (exposure in c("Recent", "Ever")) {
  for (analysis in analyses) {
    perform_all_violence_analysis_rho2(fsw_data_art_adherence, analysis, exposure)
  }
}

# run for recent and ever violence, art use and rho = 0.4
for (exposure in c("Recent", "Ever")) {
  for (analysis in analyses) {
    perform_all_violence_analysis_rho1(fsw_data_art_uptake, analysis, exposure)
  }
}

# run for recent and ever violence, art use and rho = 0.8
for (exposure in c("Recent", "Ever")) {
  for (analysis in analyses) {
    perform_all_violence_analysis_rho2(fsw_data_art_uptake, analysis, exposure)
  }
}