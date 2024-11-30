#  load packages
pacman::p_load("tidyverse")

# Cleaning functions ------------------------------------------------------

# convert to numeric
convert_to_numeric <- function(df) {
  df <- transform(df, 
                  unadj_est = as.numeric(unadj_est), 
                  un_lower = as.numeric(un_lower),
                  un_upper = as.numeric(un_upper),
                  adj_est = as.numeric(adj_est), 
                  adj_lower = as.numeric(adj_lower),
                  adj_upper = as.numeric(adj_upper),
                  effect_best = as.numeric(effect_best), 
                  effect_best_lower = as.numeric(effect_best_lower),
                  effect_best_upper = as.numeric(effect_best_upper))
  return(df)
}

# log transform
log_transform <- function(df) {
  df <- transform(df, 
                  unadj_est_ln = log(unadj_est),
                  un_lower_ln = log(un_lower),
                  un_upper_ln = log(un_upper),
                  adj_est_ln = log(adj_est),
                  adj_lower_ln = log(adj_lower),
                  adj_upper_ln = log(adj_upper),
                  effect_best_ln = log(effect_best),
                  effect_best_lower_ln = log(effect_best_lower),
                  effect_best_upper_ln = log(effect_best_upper))
  return(df)
}

# log transform variance ORs
log_transform_var <- function(df) {
  df <- transform(df, 
                  unadj_var_ln = ((un_upper_ln - un_lower_ln) / (2 * 1.96))^2,
                  adj_var_ln = ((adj_upper_ln - adj_lower_ln) / (2 * 1.96))^2,
                  effect_best_var_ln = ((effect_best_upper_ln - effect_best_lower_ln) / (2 * 1.96))^2)
  return(df)
}

# concatenate unadjusted estimates
format_violence_data_unadj <- function(df) {
  df$effect_unadj_str <- sprintf("%.2f", df$unadj_est)
  df$unaj_lower_str <- sprintf("%.2f", df$un_lower)
  df$unaj_upper_str <- sprintf("%.2f", df$un_upper)
  df$unadj_or_95 <- paste0(df$unaj_lower_str, "–", df$unaj_upper_str)
  df$unadj_or_95_2 <- paste0(df$effect_unadj_str, " (", df$unadj_or_95, ")")
  return(df)
}

# concatenate adjusted estimates
format_violence_data_adj <- function(df) {
  df$effect_adj_str <- sprintf("%.2f", df$adj_est)
  df$adj_lower_str <- sprintf("%.2f", df$adj_lower)
  df$adj_upper_str <- sprintf("%.2f", df$adj_lower)
  df$adj_or_95 <- paste0(df$adj_lower_str, "–", df$adj_upper_str)
  df$adj_or_95_2 <- paste0(df$effect_adj__str, " (", df$adj_or_95, ")")
  return(df)
}

# concatenate best estimates
format_violence_data <- function(df) {
  df$effect_best_str <- sprintf("%.2f", df$effect_best)
  df$lower_str <- sprintf("%.2f", df$effect_best_lower)
  df$upper_str <- sprintf("%.2f", df$effect_best_upper)
  df$or_95 <- paste0(df$lower_str, "–", df$upper_str)
  df$or_95_2 <- paste0(df$effect_best_str, " (", df$or_95, ")")
  return(df)
}

# Overall HIV infection analysis using most exposed variable functions ------------------------------------------------------

# overall unadjusted meta analyses
meta_analysis_unadj <- function(df, df_name) {
  # Filter the dataframe for studies where the outcome is "HIV prevalence"
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  print(paste("Rows after outcome filter:", nrow(filtered_df)))
  
  # Filter where unadjusted effect is missing
  filtered_df <- filtered_df %>% filter(una_effect != "NR")
  print(paste("Rows after una_effect filter:", nrow(filtered_df)))
  
  # Ensure numeric columns for meta-analysis
  filtered_df <- filtered_df %>% 
    mutate(across(c(unadj_est_ln, un_lower_ln, un_upper_ln), as.numeric))
  
  # Remove rows with missing effect estimates
  filtered_df <- filtered_df %>% filter(!is.na(unadj_est_ln) & !is.na(un_lower_ln) & !is.na(un_upper_ln))
  print(paste("Rows after removing NAs:", nrow(filtered_df)))
  
  # Meta-analysis
  if (nrow(filtered_df) == 0) {
    stop("No studies available for meta-analysis.")
  }
  
  result <- metagen(TE = unadj_est_ln,
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
  
  # Save the meta-analysis result using a unique name
  assign(paste0(df_name, "_meta"), result, envir = .GlobalEnv)
  print(summary(result))
  
  # Generate forest plot
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_unadj.png")
  png(filename = filename, width = 50, height = 18, units = "cm", res = 600)
  
  # Check exposure_time_frame
  if (any(filtered_df$exposure_time_frame == "Ever")) {
    leftcols <- c("study", "exposure_definition_short")
    leftlabs <- c("Study", "Exposure definition")
  } else {
    leftcols <- c("study", "exposure_definition_short", "exposure_time_frame")
    leftlabs <- c("Study", "Exposure definition", "Exposure time frame")
  }
  
  forest(result,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols, 
         leftlabs = leftlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  dev.off()
}

# overall adjusted meta analyses

meta_analysis_adj <- function(df, df_name) {
  # Filter the dataframe for studies where the outcome is "HIV prevalence"
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  print(paste("Rows after outcome filter:", nrow(filtered_df)))
  
  # Filter where adjusted effect is missing
  filtered_df <- filtered_df %>% filter(adj_effect != "NR")
  print(paste("Rows after una_effect filter:", nrow(filtered_df)))
  
  # Ensure numeric columns for meta-analysis
  filtered_df <- filtered_df %>% 
    mutate(across(c(adj_est_ln, adj_lower_ln, adj_upper_ln), as.numeric))
  
  # Remove rows with missing effect estimates
  filtered_df <- filtered_df %>% filter(!is.na(adj_est_ln) & !is.na(adj_lower_ln) & !is.na(adj_upper_ln))
  print(paste("Rows after removing NAs:", nrow(filtered_df)))
  
  # Meta-analysis
  if (nrow(filtered_df) == 0) {
    stop("No studies available for meta-analysis.")
  }
  
  result <- metagen(TE = adj_est_ln,
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
  
  # Save the meta-analysis result using a unique name
  assign(paste0(df_name, "_meta"), result, envir = .GlobalEnv)
  print(summary(result))
  
  # Generate forest plot
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_adj.png")
  png(filename = filename, width = 50, height = 18, units = "cm", res = 600)
  
  # Check exposure_time_frame
  if (any(filtered_df$exposure_time_frame == "Ever")) {
    leftcols <- c("study", "exposure_definition_short")
    leftlabs <- c("Study", "Exposure definition")
  } else {
    leftcols <- c("study", "exposure_definition_short", "exposure_time_frame")
    leftlabs <- c("Study", "Exposure definition", "Exposure time frame")
  }
  
  forest(result,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols, 
         leftlabs = leftlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  dev.off()
}

# overall combined meta analyses

meta_analysis_best <- function(df, df_name) {
  # Filter the dataframe for studies where the outcome is "HIV prevalence"
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  print(paste("Rows after outcome filter:", nrow(filtered_df)))
  
  # Ensure numeric columns for meta-analysis
  filtered_df <- filtered_df %>% 
    mutate(across(c(effect_best_ln, effect_best_lower_ln, effect_best_upper_ln), as.numeric))
  
  # Remove rows with missing effect estimates
  filtered_df <- filtered_df %>% filter(!is.na(effect_best_ln) & !is.na(effect_best_lower_ln) & !is.na(effect_best_upper_ln))
  print(paste("Rows after removing NAs:", nrow(filtered_df)))
  
  # Meta-analysis
  if (nrow(filtered_df) == 0) {
    stop("No studies available for meta-analysis.")
  }
  
  result <- metagen(TE = effect_best_ln,
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
  
  # Save the meta-analysis result using a unique name
  assign(paste0(df_name, "_meta"), result, envir = .GlobalEnv)
  print(summary(result))
  
  # Generate forest plot
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_best.png")
  png(filename = filename, width = 50, height = 18, units = "cm", res = 600)
  
  # Check exposure_time_frame
  if (any(filtered_df$exposure_time_frame == "Ever")) {
    leftcols <- c("study", "exposure_definition_short")
    leftlabs <- c("Study", "Exposure definition")
  } else {
    leftcols <- c("study", "exposure_definition_short", "exposure_time_frame")
    leftlabs <- c("Study", "Exposure definition", "Exposure time frame")
  }
  
  forest(result,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols, 
         leftlabs = leftlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  dev.off()
}

# Overall HIV infection analysis multi-level functions ------------------------------------------------------
mlm_analysis_unadj <- function(df, df_name) {
  # Filter the dataframe for studies where the outcome is "HIV prevalence"
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  print(paste("Rows after outcome filter:", nrow(filtered_df)))
  
  # Filter where unadjusted effect is missing
  filtered_df <- filtered_df %>% filter(una_effect != "NR")
  print(paste("Rows after una_effect filter:", nrow(filtered_df)))
  
  # Ensure numeric columns for meta-analysis
  filtered_df <- filtered_df %>% 
    mutate(across(c(unadj_est_ln, adj_var_ln), as.numeric))
  
  result <- rma.mv(
  yi = unadj_est_ln,                   # Log-transformed odds ratios
  V = unadj_var_ln,                    # Variance of log-transformed odds ratios
  slab = study,                        # Study labels
  data = filtered_df,           # Data
  random = ~ 1 | study/study_num,      # Random effects structure
  test = "t",                          # Use t-test for significance
  method = "REML"                      # Restricted maximum likelihood
)

summary(result)
exp(coef(result))

# Save the meta-analysis result using a unique name
#assign(paste0(df_name, "_meta"), result, envir = .GlobalEnv)
print(summary(result))
print(exp(coef(result)))

}

# Overall HIV infection analysis multi-level functions ------------------------------------------------------
mlm_analysis_adj <- function(df, df_name) {
  # Filter the dataframe for studies where the outcome is "HIV prevalence"
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  print(paste("Rows after outcome filter:", nrow(filtered_df)))
  
  # Filter where unadjusted effect is missing
  filtered_df <- filtered_df %>% filter(una_effect != "NR")
  print(paste("Rows after una_effect filter:", nrow(filtered_df)))
  
  # Ensure numeric columns for meta-analysis
  filtered_df <- filtered_df %>% 
    mutate(across(c(adj_est_ln, adj_var_ln), as.numeric))
  
  result <- rma.mv(
    yi = adj_est_ln,                   # Log-transformed odds ratios
    V = adj_var_ln,                    # Variance of log-transformed odds ratios
    slab = study,                        # Study labels
    data = filtered_df,           # Data
    random = ~ 1 | study/study_num,      # Random effects structure
    test = "t",                          # Use t-test for significance
    method = "REML"                      # Restricted maximum likelihood
  )
  
  summary(result)
  exp(coef(result))
  
  # Save the meta-analysis result using a unique name
  #assign(paste0(df_name, "_meta"), result, envir = .GlobalEnv)
  print(summary(result))
  print(exp(coef(result)))
  
}


# Overall HIV infection analysis using most exposed variable functions ------------------------------------------------------


# forest plot for all adjusted estimates

all_adjusted_forest <- function(df, df_name) {
  
  # filter the dataframe
  filtered_df <- df %>% 
    filter(outcome == "HIV prevalence") %>% 
    filter(adj == "Adjusted")
  
  # Meta-analysis
  all_hiv_adj <- metagen(TE = adj_est_ln,
                            lower = adj_lower_ln,
                            upper = adj_upper_ln,
                            studlab = study,
                            data = filtered_df,
                            sm = "OR",
                            method.tau = "DL",
                            comb.fixed = FALSE,
                            comb.random = FALSE, 
                            backtransf = TRUE,
                            text.random = "Overall")
  
  # Print summary
  print(summary(all_hiv_adj)) 
  
  # Create filename using dataframe name
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_adj_plot.png")
  
  # Create PNG for forest plot
  png(filename = filename, 
      width = 45, height = 28, units = "cm", res = 500)
  
  # Generate forest plot
  forest.meta(all_hiv_adj, 
              sortvar = TE,
              xlim = c(0.2, 4),             
              leftcols = c("study", "exposure_definition_short", "exposed_perc_string", "exposure_time_frame"), 
              leftlabs = c("Study", "Exposure definition", "Exposed (%)", "Time frame"),
              rightcols = c("adj_or_95_2" , "who_region", "perpetrator"), 
              rightlabs = c("OR (95% CI)","WHO region", "Perpetrator"), 
              pooled.totals = FALSE,
              xintercept = 1,
              addrow.overall = TRUE,
              test.subgroup = FALSE,
              overall.hetstat = FALSE,
              overall = FALSE,
              labeltext = TRUE,
              col.subgroup = "black",
              print.subgroup.name = FALSE)
  
  # Close the PNG device
  dev.off()
}

# forest plot for all unadjusted estimates

all_unadjusted_forest <- function(df, df_name) {
  
  # filter the dataframe
  filtered_df <- df %>% 
    filter(outcome == "HIV prevalence") %>% 
    filter(adj == "Unadjusted")
  
  # Meta-analysis
  all_hiv_unadj <- metagen(TE = unadj_est_ln,
                         lower = un_lower_ln,
                         upper = un_upper_ln,
                         studlab = study,
                         data = filtered_df,
                         sm = "OR",
                         method.tau = "DL",
                         comb.fixed = FALSE,
                         comb.random = FALSE, 
                         backtransf = TRUE,
                         text.random = "Overall")
  
  # Print summary
  print(summary(all_hiv_unadj)) 
  
  # Create filename using dataframe name
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_unadj_plot.png")
  
  # Create PNG for forest plot
  png(filename = filename, 
      width = 45, height = 28, units = "cm", res = 500)
  
  # Generate forest plot
  forest.meta(all_hiv_unadj, 
              sortvar = TE,
              xlim = c(0.2, 4),             
              leftcols = c("study", "exposure_definition_short", "exposed_perc_string", "exposure_time_frame"), 
              leftlabs = c("Study", "Exposure definition", "Exposed (%)", "Time frame"),
              rightcols = c("unadj_or_95_2" , "who_region", "perpetrator"), 
              rightlabs = c("OR (95% CI)","WHO region", "Perpetrator"), 
              pooled.totals = FALSE,
              xintercept = 1,
              addrow.overall = TRUE,
              test.subgroup = FALSE,
              overall.hetstat = FALSE,
              overall = FALSE,
              labeltext = TRUE,
              col.subgroup = "black",
              print.subgroup.name = FALSE)
  
  # Close the PNG device
  dev.off()
}

# forest plot for all best estimates

all_best_forest <- function(df, df_name) {
  
  # filter the dataframe
  filtered_df <- df %>% 
    filter(outcome == "HIV prevalence")
  
  # Meta-analysis
  all_hiv_best <- metagen(TE = effect_best_ln,
                           lower = effect_best_lower_ln,
                           upper = effect_best_upper_ln,
                           studlab = study,
                           data = filtered_df,
                           sm = "OR",
                           method.tau = "DL",
                           comb.fixed = FALSE,
                           comb.random = FALSE, 
                           backtransf = TRUE,
                           text.random = "Overall")
  
  # Print summary
  print(summary(all_hiv_best)) 
  
  # Create filename using dataframe name
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_best_plot.png")
  
  # Create PNG for forest plot
  png(filename = filename, 
      width = 45, height = 28, units = "cm", res = 500)
  
  # Generate forest plot
  forest.meta(all_hiv_best, 
              sortvar = TE,
              xlim = c(0.2, 4),             
              leftcols = c("study", "exposure_definition_short", "exposed_perc_string", "exposure_time_frame"), 
              leftlabs = c("Study", "Exposure definition", "Exposed (%)", "Time frame"),
              rightcols = c("adj" ,"or_95_2" , "who_region", "perpetrator"), 
              rightlabs = c("Model" ,"OR (95% CI)","WHO region", "Perpetrator"), 
              pooled.totals = FALSE,
              xintercept = 1,
              addrow.overall = TRUE,
              test.subgroup = FALSE,
              overall.hetstat = FALSE,
              overall = FALSE,
              labeltext = TRUE,
              col.subgroup = "black",
              print.subgroup.name = FALSE)
  
  # Close the PNG device
  dev.off()
}

# forest plot for ever best

meta_analysis_vs <- function(df, df_name) {
  # Filter the dataframe for studies where the outcome is "HIV prevalence"
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  
  # Perform the meta-analysis
  result <- metagen(TE = effect_best_ln,
                    lower = effect_best_lower_ln,
                    upper = effect_best_upper_ln,
                    studlab = study,
                    data = filtered_df,
                    sm = "OR",
                    method.tau = "REML",
                    comb.fixed = FALSE,
                    comb.random = TRUE, 
                    backtransf = TRUE,
                    text.random = "Overall")
  
  
  assign(paste0(df_name, "_forest"), result, envir = .GlobalEnv)
  print(summary(result))
}



# HIV infection subgroup functions ---------------------------------------------

# LIC subgroup 
meta_analysis_best_ldc <- function(df, df_name, suffix = "_forest") {
  
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  result <- metagen(TE = effect_best_ln,
                    lower = effect_best_lower_ln,
                    upper = effect_best_upper_ln,
                    studlab = study,
                    data = filtered_df,
                    sm = "OR",
                    method.tau = "REML",
                    comb.fixed = FALSE,
                    comb.random = TRUE, 
                    backtransf = TRUE,
                    byvar = ldc_bin,
                    text.random = "Overall") 
  assign(paste0(df_name, suffix), result, envir = .GlobalEnv)
  print(summary(result))  
  
  # Create filename using dataframe name
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_lic_plot_best.png")
  
  # Create PNG for forest plot
  png(filename = filename, 
      width = 20, height = 45, units = "cm", res = 600)
  
  result_forest <- forest.meta(result, 
                               sortvar = study,
                               xlim=c(0.2, 4),             
                               leftcols = c("study", "adj"), 
                               leftlabs = c("Study", "Model"),
                               pooled.totals = F,
                               xintercept=1,
                               addrow.overall = T,
                               test.subgroup = T,
                               overall.hetstat = F,
                               overall = F,
                               labeltext = TRUE,
                               col.subgroup = "black",
                               print.subgroup.name = FALSE) 
  dev.off()
}

# WHO subgroup 
meta_analysis_best_who_sv <- function(df, df_name, suffix = "_forest") {
  
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  result <- metagen(TE = effect_best_ln,
                    lower = effect_best_lower_ln,
                    upper = effect_best_upper_ln,
                    studlab = study,
                    data = filtered_df,
                    sm = "OR",
                    method.tau = "REML",
                    comb.fixed = FALSE,
                    comb.random = TRUE, 
                    backtransf = TRUE,
                    byvar = who_region,
                    text.random = "Overall") 
  assign(paste0(df_name, suffix), result, envir = .GlobalEnv)
  print(summary(result))  
  
  # Create filename using dataframe name
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_who_plot_best.png")
  
  # Create PNG for forest plot
  png(filename = filename, 
      width = 20, height = 45, units = "cm", res = 600)
  
  result_forest <- forest.meta(result, 
                               sortvar = study,
                               xlim=c(0.2, 4),             
                               leftcols = c("study", "adj"), 
                               leftlabs = c("Study", "Model"),
                               pooled.totals = F,
                               xintercept=1,
                               addrow.overall = T,
                               test.subgroup = T,
                               overall.hetstat = F,
                               overall = F,
                               labeltext = TRUE,
                               col.subgroup = "black",
                               print.subgroup.name = FALSE) 
  dev.off()
}

# recruitment method subgroup
meta_analysis_best_recruit <- function(df, df_name, suffix = "_forest") {
  
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  result <- metagen(TE = effect_best_ln,
                    lower = effect_best_lower_ln,
                    upper = effect_best_upper_ln,
                    studlab = study,
                    data = filtered_df,
                    sm = "OR",
                    method.tau = "REML",
                    comb.fixed = FALSE,
                    comb.random = TRUE, 
                    backtransf = TRUE,
                    byvar = recruitment,
                    text.random = "Overall") 
  assign(paste0(df_name, suffix), result, envir = .GlobalEnv)
  print(summary(result))  
  
  # Create filename using dataframe name
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_recruit_plot_best.png")
  
  # Create PNG for forest plot
  png(filename = filename, 
      width = 20, height = 45, units = "cm", res = 600)
  
  result_forest <- forest.meta(result, 
                               sortvar = study,
                               xlim=c(0.2, 4),             
                               leftcols = c("study", "adj"), 
                               leftlabs = c("Study", "Model"),
                               pooled.totals = F,
                               xintercept=1,
                               addrow.overall = T,
                               test.subgroup = T,
                               overall.hetstat = F,
                               overall = F,
                               labeltext = TRUE,
                               col.subgroup = "black",
                               print.subgroup.name = FALSE) 
  dev.off()
}

# study setting subgroup
meta_analysis_best_setting <- function(df, df_name, suffix = "_forest") {
  
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  result <- metagen(TE = effect_best_ln,
                    lower = effect_best_lower_ln,
                    upper = effect_best_upper_ln,
                    studlab = study,
                    data = filtered_df,
                    sm = "OR",
                    method.tau = "REML",
                    comb.fixed = FALSE,
                    comb.random = TRUE, 
                    backtransf = TRUE,
                    byvar = study_setting,
                    text.random = "Overall") 
  assign(paste0(df_name, suffix), result, envir = .GlobalEnv)
  print(summary(result))  
  
  # Create filename using dataframe name
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_setting_plot_best.png")
  
  # Create PNG for forest plot
  png(filename = filename, 
      width = 20, height = 45, units = "cm", res = 600)
  
  result_forest <- forest.meta(result, 
                               sortvar = study,
                               xlim=c(0.2, 4),             
                               leftcols = c("study", "adj"), 
                               leftlabs = c("Study", "Model"),
                               pooled.totals = F,
                               xintercept=1,
                               addrow.overall = T,
                               test.subgroup = T,
                               overall.hetstat = F,
                               overall = F,
                               labeltext = TRUE,
                               col.subgroup = "black",
                               print.subgroup.name = FALSE) 
  dev.off()
}

# perpetrator subgroup
meta_analysis_best_perp <- function(df, df_name, suffix = "_forest") {
  
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  result <- metagen(TE = effect_best_ln,
                    lower = effect_best_lower_ln,
                    upper = effect_best_upper_ln,
                    studlab = study,
                    data = filtered_df,
                    sm = "OR",
                    method.tau = "REML",
                    comb.fixed = FALSE,
                    comb.random = TRUE, 
                    backtransf = TRUE,
                    byvar = perpetrator,
                    text.random = "Overall") 
  assign(paste0(df_name, suffix), result, envir = .GlobalEnv)
  print(summary(result))  
  
  # Create filename using dataframe name
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_perp_plot_best.png")
  
  # Create PNG for forest plot
  png(filename = filename, 
      width = 20, height = 45, units = "cm", res = 600)
  
  result_forest <- forest.meta(result, 
                               sortvar = study,
                               xlim=c(0.2, 4),             
                               leftcols = c("study", "adj"), 
                               leftlabs = c("Study", "Model"),
                               pooled.totals = F,
                               xintercept=1,
                               addrow.overall = T,
                               test.subgroup = T,
                               overall.hetstat = F,
                               overall = F,
                               labeltext = TRUE,
                               col.subgroup = "black",
                               print.subgroup.name = FALSE) 
  dev.off()
}

# all perpetrator subgroup
meta_analysis_best_perp_all <- function(df, df_name, suffix = "_forest") {
  
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  result <- metagen(TE = effect_best_ln,
                    lower = effect_best_lower_ln,
                    upper = effect_best_upper_ln,
                    studlab = study,
                    data = filtered_df,
                    sm = "OR",
                    method.tau = "REML",
                    comb.fixed = FALSE,
                    comb.random = TRUE, 
                    backtransf = TRUE,
                    byvar = perpetrator,
                    text.random = "Overall") 
  assign(paste0(df_name, suffix), result, envir = .GlobalEnv)
  print(summary(result))  
  
  # Create filename using dataframe name
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_perp_all_plot_best.png")
  
  # Create PNG for forest plot
  png(filename = filename, 
      width = 20, height = 45, units = "cm", res = 600)
  
  result_forest <- forest.meta(result, 
                               sortvar = study,
                               xlim=c(0.2, 4),             
                               leftcols = c("study", "adj"), 
                               leftlabs = c("Study", "Model"),
                               pooled.totals = F,
                               xintercept=1,
                               addrow.overall = T,
                               test.subgroup = T,
                               overall.hetstat = F,
                               overall = F,
                               labeltext = TRUE,
                               col.subgroup = "black",
                               print.subgroup.name = FALSE) 
  dev.off()
}

# pre-2016 subgroup
meta_analysis_best_pre2016 <- function(df, df_name, suffix = "_forest") {
  
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  result <- metagen(TE = effect_best_ln,
                    lower = effect_best_lower_ln,
                    upper = effect_best_upper_ln,
                    studlab = study,
                    data = filtered_df,
                    sm = "OR",
                    method.tau = "REML",
                    comb.fixed = FALSE,
                    comb.random = TRUE, 
                    backtransf = TRUE,
                    byvar = pre_2016,
                    text.random = "Overall") 
  assign(paste0(df_name, suffix), result, envir = .GlobalEnv)
  print(summary(result))  
  
  # Create filename using dataframe name
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_pre2016_plot_best.png")
  
  # Create PNG for forest plot
  png(filename = filename, 
      width = 20, height = 45, units = "cm", res = 600)
  
  result_forest <- forest.meta(result, 
                               sortvar = study,
                               xlim=c(0.2, 4),             
                               leftcols = c("study", "adj"), 
                               leftlabs = c("Study", "Model"),
                               pooled.totals = F,
                               xintercept=1,
                               addrow.overall = T,
                               test.subgroup = T,
                               overall.hetstat = F,
                               overall = F,
                               labeltext = TRUE,
                               col.subgroup = "black",
                               print.subgroup.name = FALSE) 
  dev.off()
}

# risk of bias subgroup


# subgroup cleaning
process_violence_data <- function(sheet_name) {
  
  data <- read_excel(
    "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Violence subgroup estimates.xlsx", 
    sheet = sheet_name
  )
  
  data <- transform(data, 
                    effect_ln = log(effect),
                    lower_ln = log(lower),
                    upper_ln = log(upper))
  
  data$lower_str <- format(data$lower, nsmall = 2)
  data$upper_str <- format(data$upper, nsmall = 2)
  data$or_95 <- paste0(data$lower_str, "–", data$upper_str)
  data$or_95_2 <- paste0(data$effect, " ", "(", data$or_95, ")")
  
  df_name <- paste0("summary_violence_", gsub(" ", "_", sheet_name))
  assign(df_name, data, envir = .GlobalEnv)
  
  return(df_name)
}

# sexual violence subgroup forest plots
subgroup_plot_forest_sv <- function(df_name) {
  
  data <- get(df_name)
  data_sv <- data[data$name == "Sexual violence", ]
  
  meta_analysis_result <- metagen(TE = effect_ln,
                                  lower = lower_ln,
                                  upper = upper_ln,
                                  studlab = name,
                                  data = data_sv,
                                  sm = "OR",
                                  method.tau = "DL",
                                  comb.fixed = FALSE,
                                  comb.random = FALSE, 
                                  backtransf = TRUE,
                                  byvar = time_frame,
                                  text.random = "Overall")
  
  summary(meta_analysis_result)
  
  plot_filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/plots/", df_name, "_sv_forest.png")
  
  png(filename = plot_filename, width = 28, height = 44, units = "cm", res = 500)
  
  forest.meta(meta_analysis_result, 
              sortvar = name,
              xlim = c(0.2, 4),             
              leftcols = c("category"), 
              leftlabs = c("Category"),
              rightcols = c("or_95_2", "studies"), 
              rightlabs = c("OR (95% CI)", "Studies"), 
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
  
  return(meta_analysis_result)
}

# physical violence subgroup forest plots
subgroup_plot_forest_pv <- function(df_name) {
  
  data <- get(df_name)
  data_sv <- data[data$name == "Physical violence", ]
  
  meta_analysis_result <- metagen(TE = effect_ln,
                                  lower = lower_ln,
                                  upper = upper_ln,
                                  studlab = name,
                                  data = data_sv,
                                  sm = "OR",
                                  method.tau = "DL",
                                  comb.fixed = FALSE,
                                  comb.random = FALSE, 
                                  backtransf = TRUE,
                                  byvar = time_frame,
                                  text.random = "Overall")
  
  summary(meta_analysis_result)
  
  plot_filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/plots/", df_name, "_pv_forest.png")
  
  png(filename = plot_filename, width = 28, height = 44, units = "cm", res = 500)
  
  forest.meta(meta_analysis_result, 
              sortvar = name,
              xlim = c(0.2, 4),             
              leftcols = c("category"), 
              leftlabs = c("Category"),
              rightcols = c("or_95_2", "studies"), 
              rightlabs = c("OR (95% CI)", "Studies"), 
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
  
  return(meta_analysis_result)
}

# physical and/or sexual violence subgroup forest plots
subgroup_plot_forest_psv <- function(df_name) {
  
  data <- get(df_name)
  data_sv <- data[data$name == "Sexual and/or physical violence", ]
  
  meta_analysis_result <- metagen(TE = effect_ln,
                                  lower = lower_ln,
                                  upper = upper_ln,
                                  studlab = name,
                                  data = data_sv,
                                  sm = "OR",
                                  method.tau = "DL",
                                  comb.fixed = FALSE,
                                  comb.random = FALSE, 
                                  backtransf = TRUE,
                                  byvar = time_frame,
                                  text.random = "Overall")
  
  summary(meta_analysis_result)
  
  plot_filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/plots/", df_name, "_psv_forest.png")
  
  png(filename = plot_filename, width = 28, height = 44, units = "cm", res = 500)
  
  forest.meta(meta_analysis_result, 
              sortvar = name,
              xlim = c(0.2, 4),             
              leftcols = c("category"), 
              leftlabs = c("Category"),
              rightcols = c("or_95_2", "studies"), 
              rightlabs = c("OR (95% CI)", "Studies"), 
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
  
  return(meta_analysis_result)
}


# Overall ART and VS analysis functions ------------------------------------------

# unadjusted recent
meta_analysis_art_vs_unadj_rec <- function(df, df_name) {

filtered_df <- df %>% 
    filter(una_effect != "NR") %>% 
    filter(use != "no") %>%
    filter(exposure_tf_bin == "Recent")

# meta-analysis
result <- metagen(TE = unadj_est_ln,
                  lower = un_lower_ln,
                  upper = un_upper_ln,
                  studlab = study,
                  data = filtered_df,
                  sm = "OR",
                  method.tau = "REML",
                  comb.fixed = FALSE,
                  comb.random = TRUE, 
                  backtransf = TRUE) 

assign(paste0(df_name), result, envir = .GlobalEnv)
print(summary(result))  
  
filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_unadj_rec.png")
  
png(filename = filename, width = 50, height = 14, units = "cm",
      res = 600)
  
  # check if exposure_time_frame contains "Ever"
  if (any(filtered_df$exposure_time_frame == "Ever")) {
    leftcols <- c("study", "exposure_definition_short")
    leftlabs <- c("Study", "Exposure definition")
  } else {
    leftcols <- c("study", "exposure_definition_short", "exposure_time_frame")
    leftlabs <- c("Study", "Exposure definition", "Exposure time frame")
  }
  
  forest(result,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols, 
         leftlabs = leftlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  # close 
  dev.off()
  
}


# adjusted recent
meta_analysis_art_vs_adj_rec <- function(df, df_name) {
  
  filtered_df <- df %>% 
    filter(adj_effect != "NR") %>% 
    filter(use != "no") %>%
    filter(exposure_tf_bin == "Recent")
  
  # meta-analysis
  result <- metagen(TE = adj_est_ln,
                    lower = adj_lower_ln,
                    upper = adj_upper_ln,
                    studlab = study,
                    data = filtered_df,
                    sm = "OR",
                    method.tau = "REML",
                    comb.fixed = FALSE,
                    comb.random = TRUE, 
                    backtransf = TRUE) 
  
  assign(paste0(df_name), result, envir = .GlobalEnv)
  print(summary(result))  
  
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_adj_rec.png")
  
  png(filename = filename, width = 50, height = 14, units = "cm",
      res = 600)
  
  # check if exposure_time_frame contains "Ever"
  if (any(filtered_df$exposure_time_frame == "Ever")) {
    leftcols <- c("study", "exposure_definition_short")
    leftlabs <- c("Study", "Exposure definition")
  } else {
    leftcols <- c("study", "exposure_definition_short", "exposure_time_frame")
    leftlabs <- c("Study", "Exposure definition", "Exposure time frame")
  }
  
  forest(result,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols, 
         leftlabs = leftlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  # close 
  dev.off()
  
}


# best recent
meta_analysis_art_vs_best_rec <- function(df, df_name) {
  
  filtered_df <- df %>% 
    filter(use != "no") %>%
    filter(exposure_tf_bin == "Recent")
  
  # meta-analysis
  result <- metagen(TE = effect_best_ln,
                    lower = effect_best_lower_ln,
                    upper = effect_best_upper_ln,
                    studlab = study,
                    data = filtered_df,
                    sm = "OR",
                    method.tau = "REML",
                    comb.fixed = FALSE,
                    comb.random = TRUE, 
                    backtransf = TRUE) 
  
  assign(paste0(df_name), result, envir = .GlobalEnv)
  print(summary(result))  
  
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_best_rec.png")
  
  png(filename = filename, width = 50, height = 14, units = "cm",
      res = 600)
  
  # check if exposure_time_frame contains "Ever"
  if (any(filtered_df$exposure_time_frame == "Ever")) {
    leftcols <- c("study", "exposure_definition_short")
    leftlabs <- c("Study", "Exposure definition")
  } else {
    leftcols <- c("study", "exposure_definition_short", "exposure_time_frame")
    leftlabs <- c("Study", "Exposure definition", "Exposure time frame")
  }
  
  forest(result,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols, 
         leftlabs = leftlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  # close 
  dev.off()
  
}

# unadjusted ever
meta_analysis_art_vs_unadj_ever <- function(df, df_name) {
  
  filtered_df <- df %>% 
    filter(una_effect != "NR") %>% 
    filter(use != "no") %>%
    filter(exposure_tf_bin == "Ever")
  
  # meta-analysis
  result <- metagen(TE = unadj_est_ln,
                    lower = un_lower_ln,
                    upper = un_upper_ln,
                    studlab = study,
                    data = filtered_df,
                    sm = "OR",
                    method.tau = "REML",
                    comb.fixed = FALSE,
                    comb.random = TRUE, 
                    backtransf = TRUE) 
  
  assign(paste0(df_name), result, envir = .GlobalEnv)
  print(summary(result))  
  
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_unadj_ever.png")
  
  png(filename = filename, width = 50, height = 14, units = "cm",
      res = 600)
  
    leftcols <- c("study", "exposure_definition_short")
    leftlabs <- c("Study", "Exposure definition")

  forest(result,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols, 
         leftlabs = leftlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  # close 
  dev.off()
  
}


# adjusted ever
meta_analysis_art_vs_adj_ever <- function(df, df_name) {
  
  filtered_df <- df %>% 
    filter(adj_effect != "NR") %>% 
    filter(use != "no") %>%
    filter(exposure_tf_bin == "Ever")
  
  # meta-analysis
  result <- metagen(TE = adj_est_ln,
                    lower = adj_lower_ln,
                    upper = adj_upper_ln,
                    studlab = study,
                    data = filtered_df,
                    sm = "OR",
                    method.tau = "REML",
                    comb.fixed = FALSE,
                    comb.random = TRUE, 
                    backtransf = TRUE) 
  
  assign(paste0(df_name), result, envir = .GlobalEnv)
  print(summary(result))  
  
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_adj_ever.png")
  
  png(filename = filename, width = 50, height = 14, units = "cm",
      res = 600)
  
  leftcols <- c("study", "exposure_definition_short")
  leftlabs <- c("Study", "Exposure definition")
  
  forest(result,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols, 
         leftlabs = leftlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  # close 
  dev.off()
  
}


# best ever
meta_analysis_art_vs_best_ever <- function(df, df_name) {
  
  filtered_df <- df %>% 
    filter(use != "no") %>%
    filter(exposure_tf_bin == "Ever")
  
  # meta-analysis
  result <- metagen(TE = effect_best_ln,
                    lower = effect_best_lower_ln,
                    upper = effect_best_upper_ln,
                    studlab = study,
                    data = filtered_df,
                    sm = "OR",
                    method.tau = "REML",
                    comb.fixed = FALSE,
                    comb.random = TRUE, 
                    backtransf = TRUE) 
  
  assign(paste0(df_name), result, envir = .GlobalEnv)
  print(summary(result))  
  
  filename <- paste0("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/", df_name, "_best_ever.png")
  
  png(filename = filename, width = 50, height = 14, units = "cm",
      res = 600)
  
  leftcols <- c("study", "exposure_definition_short")
  leftlabs <- c("Study", "Exposure definition")
  
  forest(result,
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols, 
         leftlabs = leftlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  # close 
  dev.off()
  
}