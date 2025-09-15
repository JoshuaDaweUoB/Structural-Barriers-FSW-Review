# load packages 
pacman::p_load("meta", "metafor", "readxl", "tidyverse", "kableExtra", "robumeta", "clubSandwich") 

# set working directory
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence")

# settings
settings.meta(CIbracket = "(") 
settings.meta(CIseparator = "-") 

# columns
leftcols_recent <- c("study", "study_num", "effect_num", "exposure_definition_short", "exposure_time_frame", "outcome_time_frame", "perpetrator", "country")
leftlabs_recent <- c("Study", "Study number", "Effect number", "Exposure definition", "Exposure time frame", "Testing timeframe", "Perpetrator", "Country")
leftcols_lifetime <- c("study", "study_num", "effect_num", "exposure_definition_short",  "outcome_time_frame", "perpetrator", "country")
leftlabs_lifetime <- c("Study", "Study number", "Effect number", "Exposure definition", "Testing timeframe", "Perpetrator", "Country")
rightcols <- c("effect", "ci")
rightlabs = c("Estimate", "95% CI")

# lists for loops and functions
analyses <- c("unadj", "adj", "best")
exposures <- c("recent", "ever")
violence_types <- c("Physical", "Sexual", "Physical or sexual")

## modelling

# constant sampling correlation
rho <- 0.6

# list of model types
analyses <- c("unadj", "adj", "best")

# variable names for each analysis
var_names <- list(
  unadj = list(est = "unadj_est_ln", var = "unadj_var_ln", lower = "un_lower_ln", upper = "un_upper_ln"),
  adj = list(est = "adj_est_ln", var = "adj_var_ln", lower = "adj_lower_ln", upper = "adj_upper_ln"),
  best = list(est = "effect_best_ln", var = "effect_best_var_ln", lower = "effect_best_lower_ln", upper = "effect_best_upper_ln")
)

# labels for plots
analysis_labels <- c(
  best = "Combined",
  unadj = "Unadjusted",
  adj = "Adjusted"
)
exposure_labels <- c(
  recent = "Recent",
  ever = "Ever"
)

violence_type_labels <- c(
  physical = "Physical",
  sexual = "Sexual",
  physical_sexual = "Physical or sexual",
  other = "Other violence"
)

## recent violence data 

# function for recent violence
perform_analysis_recent <- function(df, analysis) {
 
  # filter
  filtered_df <- df %>% filter(exposure_tf_bin == "Recent", !is.na(.[[var_names[[analysis]]$est]]))
  
  # study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)

  # covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # multilevel random effects
  result <- rma.mv(filtered_df[[var_names[[analysis]]$est]], 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE)       
  
  print(paste("RECENT", analysis, "analysis:"))
  print(result)
  print(exp(coef(result)))
  
  # analysis
  result2 <- metagen(
    TE = filtered_df[[var_names[[analysis]]$est]],
    lower = filtered_df[[var_names[[analysis]]$lower]],
    upper = filtered_df[[var_names[[analysis]]$upper]],
    studlab = filtered_df$study,
    data = filtered_df,
    sm = "OR",
    method.tau = "REML",
    common = FALSE,
    random = TRUE, 
    backtransf = TRUE,
    byvar = filtered_df$outcome_tf_bin,
    text.random = "Overall",
    print.byvar = FALSE
  )
  
  print(summary(result2))
  
  result2$TE.random <- result$b
  result2$lower.random <- result$ci.lb
  result2$upper.random <- result$ci.ub
  
  filename <- paste0("Plots/testing/", analysis, "_recent_testing.png")
  png(filename = filename, width = 55, height = 14, units = "cm", res = 600)
  
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
         col.subgroup = "black",
         byvar = filtered_df$outcome_tf_bin
  )
  
  dev.off()

  # eggers test
  eggers <- metabias(result2, method.bias = "linreg")
  eggers_p <- if (!is.null(eggers$p.value)) eggers$p.value else NA
  eggers_p_str <- if (!is.na(eggers_p)) sprintf("p = %.3f", eggers_p) else ""

  # funnel plot
  funnel_label <- paste0(analysis_labels[[analysis]], " - ", exposure_labels[["recent"]])
  sanitized_label <- gsub("[/\\?<>\\:*|\"]", "-", funnel_label)
  funnel_filename <- paste0("Plots/testing/funnel plots/", sanitized_label, ".png")

  png(filename = funnel_filename, width = 15, height = 15, units = "cm", res = 300)
  funnel(result2, main = paste0(funnel_label, "\nEgger's test ", eggers_p_str))
  dev.off()
}

# loop over each analysis
for (analysis in analyses) {
  perform_analysis_recent(fsw_data_test, analysis)
}

## lifetime violence data

# function for lifetime violence
perform_analysis_ever <- function(df, analysis) {
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
  
  print(paste("EVER", analysis, "analysis:"))
  print(result)
  print(exp(coef(result)))
  
  result2 <- metagen(
    TE = filtered_df[[var_names[[analysis]]$est]],
    lower = filtered_df[[var_names[[analysis]]$lower]],
    upper = filtered_df[[var_names[[analysis]]$upper]],
    studlab = filtered_df$study,
    data = filtered_df,
    sm = "OR",
    method.tau = "REML",
    common = FALSE,
    random = TRUE, 
    backtransf = TRUE,
    byvar = filtered_df$outcome_tf_bin, 
    text.random = "Overall",
    print.byvar = FALSE    
  )
  
  print(summary(result2))
  
  result2$TE.random <- result$b
  result2$lower.random <- result$ci.lb
  result2$upper.random <- result$ci.ub
  
  filename <- paste0("Plots/testing/", analysis, "_ever_testing.png")
  png(filename = filename, width = 45, height = 18, units = "cm", res = 600)
  
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
         col.subgroup = "black",
         byvar = filtered_df$outcome_tf_bin, # <-- add subgroup here
         print.byvar = FALSE                 # <-- hide variable name in subgroup header
  )
  
  dev.off()

  # eggers test
  eggers <- metabias(result2, method.bias = "linreg")
  eggers_p <- if (!is.null(eggers$p.value)) eggers$p.value else NA
  eggers_p_str <- if (!is.na(eggers_p)) sprintf("p = %.3f", eggers_p) else ""

  # funnel plot
  funnel_label <- paste0(analysis_labels[[analysis]], " - ", exposure_labels[["ever"]])
  sanitized_label <- gsub("[/\\?<>\\:*|\"]", "-", funnel_label)
  funnel_filename <- paste0("Plots/testing/funnel plots/", sanitized_label, ".png")

  png(filename = funnel_filename, width = 15, height = 15, units = "cm", res = 300)
  funnel(result2, main = paste0(funnel_label, "\nEgger's test ", eggers_p_str))
  dev.off()
}

# loop over each analysis
for (analysis in analyses) {
  perform_analysis_ever(fsw_data_test, analysis)
}

# function for combining funnel plots
get_testing_funnel_path <- function(analysis, exposure) {
  paste0(
    "Plots/testing/funnel plots/",
    analysis_labels[[analysis]], " - ", exposure_labels[[exposure]], ".png"
  )
}

# funnel plots list
all_funnel_imgs <- vector("list", length = length(exposures) * length(analyses))
idx <- 1
for (i in seq_along(exposures)) {
  for (j in seq_along(analyses)) {
    file <- get_testing_funnel_path(analyses[j], exposures[i])
    if (file.exists(file)) {
      all_funnel_imgs[[idx]] <- rasterGrob(readPNG(file), interpolate = TRUE)
    } else {
      all_funnel_imgs[[idx]] <- nullGrob()
    }
    idx <- idx + 1
  }
}

# combine figures 
combined_filename <- "Plots/testing/funnel plots/testing_funnel_grid.png"
png(combined_filename, width = 1800, height = 1200, res = 150)
grid.arrange(
  grobs = all_funnel_imgs,
  nrow = length(exposures),
  ncol = length(analyses),
  top = "Testing: Funnel plots"
)
dev.off()

## subgroup analysis

filtered_df <- fsw_data_test %>%
  filter(exposure_tf_bin == "Recent")

# function for "recent violence"
process_and_plot(
  data = filtered_df,
  data_name = "filtered_df",
  output_plot_filename = "Plots/subgroups/recent_test_subgroup.png"
)

filtered_df <- fsw_data_test %>%
  filter(exposure_tf_bin == "Ever")

# function for "ever violence"
process_and_plot(
  data = filtered_df,
  data_name = "filtered_df",
  output_plot_filename = "Plots/subgroups/ever_test_subgroup.png"
)

## sensitivity analysis

# run for recent and ever violence and rho = 0.4
for (exposure in c("Recent", "Ever")) {
  for (analysis in analyses) {
    perform_all_violence_analysis_rho1(fsw_data_test, analysis, exposure)
  }
}

# run for recent and ever violence and rho = 0.8
for (exposure in c("Recent", "Ever")) {
  for (analysis in analyses) {
    perform_all_violence_analysis_rho2(fsw_data_test, analysis, exposure)
  }
}
