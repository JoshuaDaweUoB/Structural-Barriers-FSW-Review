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

# dataframes for each exposure and type of violence
dataframes <- list(
  physical = list(recent = "fsw_data_pv_recent", ever = "fsw_data_pv_ever"),
  sexual = list(recent = "fsw_data_sv_recent", ever = "fsw_data_sv_ever"),
  physical_sexual = list(recent = "fsw_data_psv_recent", ever = "fsw_data_psv_ever")
)

# filenames for each analysis, exposure, and type of violence
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

# left columns and labels for each exposure
leftcols <- list(
  recent = leftcols_recent,
  ever = leftcols_lifetime
)
leftlabs <- list(
  recent = leftlabs_recent,
  ever = leftlabs_lifetime
)

## any violence forest plots

# function to analyse violence overall with subgroups for exposure_type
perform_all_violence_analysis <- function(df, analysis, exposure) {
  
  # filter to infection
  filtered_df <- df %>%
    filter(outcome == "HIV prevalence", exposure_tf_bin == exposure) %>%
    filter(!is.na(.data[[var_names[[analysis]]$est]]))

  # study_num and effect_num columns
  filtered_df <- filtered_df %>%
    arrange(study) %>%
    mutate(
      study_num = cumsum(!duplicated(title)),
      effect_num = row_number()
    ) %>%
    ungroup()
    
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
    text.random = "Overall"
  )
  
  print(summary(result2))
  
  result2$TE.random <- result$b
  result2$lower.random <- result$ci.lb
  result2$upper.random <- result$ci.ub
  
  # create folder 
  filename <- paste0("Plots/prevalence/all violence/all_violence_", tolower(exposure), "_", analysis, ".png")
  
  png(filename = filename, width = 60, height = 55, units = "cm", res = 300)
  
  forest(
    result2,
    sortvar = filtered_df$study,
    xlim = c(0.2, 4),             
    leftcols = leftcols[[tolower(exposure)]], 
    leftlabs = leftlabs[[tolower(exposure)]],
    rightcols = rightcols,
    rightlabs = rightlabs,
    pooled.totals = TRUE,
    xintercept = 1,
    addrow.overall = TRUE,
    overall.hetstat = TRUE,
    overall = TRUE,
    labeltext = TRUE,
    col.subgroup = "black",
  )
}

# run for recent and ever
for (exposure in c("Recent", "Ever")) {
  for (analysis in analyses) {
    perform_all_violence_analysis(fsw_data_prev, analysis, exposure)
  }
}


  # eggers test
  eggers <- metabias(result2, method.bias = "linreg")
  eggers_p <- if (!is.null(eggers$p.value)) eggers$p.value else NA
  eggers_p_str <- if (!is.na(eggers_p)) sprintf("p = %.3f", eggers_p) else ""

  # funnel plot
  funnel_filename <- paste0("Plots/prevalence/all violence/funnel plots/all_violence_", tolower(exposure), "_", analysis, "_funnel.png")
  png(filename = funnel_filename, width = 15, height = 15, units = "cm", res = 300)
  funnel(result2, main = paste(analysis, "-", exposure, "\nEgger's test", eggers_p_str))
  dev.off()
## violence by type

# function to analyse by violent type and create forest plots
perform_analysis <- function(df, analysis, exposure, violence_type) {
  # filter the dataframe
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  filtered_df <- filtered_df %>% filter(!is.na(filtered_df[[var_names[[analysis]]$est]]))

  # study_num and effect_num columns
  filtered_df <- filtered_df %>%
    arrange(study) %>%
    mutate(
      study_num = cumsum(!duplicated(title)),
      effect_num = row_number()
    ) %>%
    ungroup()
  
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
  
  # create folder
  dir.create("Plots/prevalence/violence by type", recursive = TRUE, showWarnings = FALSE)
  
  # filename
  filename <- paste0(
    "Plots/prevalence/violence by type/",
    violence_type, "_", exposure, "_", analysis, ".png"
  )
  
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

# loop over each type of violence and create forest plots
for (violence_type in names(dataframes)) {
  for (analysis in analyses) {
    for (exposure in exposures) {
      df <- get(dataframes[[violence_type]][[exposure]])
      perform_analysis(df, analysis, exposure, violence_type)
    }
  }
}

## other violence types

# function for other violence types
create_other_violence_plots <- function(data, outcome_filter = "HIV prevalence") {
  for (exposure_tf in c("Recent", "Ever")) {
    for (analysis in analyses) {
      # filter
      filtered_df <- data %>%
        filter(outcome == outcome_filter, exposure_tf_bin == exposure_tf) %>%
        filter(!is.na(.data[[var_names[[analysis]]$est]]))
      
      # skip if number of studies is less than 2
      if (nrow(filtered_df) < 2) {
        message(paste("Not enough data for", analysis, "analysis,", exposure_tf, "exposure (n =", nrow(filtered_df), ")"))
        next
      }

      # study_num and effect_num columns
      filtered_df <- create_study_effect_nums(filtered_df)
      
      # covariance matrix assuming constant sampling correlation
      V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                        cluster = filtered_df$study_num,
                                        r = rho,
                                        smooth_vi = TRUE)
      
      # multilevel random effects model
      result <- rma.mv(filtered_df[[var_names[[analysis]]$est]], 
                       V = V_mat, 
                       random = ~ 1 | study_num / effect_num,
                       data = filtered_df,   
                       sparse = TRUE)       
      
      print(paste("Results for", analysis, "analysis,", exposure_tf, "exposure:"))
      print(result)
      print(exp(coef(result)))
      
      # forest plot
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
        text.random = "Overall"
      )
      
      print(summary(result2))
      
      # update with rma.mv results
      result2$TE.random <- result$b
      result2$lower.random <- result$ci.lb
      result2$upper.random <- result$ci.ub
      
      # forest plot filename
      filename <- paste0("Plots/prevalence/other violence/hiv_infection_other_violence_", tolower(exposure_tf), "_", analysis, ".png")
      
      # forest plot
      png(filename = filename, width = 45, height = 14, units = "cm", res = 600)
      
      forest(
        result2,
        sortvar = filtered_df$study,
        xlim = c(0.2, 4),             
        leftcols = c("study", "study_num", "effect_num", "exposure_definition_short", "exposure_time_frame", "perpetrator", "country"),
        leftlabs = c("Study", "Study number", "Effect number", "Exposure definition", "Exposure time frame", "Perpetrator", "Country"),
        rightcols = rightcols,
        rightlabs = rightlabs,
        pooled.totals = TRUE,
        xintercept = 1,
        addrow.overall = TRUE,
        overall.hetstat = TRUE,
        overall = TRUE,
        labeltext = TRUE,
        col.subgroup = "black"
      )
      
      dev.off()
      message(paste("Forest plot saved as:", filename))
    }
  }
}

# function to create all plots for both Recent and Ever
create_other_violence_plots(fsw_data_other)

## subgroup analysis 

# function for "recently exposed to any violence"
process_and_plot(
  data = fsw_data_prev_recent,
  data_name = "fsw_data_prev_recent",
  output_plot_filename = "Plots/subgroups/recent_any_violence_subgroup.png"
)

# function for "ever exposed to any violence"
process_and_plot(
  data = fsw_data_prev_ever,
  data_name = "fsw_data_prev_ever",
  output_plot_filename = "Plots/subgroups/ever_any_violence_subgroup.png"
)

## sensitivity analysis

# run for recent and ever violence and rho = 0.4
for (exposure in c("Recent", "Ever")) {
  for (analysis in analyses) {
    perform_all_violence_analysis_rho1(fsw_data_prev, analysis, exposure)
  }
}

# run for recent and ever violence and rho = 0.8
for (exposure in c("Recent", "Ever")) {
  for (analysis in analyses) {
    perform_all_violence_analysis_rho2(fsw_data_prev, analysis, exposure)
  }
}
