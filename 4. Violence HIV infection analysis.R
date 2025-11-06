# load packages 
pacman::p_load("meta", "metafor", "readxl", "openxlsx", "tidyverse", "kableExtra", "robumeta", "clubSandwich", "grid", "png", "gridExtra") 

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
violence_types <- c("Physical", "Sexual", "Physical or sexual", "Other violence")

## modelling

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

violence_type_labels <- c(
  physical = "Physical",
  sexual = "Sexual",
  physical_sexual = "Physical or sexual",
  other = "Other violence"
)

# dataframes for each exposure and type of violence
dataframes <- list(
  physical = list(recent = "fsw_data_pv_recent", ever = "fsw_data_pv_ever"),
  sexual = list(recent = "fsw_data_sv_recent", ever = "fsw_data_sv_ever"),
  physical_sexual = list(recent = "fsw_data_psv_recent", ever = "fsw_data_psv_ever"),
  other = list(recent = "fsw_data_other_recent", ever = "fsw_data_other_ever")
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
  ),
  other = list(
    unadj = list(recent = "Plots/other_recent_unadj_csc.png", ever = "Plots/other_ever_unadj_csc.png"),
    adj = list(recent = "Plots/other_recent_adj_csc.png", ever = "Plots/other_ever_adj_csc.png"),
    best = list(recent = "Plots/other_recent_best_csc.png", ever = "Plots/other_ever_best_csc.png")
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
  png(filename = filename, width = 80, height = 60, units = "cm", res = 600) 
  
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

  dev.off()

  # eggers test
  eggers <- metabias(result2, method.bias = "linreg")
  eggers_p <- if (!is.null(eggers$p.value)) eggers$p.value else NA
  eggers_p_str <- if (!is.na(eggers_p)) sprintf("p = %.3f", eggers_p) else ""

  # funnel plot
  funnel_filename <- paste0("Plots/prevalence/all violence/funnel plots/all_violence_", tolower(exposure), "_", analysis, "_funnel.png")
  funnel_label <- paste0(analysis_labels[[analysis]], " - ", exposure_labels[[tolower(exposure)]])
  funnel_filename <- paste0("Plots/prevalence/all violence/funnel plots/", funnel_label, ".png")

  png(filename = funnel_filename, width = 15, height = 15, units = "cm", res = 300)
  funnel(result2, main = paste0(funnel_label, "\nEgger's test ", eggers_p_str))
  dev.off()
}

# run for recent and ever
for (exposure in c("Recent", "Ever")) {
  for (analysis in analyses) {
    perform_all_violence_analysis(fsw_data_prev, analysis, exposure)
  }
}

# function for combining funnel plots
get_all_violence_funnel_path <- function(analysis, exposure) {
  paste0(
    "Plots/prevalence/all violence/funnel plots/",
    analysis_labels[[analysis]], " - ", exposure_labels[[exposure]], ".png"
  )
}

# funnel plots list
all_funnel_imgs <- vector("list", length = length(exposures) * length(analyses))
idx <- 1
for (i in seq_along(exposures)) {
  for (j in seq_along(analyses)) {
    file <- get_all_violence_funnel_path(analyses[j], exposures[i])
    if (file.exists(file)) {
      all_funnel_imgs[[idx]] <- rasterGrob(readPNG(file), interpolate = TRUE)
    } else {
      all_funnel_imgs[[idx]] <- nullGrob()
    }
    idx <- idx + 1
  }
}

# combine figures 
combined_filename <- "Plots/prevalence/all violence/funnel plots/all_violence_funnel_grid.png"
png(combined_filename, width = 1800, height = 1200, res = 150)
grid.arrange(
  grobs = all_funnel_imgs,
  nrow = length(exposures),
  ncol = length(analyses),
  top = "All violence: Funnel plots"
)
dev.off()

## violence by type

# function to analyse by violent type and create forest plots
perform_analysis <- function(df, analysis, exposure, violence_type) {

  # filter the dataframe
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  filtered_df <- filtered_df %>% filter(!is.na(filtered_df[[var_names[[analysis]]$est]]))

  # skip if k <= 1
  if (nrow(filtered_df) <= 1) {
    message(paste("Skipping", analysis, exposure, "- not enough studies (k <=", nrow(filtered_df), ")"))
    return(NULL)
  }

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

eggers <- metabias(result2, method.bias = "linreg")
eggers_p <- if (!is.null(eggers$p.value)) eggers$p.value else NA
eggers_p_str <- if (!is.na(eggers_p)) sprintf("p = %.3f", eggers_p) else ""

# labels for funnel plot
violence_type_label <- violence_type_labels[[violence_type]]
analysis_label <- analysis_labels[[analysis]]
exposure_label <- exposure_labels[[tolower(exposure)]]
funnel_label <- paste0(violence_type_label, " - ", analysis_label, " - ", exposure_label)
funnel_filename <- paste0("Plots/prevalence/violence by type/funnel plots/", funnel_label, ".png")

png(filename = funnel_filename, width = 15, height = 15, units = "cm", res = 300)
funnel(result2, main = paste0(funnel_label, "\nEgger's test ", eggers_p_str))
dev.off()
}

# loop to create forest plots
for (violence_type in names(dataframes)) {
  for (analysis in analyses) {
    for (exposure in exposures) {
      df <- get(dataframes[[violence_type]][[exposure]])
      perform_analysis(df, analysis, exposure, violence_type)
    }
  }
}

## combining plots
get_funnel_path <- function(violence, analysis, exposure) {
  paste0(
    "Plots/prevalence/violence by type/funnel plots/",
    violence, " - ", analysis_labels[[analysis]], " - ", exposure_labels[[exposure]], ".png"
  )
}

for (violence in violence_types) {
  funnel_imgs <- vector("list", length = length(exposures) * length(analyses))
  idx <- 1
  for (i in seq_along(exposures)) { 
    for (j in seq_along(analyses)) { 
      file <- get_funnel_path(violence, analyses[j], exposures[i])
      if (file.exists(file)) {
        funnel_imgs[[idx]] <- rasterGrob(readPNG(file), interpolate = TRUE)
      } else {
        funnel_imgs[[idx]] <- nullGrob()
      }
      idx <- idx + 1
    }
  }
  
  combined_filename <- paste0("Plots/prevalence/violence by type/funnel plots/", violence, "_funnel_grid_recent_ever.png")
  png(combined_filename, width = 1800, height = 1200, res = 150)
  grid.arrange(
    grobs = funnel_imgs,
    nrow = length(exposures),
    ncol = length(analyses),
    top = violence
  )
  dev.off()
}











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
