# load packages 
pacman::p_load("meta", "metafor", "readxl", "openxlsx", "tidyverse", "kableExtra", "robumeta", "clubSandwich", "grid", "png", "gridExtra") 

# set working directory
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Analysis/Structural-Barriers-FSW-Review")

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
  
  filtered_df <- df %>%
    filter(outcome == "HIV prevalence", exposure_tf_bin == exposure) %>%
    filter(!is.na(.data[[var_names[[analysis]]$est]]))

  filtered_df <- filtered_df %>%
    arrange(study) %>%
    mutate(
      study_num = cumsum(!duplicated(title)),
      effect_num = row_number()
    ) %>%
    ungroup()
    
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
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
  
  filename <- paste0("Plots/prevalence/all violence/all_violence_", tolower(exposure), "_", analysis, ".png")
  png(filename = filename, width = 80, height = 60, units = "cm", res = 300) 
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

  eggers <- metabias(result2, method.bias = "linreg")
  eggers_p <- if (!is.null(eggers$p.value)) eggers$p.value else NA
  eggers_p_str <- if (!is.na(eggers_p)) sprintf("p = %.3f", eggers_p) else ""

  funnel_label <- paste0(analysis_labels[[analysis]], " - ", exposure_labels[[tolower(exposure)]])
  funnel_filename <- paste0("Plots/prevalence/all violence/funnel plots/", funnel_label, ".png")

  png(filename = funnel_filename, width = 15, height = 15, units = "cm", res = 300)
  funnel(result2, main = paste0(funnel_label, "\nEgger's test ", eggers_p_str))
  dev.off()

  return(result2)
}

# run main analysis and store result2 objects
all_violence_results <- list()
for (exposure in c("Recent", "Ever")) {
  for (analysis in analyses) {
    all_violence_results[[paste(exposure, analysis)]] <- perform_all_violence_analysis(fsw_data_prev, analysis, exposure)
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

# trim and fill
perform_all_violence_trimfill <- function(result2, analysis, exposure) {
  
  funnel_label <- paste0(analysis_labels[[analysis]], " - ", exposure_labels[[tolower(exposure)]])
  
  tf_result <- trimfill(result2)
  print(summary(tf_result))
  print(paste("Studies trimmed and filled:", tf_result$k0))
  print("Original random effect estimate (OR):")
  print(exp(result2$TE.random))
  print("Trim and fill adjusted estimate (OR):")
  print(exp(tf_result$TE.random))
  
  tf_funnel_label <- paste0(funnel_label, " - Trim and Fill")
  tf_funnel_filename <- paste0("Plots/prevalence/all violence/funnel plots/", tf_funnel_label, ".png")
  
  png(filename = tf_funnel_filename, width = 15, height = 15, units = "cm", res = 300)
  funnel(tf_result, main = paste0(tf_funnel_label, "\nStudies trimmed: ", tf_result$k0))
  dev.off()

  trim_fill_results <<- rbind(trim_fill_results, data.frame(
    analysis_type = "All violence",
    violence_type = NA,
    exposure = exposure,
    analysis = analysis,
    original_or = exp(result2$TE.random),
    original_ci_lower = exp(result2$lower.random),
    original_ci_upper = exp(result2$upper.random),
    tf_or = exp(tf_result$TE.random),
    tf_ci_lower = exp(tf_result$lower.random),
    tf_ci_upper = exp(tf_result$upper.random),
    studies_trimmed = tf_result$k0
  ))
}

# create dataframe to store trim and fill results
trim_fill_results <- data.frame(
  analysis_type = character(),
  violence_type = character(),
  exposure = character(),
  analysis = character(),
  original_or = numeric(),
  original_ci_lower = numeric(),
  original_ci_upper = numeric(),
  tf_or = numeric(),
  tf_ci_lower = numeric(),
  tf_ci_upper = numeric(),
  studies_trimmed = numeric(),
  stringsAsFactors = FALSE
)

# run trim and fill
for (exposure in c("Recent", "Ever")) {
  for (analysis in analyses) {
    perform_all_violence_trimfill(all_violence_results[[paste(exposure, analysis)]], analysis, exposure)
  }
}

get_all_violence_tf_funnel_path <- function(analysis, exposure) {
  paste0(
    "Plots/prevalence/all violence/funnel plots/",
    analysis_labels[[analysis]], " - ", exposure_labels[[exposure]], " - Trim and Fill.png"
  )
}

# trim and fill funnel plots list
all_tf_funnel_imgs <- vector("list", length = length(exposures) * length(analyses))
idx <- 1
for (i in seq_along(exposures)) {
  for (j in seq_along(analyses)) {
    file <- get_all_violence_tf_funnel_path(analyses[j], exposures[i])
    if (file.exists(file)) {
      all_tf_funnel_imgs[[idx]] <- rasterGrob(readPNG(file), interpolate = TRUE)
    } else {
      all_tf_funnel_imgs[[idx]] <- nullGrob()
    }
    idx <- idx + 1
  }
}

# combine trim and fill figures
tf_combined_filename <- "Plots/prevalence/all violence/funnel plots/all_violence_funnel_grid_TrimFill.png"
png(tf_combined_filename, width = 1800, height = 1200, res = 150)
grid.arrange(
  grobs = all_tf_funnel_imgs,
  nrow = length(exposures),
  ncol = length(analyses),
  top = "All violence: Funnel plots - Trim and Fill"
)
dev.off()

## violence by type

# function to analyse by violent type and create forest plots
perform_analysis <- function(df, analysis, exposure, violence_type) {

  # filter the dataframe
  filtered_df <- df %>% filter(outcome == "HIV prevalence")
  filtered_df <- filtered_df %>% filter(!is.na(filtered_df[[var_names[[analysis]]$est]]))

  # make sure no dataframes incorrectly filtered
  if (nrow(filtered_df) == 0) {
    message("Skipping: ", violence_type, " ", exposure, " ", analysis, " (no data)")
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
  result <- rma.mv(
    filtered_df[[var_names[[analysis]]$est]],
    V = V_mat,
    random = ~ 1 | study_num / effect_num,
    data = filtered_df,
    sparse = TRUE,
    control = list(
      optimizer = "nlminb",
      iter.max = 10000,
      eval.max = 10000,
      rel.tol = 1e-8
    )
  )    
  
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
  filename <- paste0(
    "Plots/prevalence/violence by type/",
    violence_type, "_", exposure, "_", analysis, ".png"
  )
  
  png(filename = filename, width = 45, height = 22, units = "cm", res = 600)
  
  forest(result2,
         sortvar = filtered_df$study,
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

  # trim and fill analysis
  tf_result <- trimfill(result2)
  print(summary(tf_result))
  print(paste("Studies trimmed and filled:", tf_result$k0))
  
  # compare original vs adjusted estimates
  print("Original random effect estimate (OR):")
  print(exp(result2$TE.random))
  print("Trim and fill adjusted estimate (OR):")
  print(exp(tf_result$TE.random))
  
  # trim and fill funnel plot
  tf_funnel_label <- paste0(violence_type_label, " - ", analysis_label, " - ", exposure_label, " - Trim and Fill")
  tf_funnel_filename <- paste0("Plots/prevalence/violence by type/funnel plots/", tf_funnel_label, ".png")
  
  png(filename = tf_funnel_filename, width = 15, height = 15, units = "cm", res = 300)
  funnel(tf_result, main = paste0(tf_funnel_label, "\nStudies trimmed: ", tf_result$k0))
  dev.off()

  trim_fill_results <<- rbind(trim_fill_results, data.frame(
  analysis_type = "By violence type",
  violence_type = violence_type,
  exposure = exposure,
  analysis = analysis,
  original_or = exp(result2$TE.random),
  original_ci_lower = exp(result2$lower.random),
  original_ci_upper = exp(result2$upper.random),
  tf_or = exp(tf_result$TE.random),
  tf_ci_lower = exp(tf_result$lower.random),
  tf_ci_upper = exp(tf_result$upper.random),
  studies_trimmed = tf_result$k0
))
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

get_tf_funnel_path <- function(violence, analysis, exposure) {
  paste0(
    "Plots/prevalence/violence by type/funnel plots/",
    violence, " - ", analysis_labels[[analysis]], " - ", exposure_labels[[exposure]], " - Trim and Fill.png"
  )
}

for (violence in violence_types) {
  # original funnel plots
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
  
  # trim and fill combined plots
  tf_funnel_imgs <- vector("list", length = length(exposures) * length(analyses))
  idx <- 1
  for (i in seq_along(exposures)) { 
    for (j in seq_along(analyses)) { 
      file <- get_tf_funnel_path(violence, analyses[j], exposures[i])
      if (file.exists(file)) {
        tf_funnel_imgs[[idx]] <- rasterGrob(readPNG(file), interpolate = TRUE)
      } else {
        tf_funnel_imgs[[idx]] <- nullGrob()
      }
      idx <- idx + 1
    }
  }
  
  tf_combined_filename <- paste0("Plots/prevalence/violence by type/funnel plots/", violence, "_funnel_grid_recent_ever_TrimFill.png")
  png(tf_combined_filename, width = 1800, height = 1200, res = 150)
  grid.arrange(
    grobs = tf_funnel_imgs,
    nrow = length(exposures),
    ncol = length(analyses),
    top = paste0(violence, " - Trim and Fill")
  )
  dev.off()
}

# format table for output
trim_fill_table <- trim_fill_results %>%
  mutate(
    # Create a display name for violence type
    violence_display = ifelse(is.na(violence_type), "All violence", violence_type),
    original_effect = paste0(
      sprintf("%.2f", original_or), 
      " (", sprintf("%.2f", original_ci_lower), 
      "-", sprintf("%.2f", original_ci_upper), ")"
    ),
    tf_effect = paste0(
      sprintf("%.2f", tf_or), 
      " (", sprintf("%.2f", tf_ci_lower), 
      "-", sprintf("%.2f", tf_ci_upper), ")"
    )
  ) %>%
  select(
    analysis_type, violence_display, exposure, analysis,
    original_effect, tf_effect, studies_trimmed
  ) %>%
  rename(
    "Analysis Type" = analysis_type,
    "Violence Type" = violence_display,
    "Exposure" = exposure,
    "Model" = analysis,
    "Original OR (95% CI)" = original_effect,
    "Trim & Fill OR (95% CI)" = tf_effect,
    "Studies Trimmed" = studies_trimmed
  )

# create workbook
wb <- createWorkbook()
addWorksheet(wb, "Trim and Fill Results")

# add title
writeData(wb, sheet = 1, x = "Comparison of Original vs Trim and Fill Pooled Effects", startRow = 1)

# add table data starting at row 3
writeData(wb, sheet = 1, x = trim_fill_table, startRow = 3)

# format header row (no background color)
headerStyle <- createStyle(textDecoration = "bold", 
                           halign = "center", valign = "center", wrapText = TRUE)
for (col in 1:ncol(trim_fill_table)) {
  addStyle(wb, sheet = 1, style = headerStyle, rows = 3, cols = col)
}

# auto-fit column widths
setColWidths(wb, sheet = 1, cols = 1:ncol(trim_fill_table), widths = "auto")

# save workbook
saveWorkbook(wb, "Plots/trim_fill_comparison_table.xlsx", overwrite = TRUE)

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

## standalone forest plot: physical and/or sexual violence, recent exposure, excluding Budhwani


analysis <- "best"

filtered_df <- fsw_data_psv_recent_nobud %>%
  filter(outcome == "HIV prevalence") %>%
  filter(!is.na(.data[[var_names[[analysis]]$est]])) %>%
  arrange(study) %>%
  mutate(
    study_num = cumsum(!duplicated(title)),
    effect_num = row_number()
  ) %>%
  ungroup()

V_mat <- impute_covariance_matrix(
  filtered_df[[var_names[[analysis]]$var]],
  cluster = filtered_df$study_num,
  r = rho,
  smooth_vi = TRUE
)

result <- rma.mv(
  filtered_df[[var_names[[analysis]]$est]],
  V = V_mat,
  random = ~ 1 | study_num / effect_num,
  data = filtered_df,
  sparse = TRUE,
  control = list(
    optimizer = "nlminb",
    iter.max = 10000,
    eval.max = 10000,
    rel.tol = 1e-8
  )
)

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

result2$TE.random <- result$b
result2$lower.random <- result$ci.lb
result2$upper.random <- result$ci.ub

filename <- "Plots/prevalence/violence by type/psv_recent_nobud_best.png"
png(filename = filename, width = 45, height = 22, units = "cm", res = 600)
forest(
  result2,
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
  col.subgroup = "black"
)
dev.off()
