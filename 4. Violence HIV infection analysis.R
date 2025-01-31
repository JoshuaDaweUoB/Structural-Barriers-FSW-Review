# load packages 
pacman::p_load("meta", "metafor", "readxl", "openxlsx", "tidyverse", "kableExtra", "robumeta", "clubSandwich") 

# set working directory
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence")

# settings
settings.meta(CIbracket = "(") 
settings.meta(CIseparator = "-") 

# columns
leftcols_recent <- c("study", "exposure_definition_short", "exposure_time_frame", "perpetrator", "country")
leftlabs_recent <- c("Study", "Exposure definition", "Exposure time frame", "Perpetrator", "Country")
leftcols_lifetime <- c("study", "exposure_definition_short", "perpetrator", "country")
leftlabs_lifetime <- c("Study", "Exposure definition", "Perpetrator", "Country")
rightcols <- c("effect", "ci")
rightlabs = c("Estimate", "95% CI")

# lists for loops and functions
subgroup_columns <- c("ldc_bin", "pre_2016", "recruitment", "perpetrator", "who_region")
dataframes <- list(fsw_data_pv_recent, fsw_data_sv_recent, fsw_data_psv_recent)
dataframe_names <- c("fsw_data_pv_recent", "fsw_data_sv_recent", "fsw_data_psv_recent")
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
  png(filename = filename, width = 38, height = 18, units = "cm", res = 600)
  
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

#### subgroup analyses

## recent violence

# List to store results
results_list <- list()
skipped_results_list <- list()

# Loop through each dataframe
for (i in 1:length(dataframes)) {
  
  # Get the current dataframe and name
  current_df <- dataframes[[i]]
  current_df_name <- dataframe_names[i]
  
  # Filter the dataframe for outcome = "HIV prevalence" and effect_best != "NA"
  filtered_df <- current_df %>% 
    filter(outcome == "HIV prevalence") %>% 
    filter(effect_best != "NA")
  
  # Loop through each subgroup column
  for (column in subgroup_columns) {
    # Get unique levels of the current column
    unique_levels <- unique(filtered_df[[column]])
    
    # Loop through each unique level in the current column
    for (level in unique_levels) {
      # Filter data for the current level
      subgroup_df <- filtered_df[filtered_df[[column]] == level, ]
      
      # Handle cases with only one study
      if (nrow(subgroup_df) == 1) {
        message(paste("Only one study for", current_df_name, column, "=", level))
        results_list[[length(results_list) + 1]] <- data.frame(
          dataframe = current_df_name,
          column = column,
          level = level,
          pooled_OR = exp(subgroup_df$effect_best_ln),
          lower_CI = exp(subgroup_df$effect_best_lower_ln),
          upper_CI = exp(subgroup_df$effect_best_upper_ln),
          estimates = nrow(subgroup_df),
          studies = length(unique(subgroup_df$study_num))
        )
        next
      }
      
      # Create a covariance matrix for the subgroup
      V_mat <- impute_covariance_matrix(subgroup_df$effect_best_var_ln,
                                        cluster = subgroup_df$study_num,
                                        r = rho,
                                        smooth_vi = TRUE)
      
      # Fit the multilevel random effects model for the subgroup
      result <- rma.mv(effect_best_ln, 
                       V = V_mat, 
                       random = ~ 1 | study_num / effect_num,
                       data = subgroup_df,   
                       control = list(rel.tol = 1e-8),
                       sparse = TRUE)
      
      # Use `result` to update the `metagen` object
      result2 <- metagen(TE = subgroup_df$effect_best_ln,
                         lower = subgroup_df$effect_best_lower_ln,
                         upper = subgroup_df$effect_best_upper_ln,
                         studlab = subgroup_df$study,
                         data = subgroup_df,
                         sm = "OR",
                         method.tau = "REML",
                         common = FALSE,
                         random = TRUE, 
                         backtransf = TRUE,
                         text.random = "Overall")
      
      # Update result2 with values from the multilevel model
      result2$TE.random <- result$b
      result2$lower.random <- result$ci.lb
      result2$upper.random <- result$ci.ub
      
      # Store the results in the list
      results_list[[length(results_list) + 1]] <- data.frame(
        dataframe = current_df_name,
        column = column,
        level = level,
        pooled_OR = exp(result$b),
        lower_CI = exp(result$ci.lb),
        upper_CI = exp(result$ci.ub),
        estimates = nrow(subgroup_df),
        studies = length(unique(subgroup_df$study_num))
      )
      
      # Create a unique filename for the forest plot
      filename <- paste0("Plots/subgroups/", current_df_name, "_", column, "_", level, ".png")
      
      # Save the forest plot as a PNG file
      png(filename = filename, width = 38, height = 18, units = "cm", res = 600)
      forest(result2,
             sortvar = subgroup_df$study,
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
}

# Combine the results list and skipped results list
all_results_list <- c(results_list, skipped_results_list)

# Convert the combined results list to a dataframe
results_df <- do.call(rbind, all_results_list)

# transform ORs and CIs to log scale
results_df <- results_df %>%
    mutate(
      pooled_OR = log(pooled_OR),
      lower_CI = log(lower_CI),
      upper_CI = log(upper_CI)
    )

# Recode the values in the column variable
results_df <- results_df %>%
    mutate(column = recode(column,
                           ldc_bin = "Least developed country",
                           pre_2016 = "Published before 2016",
                           recruitment = "Recruitment strategy",
                           who_region = "WHO region",
                           perpetrator = "Perpetrator"))

# Recode the values in the column variable
results_df <- results_df %>%
    mutate(level = recode(level,
                           no = "No",
                           yes = "Yes",
                           "FALSE" = "No",
                           "TRUE" = "Yes"))                           

# Split the results_df into three dataframes
fsw_data_pv_recent_results <- results_df[results_df$dataframe == "fsw_data_pv_recent", ]
fsw_data_sv_recent_results <- results_df[results_df$dataframe == "fsw_data_sv_recent", ]
fsw_data_psv_recent_results <- results_df[results_df$dataframe == "fsw_data_psv_recent", ]

# Combine the three dataframes into one
combined_results_df <- bind_rows(
  fsw_data_pv_recent_results,
  fsw_data_sv_recent_results,
  fsw_data_psv_recent_results
)

# Function to create a forest plot
forest_plot <- function(df, title) {
  # Create a metagen object
  meta_analysis <- metagen(
    TE = df$pooled_OR,
    lower = df$lower_CI,
    upper = df$upper_CI,
    studlab = df$level,
    data = df,
    sm = "OR",
    method.tau = "REML",
    common = FALSE,
    random = FALSE,
    backtransf = TRUE,
    subgroup = df$column
  )
  
  # Create the forest plot
  forest(meta_analysis,
         sortvar = df$level,
         xlim = c(0.2, 4),
         leftcols = c("studlab"),
         leftlabs = c("Subgroup level"),
         rightcols = c("effect", "ci"),
         rightlabs = c("Estimate", "95% CI"),
         pooled.totals = FALSE,
         addfit = FALSE,
         xintercept = 1,
         addrow.overall = FALSE,
         overall.hetstat = FALSE,
         overall = FALSE,
         addrow.subgroups = FALSE,
         labeltext = TRUE,
         col.subgroup = "black",
         print.subgroup.name = FALSE,         
         main = title)
}

# Create forest plots for each dataframe
png(filename = "Plots/fsw_data_pv_recent_forest_plot.png", width = 25, height = 18, units = "cm", res = 600)
forest_plot(fsw_data_pv_recent_results, "FSW Data PV Recent Results")
dev.off()

png(filename = "Plots/fsw_data_sv_recent_forest_plot.png", width = 25, height = 18, units = "cm", res = 600)
forest_plot(fsw_data_sv_recent_results, "FSW Data SV Recent Results")
dev.off()

png(filename = "Plots/fsw_data_psv_recent_forest_plot.png", width = 25, height = 18, units = "cm", res = 600)
forest_plot(fsw_data_psv_recent_results, "FSW Data PSV Recent Results")
dev.off()


