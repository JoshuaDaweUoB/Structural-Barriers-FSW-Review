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
subgroup_columns <- c("ldc_bin", "lmic_bin", "hiv_decrim", "gbv_law", "pre_2016", "recruitment", "perpetrator", "who_region", "rob_score")
dataframes <- list(fsw_data_pv_recent, fsw_data_sv_recent, fsw_data_psv_recent, fsw_data_pv_ever, fsw_data_sv_ever, fsw_data_psv_ever)  
dataframe_names <- c("fsw_data_pv_recent", "fsw_data_sv_recent", "fsw_data_psv_recent", "fsw_data_pv_ever", "fsw_data_sv_ever", "fsw_data_psv_ever")
analyses <- c("unadj", "adj", "best")
exposures <- c("recent", "ever")

# filtered dataframe
filtered_df <- fsw_data_pv_recent %>% 
  filter(outcome == "HIV prevalence") %>% 
  filter(effect_best != "NA")
  
  # Create study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)

# Perform meta-analysis
  pv_ldc_recent <- metagen(TE = effect_best_ln,
                                  lower = effect_best_lower_ln,
                                  upper = effect_best_upper_ln,
                                  data = filtered_df,
                                  sm = "OR",
                                  method.tau = "DL",
                                  comb.fixed = FALSE,
                                  comb.random = FALSE, 
                                  backtransf = TRUE,
                                  byvar = ldc_bin,
                                  text.random = "Overall")
  
  # Print summary
  print(summary(pv_ldc_recent))
  
  # Save forest plot
  png(filename = "Plots/subgroups/ldc_pv_subgroup.png", width = 30, height = 28, units = "cm", res = 600)
  forest(pv_ldc_recent, 
         sortvar = study,
         xlim = c(0.2, 4),             
         leftcols = leftcols_recent, 
         leftlabs = leftlabs_recent,
         rightcols = rightcols, 
         rightlabs = rightlabs, 
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












# Filter the dataframe for outcome = "HIV prevalence" and effect_best != "NA"
  filtered_df <- fsw_data_pv_recent %>% 
    filter(outcome == "HIV prevalence") %>% 
    filter(effect_best != "NA")
  
  # Debug: Print the number of rows after filtering
  message(paste("Filtered dataframe rows:", nrow(filtered_df)))
  
  # Loop through each subgroup column
  for (column in subgroup_columns) {
    # Debug: Print the current column being processed
    message(paste("Processing column:", column))
    
    # Get unique levels of the current column
    unique_levels <- unique(filtered_df[[column]])
    
    # Debug: Print the unique levels
    print(unique_levels)