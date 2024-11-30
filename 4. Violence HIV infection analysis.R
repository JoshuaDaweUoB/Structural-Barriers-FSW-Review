# load packages 
pacman::p_load("meta", "metafor", "readxl", "tidyverse") 

# settings
settings.meta(CIbracket = "(") 
settings.meta(CIseparator = "-") 

# overall plots

# overall dataframes
dfs_overall <- c("fsw_data_pv_ever", "fsw_data_pv_recent", "fsw_data_sv_ever", "fsw_data_sv_recent", "fsw_data_psv_ever", "fsw_data_psv_recent")

# unadjusted 
for (var in dfs_overall) {
  df <- get(var)
  meta_analysis_unadj(df, var)
}

# adjusted 
for (var in dfs_overall) {
  df <- get(var)
  meta_analysis_adj(df, var)
}

# combined/best estimate 
for (var in dfs_overall) {
  df <- get(var)
  meta_analysis_best(df, var)
}

# forest plot for all adjusted estimates
for (var in dfs_all) {
  df <- get(var)
  all_adjusted_forest(df, var)
}

# forest plot for all unadjusted estimates
for (var in dfs_all) {
  df <- get(var)
  all_unadjusted_forest(df, var)
}

# forest plot for all best estimates
for (var in dfs_all) {
  df <- get(var)
  all_best_forest(df, var)
}

# meta analysis using multi level model
for (var in dfs_overall) {
  df <- get(var)
  mlm_analysis_unadj(df, var)
}


















