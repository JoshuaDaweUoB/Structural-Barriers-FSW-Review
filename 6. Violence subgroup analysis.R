# load packages 
pacman::p_load("meta", "metafor", "readxl", "tidyverse") 

# settings
settings.meta(CIbracket = "(") 
settings.meta(CIseparator = "-") 

# who subgroups
for (var in dfs_subgroup) {
  df <- get(var)
  meta_analysis_best_who(df, var)
}

for (var in dfs_subgroup) {
  df <- get(var)
  meta_analysis_best_who_sv(df, var)
}

# lic subgroups
for (var in dfs_subgroup) {
  df <- get(var)
  meta_analysis_best_ldc(df, var)
}

# recruitment subgroups
for (var in dfs_subgroup) {
  df <- get(var)
  meta_analysis_best_recruit(df, var)
}

# study setting subgroups
for (var in dfs_subgroup) {
  df <- get(var)
  meta_analysis_best_setting(df, var)
}

# perpetrator subgroups
for (var in dfs_subgroup) {
  df <- get(var)
  meta_analysis_best_perp(df, var)
}

# perpetrator subgroups
for (var in dfs_subgroup_perp) {
  df <- get(var)
  meta_analysis_best_perp_all(df, var)
}

# pre 2016 subgroups
for (var in dfs_subgroup) {
  df <- get(var)
  meta_analysis_best_pre2016(df, var)
}

# create sexual violence subgroup forest plots
for (df_name in dfs_subgroups) {
  subgroup_plot_forest_sv(df_name)
}

# create physical violence subgroup forest plots
for (df_name in dfs_subgroups) {
  subgroup_plot_forest_pv(df_name)
}

# create sexual or physical violence subgroup forest plots
for (df_name in dfs_subgroups) {
  subgroup_plot_forest_psv(df_name)
}
