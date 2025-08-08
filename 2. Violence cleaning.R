# load packages
pacman::p_load("readxl", "tidyverse")

## load data ##

# set wd
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence")

## load all data
fsw_data_all <- read_excel("Final data extraction.xlsx", "Final data extraction") 

# testing data
fsw_data_test <- fsw_data_all %>% filter(outcome == "HIV testing")

# prevalence data
fsw_data_prev <- fsw_data_all %>% filter(outcome == "HIV prevalence")

# Physical violence
fsw_data_pv_ever   <- fsw_data_all %>% filter(outcome == "HIV prevalence", exposure_tf_bin == "Ever", exposure_type == "Physical violence")
fsw_data_pv_recent <- fsw_data_all %>% filter(outcome == "HIV prevalence", exposure_tf_bin == "Recent", exposure_type == "Physical violence")

# Sexual violence
fsw_data_sv_ever   <- fsw_data_all %>% filter(outcome == "HIV prevalence", exposure_tf_bin == "Ever", exposure_type == "Sexual violence")
fsw_data_sv_recent <- fsw_data_all %>% filter(outcome == "HIV prevalence", exposure_tf_bin == "Recent", exposure_type == "Sexual violence")

# Physical or sexual violence
fsw_data_psv_ever   <- fsw_data_all %>% filter(outcome == "HIV prevalence", exposure_tf_bin == "Ever", exposure_type == "Physical and/or sexual violence")
fsw_data_psv_recent <- fsw_data_all %>% filter(outcome == "HIV prevalence", exposure_tf_bin == "Recent", exposure_type == "Physical and/or sexual violence")

# Other violence data
fsw_data_other      <- fsw_data_all %>% filter(outcome == "HIV prevalence", exposure_type == "Other violence")
fsw_data_other_ever <- fsw_data_all %>% filter(outcome == "HIV prevalence", exposure_type == "Other violence", exposure_tf_bin == "Ever")
fsw_data_other_rec  <- fsw_data_all %>% filter(outcome == "HIV prevalence", exposure_type == "Other violence", exposure_tf_bin == "Recent")

# ART data
fsw_data_art <- fsw_data_all %>% filter(outcome %in% c("ART use", "ART adherence"))

# VS data
fsw_data_vs <- fsw_data_all %>% filter(outcome == "Viral suppression")

# define dataframes
dfs <- c("fsw_data_test", "fsw_data_pv_ever", "fsw_data_pv_recent", "fsw_data_sv_ever", "fsw_data_sv_recent", "fsw_data_psv_ever", "fsw_data_psv_recent", "fsw_data_art", "fsw_data_vs", "fsw_data_other", "fsw_data_other_rec", "fsw_data_other_ever")

# data cleaning

# convert to numeric and log transform
for (var in dfs) {
  df <- get(var)
  df <- convert_to_numeric(df)
  df <- log_transform(df)
  df <- log_transform_var(df)
  assign(var, df)
}

# keep best estimates
for (var in dfs) {
  df <- get(var) %>%
    filter(!is.na(study))
  assign(var, df)
}

# clean unadjusted estimates
for (var in dfs) {
  df <- get(var)
  df <- format_violence_data_unadj(df)
  assign(var, df)
}

# clean adjusted estimates
for (var in dfs) {
  df <- get(var)
  df <- format_violence_data_adj(df)
  assign(var, df)
}

# clean best estimates
for (var in dfs) {
  df <- get(var)
  df <- format_violence_data(df)
  assign(var, df)
}
