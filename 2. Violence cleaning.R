# load packages
pacman::p_load("readxl", "tidyverse")

## load data ##

# set wd
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence")

# load testing data
fsw_data_test <- read_excel("All violence studies.xlsx", "Testing data - All") 

# Load physical violence data
fsw_data_pv_ever <- read_excel("All violence studies.xlsx", "Physical violence - Ever") %>%
  filter(outcome == "HIV prevalence")
fsw_data_pv_recent <- read_excel("All violence studies.xlsx", "Physical violence - Recent") %>%
  filter(outcome == "HIV prevalence")

# Load sexual violence data
fsw_data_sv_ever <- read_excel("All violence studies.xlsx", "Sexual violence - Ever") %>%
  filter(outcome == "HIV prevalence")
fsw_data_sv_recent <- read_excel("All violence studies.xlsx", "Sexual violence - Recent") %>%
  filter(outcome == "HIV prevalence")

# Load physical & sexual violence data
fsw_data_psv_ever <- read_excel("All violence studies.xlsx", "Physical or sexual - Ever") %>%
  filter(outcome == "HIV prevalence")
fsw_data_psv_recent <- read_excel("All violence studies.xlsx", "Physical or sexual - Recent") %>%
  filter(outcome == "HIV prevalence")

# load art data
fsw_data_art <- read_excel("All violence studies.xlsx", "ART data - All")

# load vs data
fsw_data_vs <- read_excel("All violence studies.xlsx", "VS data - All")

# load other violence data
fsw_data_other <- read_excel("All violence studies.xlsx", "Other violence - All")

# define dataframes
dfs <- c("fsw_data_test", "fsw_data_pv_ever", "fsw_data_pv_recent", "fsw_data_sv_ever", "fsw_data_sv_recent", "fsw_data_psv_ever", "fsw_data_psv_recent", "fsw_data_art", "fsw_data_vs", "fsw_data_other")

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
    filter(!is.na(study)) %>%
    filter(use == "yes")
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
