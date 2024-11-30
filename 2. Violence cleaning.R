# load packages
pacman::p_load("readxl", "tidyverse") 

# Load data  --------------------------------------------------------------

## set wd
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence")

# load physical violence data
fsw_data_pv_ever <- read_excel("All violence studies.xlsx", "Physical violence - Ever") 
fsw_data_pv_recent <- read_excel("All violence studies.xlsx", "Physical violence - Recent") 
# fsw_data_pv_ever_perp <- read_excel("All violence studies.xlsx", "Ever data_perp") 
# fsw_data_pv_recent_perp <- read_excel("All violence studies.xlsx", "Recent data_perp") 

# load sexual violence data
fsw_data_sv_ever <- read_excel("All violence studies.xlsx", "Sexual violence - Ever")
fsw_data_sv_recent <- read_excel("All violence studies.xlsx", "Sexual violence - Recent")
# fsw_data_sv_ever_perp <- read_excel("All violence studies.xlsx", "Ever data_perp")
# fsw_data_sv_recent_perp <- read_excel("All violence studies.xlsx", "Recent data_perp")

# load physical & sexual violence data
fsw_data_psv_ever <- read_excel("All violence studies.xlsx", "Physical or sexual - Ever")
fsw_data_psv_recent <- read_excel("All violence studies.xlsx", "Physical or sexual - Recent")
# fsw_data_psv_ever_perp <- read_excel("All violence studies.xlsx", "Ever data_perp")
# fsw_data_psv_recent_perp <- read_excel("All violence studies.xlsx", "Recent data_perp")

# load art data
fsw_data_art <- read_excel("All violence studies.xlsx", "ART data - All")
fsw_data_art_pool <- read_excel("All violence studies.xlsx", "ART data - Pool")

# load vs data
fsw_data_vs <- read_excel("All violence studies.xlsx", "VS data - All")
fsw_data_vs_pool <- read_excel("All violence studies.xlsx", "VS data - Pool")

# load sv and hiv data
fsw_data_sv_hiv <- read_excel("HIV infection data.xlsx", "Sexual violence")

# load pv and hiv data
fsw_data_pv_hiv <- read_excel("HIV infection data.xlsx", "Physical violence")

# load sv and hiv data
fsw_data_psv_hiv <- read_excel("HIV infection data.xlsx", "Physical or sexual violence")

# load all pv data
fsw_data_pv_hiv_all <- read_excel("All violence studies.xlsx", "Physical violence - All")

# load all sv data
fsw_data_sv_hiv_all <- read_excel("All violence studies.xlsx", "Sexual violence - All")

# load all pv and sv data
fsw_data_psv_hiv_all <- read_excel("All violence studies.xlsx", "Physical or sexual - All")


# Lists -------------------------------------------------------------------

# define dataframes
dfs <- c("fsw_data_pv_ever", "fsw_data_pv_recent", "fsw_data_sv_ever", "fsw_data_sv_recent", "fsw_data_psv_ever", "fsw_data_psv_recent", "fsw_data_art", "fsw_data_art_pool", "fsw_data_vs", "fsw_data_vs_pool", "fsw_data_sv_hiv", "fsw_data_pv_hiv", "fsw_data_psv_hiv", "fsw_data_pv_hiv_all", "fsw_data_sv_hiv_all", "fsw_data_psv_hiv_all")

# all estimates dataframes
dfs_all <- c("fsw_data_pv_hiv_all", "fsw_data_sv_hiv_all", "fsw_data_psv_hiv_all")

# subgroup dataframes
dfs_subgroup <- c("fsw_data_pv_ever", "fsw_data_pv_recent", "fsw_data_sv_ever", "fsw_data_sv_recent", "fsw_data_psv_ever", "fsw_data_psv_recent")
# dfs_subgroup_perp <- c("fsw_data_pv_ever_perp", "fsw_data_pv_recent_perp", "fsw_data_sv_ever_perp", "fsw_data_sv_recent_perp", "fsw_data_psv_ever_perp", "fsw_data_psv_recent_perp")

# subgroup sheet names
sheet_names_subgroups <- c("WHO", "LIC", "Perp_partial", "Perp_all", "Recruit", "Pre_2016")

# subgroup dataframes
dfs_subgroups <- c("summary_violence_WHO", "summary_violence_LIC", "summary_violence_Perp_partial", "summary_violence_Perp_all", "summary_violence_Recruit", "summary_violence_Pre_2016")
  
# Data manipulation -------------------------------------------------------

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

# convert to numeric and log transform
for (var in dfs_subgroup) {
  df <- get(var)
  df <- convert_to_numeric(df)
  df <- log_transform(df)
  assign(var, df)
}

# keep best estimates
for (var in dfs_subgroup) {
  df <- get(var) %>%
    filter(!is.na(study)) %>%
    filter(use == "yes")
  assign(var, df)
}

# convert to numeric and log transform
for (var in dfs_subgroup_perp) {
  df <- get(var)
  df <- convert_to_numeric(df)
  df <- log_transform(df)
  assign(var, df)
}

# keep best estimates
for (var in dfs_subgroup_perp) {
  df <- get(var) %>%
    filter(!is.na(study)) %>%
    filter(use == "yes")
  assign(var, df)
}

# Clean overall data ------------------------------------------------------

# load ever overall summary data
summary_violence_ever <- read_excel("Violence estimates.xlsx", "Ever") # load data

summary_violence_ever <- transform(summary_violence_ever, 
                              effect_ln = log(effect),
                              lower_ln = log(lower),
                              upper_ln = log(upper)) 

summary_violence_ever$fuck <- round(summary_violence_ever$fuck, 0)

summary_violence_ever$lower_str <- format(summary_violence_ever$lower, nsmall = 2)
summary_violence_ever$upper_str <- format(summary_violence_ever$upper, nsmall = 2)
summary_violence_ever$or_95 <- paste0(summary_violence_ever$lower_str, "–", summary_violence_ever$upper_str)

summary_violence_ever$or_95_2 <- paste0(summary_violence_ever$effect, " ", "(", summary_violence_ever$or_95, ")")

# load recent overall summary data
summary_violence_rec <- read_excel("Violence estimates.xlsx", "Recent") # load data

summary_violence_rec <- transform(summary_violence_rec, 
                                   effect_ln = log(effect),
                                   lower_ln = log(lower),
                                   upper_ln = log(upper)) 

summary_violence_rec$fuck <- round(summary_violence_rec$fuck, 0)

summary_violence_rec$lower_str <- format(summary_violence_rec$lower, nsmall = 2)
summary_violence_rec$upper_str <- format(summary_violence_rec$upper, nsmall = 2)
summary_violence_rec$or_95 <- paste0(summary_violence_rec$lower_str, "–", summary_violence_rec$upper_str)

summary_violence_rec$or_95_2 <- paste0(summary_violence_rec$effect, " ", "(", summary_violence_rec$or_95, ")")



# loop through each subgroup and clean
for (sheet in sheet_names_subgroups) {
  df_name <- process_violence_data(sheet)
}
