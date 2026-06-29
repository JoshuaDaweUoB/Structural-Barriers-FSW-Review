# load packages
pacman::p_load("readxl", "writexl", "tidyverse", "metafor")

# hiv prev studies
fsw_data_prev <- create_study_effect_nums(fsw_data_prev)

# one row per study
fsw_data_prev <- fsw_data_prev %>%
  arrange(study) %>%
  group_by(study) %>%
  mutate(sequence = row_number()) %>%
  ungroup() %>%
  filter(sequence == 1)

n_unique_studies <- n_distinct(fsw_data_prev$study)
print(n_unique_studies)

# sequence study ids
fsw_data_all <- create_study_effect_nums(fsw_data_all)

# one row per study
all_studies <- fsw_data_all %>%
  arrange(study) %>%
  group_by(study) %>%
  mutate(sequence = row_number()) %>%
  ungroup() %>%
  filter(sequence == 1)

n_unique_studies <- n_distinct(all_studies$study)
print(n_unique_studies)

# categories for hiv_perc
hiv_prevalence_studies <- fsw_data_prev %>%
  mutate(
    hiv_perc = as.numeric(hiv_perc),
    hiv_perc_cat = case_when(
      hiv_perc >= 0 & hiv_perc < 0.1 ~ "0-0.1",
      hiv_perc >= 0.1 & hiv_perc < 0.25 ~ "0.1-0.25",
      hiv_perc >= 0.25 & hiv_perc < 0.50 ~ "0.25-0.50",
      hiv_perc >= 0.50 ~ "0.50+",
      TRUE ~ NA_character_
    ),
    hiv_perc_cat = factor(hiv_perc_cat, levels = c("0-0.1", "0.1-0.25", "0.25-0.50", "0.50+"))
  ) %>%
  mutate(
    sample_size_cat = case_when(
      analytical_sample_size > 0 & analytical_sample_size < 200 ~ "<200",
      analytical_sample_size >= 200  & analytical_sample_size < 500 ~ "200-499",
      analytical_sample_size >= 500  & analytical_sample_size < 1000 ~ "500-999",
      analytical_sample_size >= 1000 ~ "1000+",
    )
  )

# hiv_perc_cat by sample_size_cat
summary_table <- hiv_prevalence_studies %>%
  group_by(sample_size_cat) %>%
  count(hiv_perc_cat) %>%
  mutate(percentage = round((n / sum(n)) * 100, 2)) %>%
  ungroup()

print(summary_table)

# sample_size_quartile
hiv_prevalence_studies <- hiv_prevalence_studies %>%
  mutate(
    sample_size_quartile = ntile(analytical_sample_size, 4)
  )

# average hiv_perc within each sample_size_quartile
average_hiv_prev <- hiv_prevalence_studies %>%
  group_by(sample_size_quartile) %>%
  summarise(
    avg_hiv_prev = mean(as.numeric(hiv_perc), na.rm = TRUE),
    n = n()
  )

print(average_hiv_prev)

# study characteristics
all_studies <- all_studies %>%
  mutate(
    hiv_num = as.numeric(hiv_num),
    hiv_perc = as.numeric(hiv_perc) * 100,
    exposed_num = as.numeric(exposed_num),
    exposed_perc = as.numeric(exposed_perc) * 100
  ) %>%
  mutate(
    hiv_num = round(hiv_num, 0),
    hiv_perc = round(hiv_perc, 0),
    exposed_num = round(exposed_num, 0),
    exposed_perc = round(exposed_perc, 0)
  ) %>%
  mutate(
    hiv_perc = paste0(hiv_perc, "%"),
    exposed_perc = paste0(exposed_perc, "%")
  ) %>%
  mutate(
    hiv_num = ifelse(is.na(hiv_num), "NR", hiv_num),
    hiv_perc = ifelse(hiv_perc == "NA%", "NR", hiv_perc),
    exposed_num = ifelse(is.na(exposed_num), "NR", exposed_num),
    exposed_perc = ifelse(exposed_perc == "NA%", "NR", exposed_perc)
  )

all_studies <- all_studies %>%
  mutate(
    location = paste0(cities, ", ", country) 
  )

# number of unique values in the country variable
num_unique_countries <- all_studies %>%
  summarise(unique_countries = n_distinct(country))

print(num_unique_countries)

# all_studies dataframe
all_studies <- all_studies %>%
  # select
  select(
    study,
    title,
    country,
    whoregion,
    study_design,
    pub_type,
    analytical_sample_size,
    overall_sample_size,
    exposed_num,
    exposed_perc,
    hiv_num,
    hiv_perc,
    recruitment,
    representative,
    adjust,
    design,
    rob_score_3cat,
    rayyan,
    income_cat
  ) %>%

  # rename
  rename(
    "study" = study,
    "Title" = title,
    "Country" = country,
    "WHO region" = whoregion,
    "Study design" = study_design,
    "Publication type" = pub_type,
    "Analytic sample size" = analytical_sample_size,
    "Overall sample size" = overall_sample_size,
    "Exposed (n)" = exposed_num,
    "Exposed (%)" = exposed_perc,
    "HIV (n)" = hiv_num,
    "HIV (%)" = hiv_perc,
    "Recruitment strategy" = recruitment,
    "Representativeness" = representative,
    "Model type" = adjust,
    "Longitudinal design" = design,
    "ROB score" = rob_score_3cat,
    "Identified via search" = rayyan,
    "Country-level income" = income_cat
  )

# save
write_xlsx(all_studies, "All studies.xlsx")

# characteristics of included studies

# adjusted and unadjusted estimates
adjust_count_table <- fsw_data_all %>%
  group_by(adjust) %>%
  summarise(n = n()) %>%
  mutate(percentage = round((n / sum(n)) * 100, 1))
print(adjust_count_table)

# number and percentage of cross_sectional studies
cross_sectional_count_table <- fsw_data_all %>%
  group_by(design) %>%
  summarise(n = n()) %>%
  mutate(percentage = round((n / sum(n)) * 100, 1))
print(cross_sectional_count_table)

# number and percentage of representative studies
n_representative <- all_studies %>%
  summarise(
    n = sum(`Representativeness` == "Random sampling", na.rm = TRUE),
    percentage = round(100 * n / n(), 1)
  )
print(n_representative)

# table for WHO region
who_region_table <- all_studies %>%
  count(`WHO region`) %>%  
  mutate(percentage = n / sum(n) * 100)  

# table for study design
study_design_table <- all_studies %>%
  count(`Study design`) %>% 
  mutate(percentage = n / sum(n) * 100) 

# number of unique countries in all_studies
num_unique_countries <- all_studies %>% summarise(unique_countries = n_distinct(Country))
print(num_unique_countries)

# table for publication type
pub_type_table <- all_studies %>%
  count(`Publication type`) %>% 
  mutate(percentage = n / sum(n) * 100)  

# table for study quality
study_quality_table <- all_studies %>%
  count(`ROB score`) %>%
  mutate(percentage = n / sum(n) * 100)

# table for income cat
income_table <- all_studies %>%
  count(`Country-level income`) %>%
  mutate(percentage = n / sum(n) * 100)

# print
print(who_region_table)
print(pub_type_table)
print(study_design_table)
print(study_quality_table)
print(income_table)

# sample size
totals_table <- all_studies %>%
  summarise(
    total_sample_size = sum(as.numeric(`Analytic sample size`), na.rm = TRUE), 
    total_hiv_n = sum(as.numeric(`HIV (n)`), na.rm = TRUE)   
  )

print(totals_table)

## total studies for violence types
violence_df <- fsw_data_all %>%
  select(study, exposure_tf_bin, outcome, exposure_type, exposed_num, exposed_perc, analytical_sample_size)

# study and estimate counts
get_counts <- function(type) {
  df <- fsw_data_all %>% filter(exposure_type == type)
  n_studies <- n_distinct(df$study)
  n_estimates <- nrow(df)
  tibble(
    violence_type = type,
    n_studies = n_studies,
    n_estimates = n_estimates
  )
}

# summary table for each violence type
violence_types <- c(
  "Physical violence",
  "Sexual violence",
  "Physical and/or sexual violence",
  "Other violence"
)

violence_counts_table <- bind_rows(lapply(violence_types, get_counts))

# missing exposure data (no exposed_num and no exposed_perc)
missing_exposure <- violence_df %>%
  mutate(
    exposed_num = as.numeric(exposed_num),
    exposed_perc = as.numeric(exposed_perc),
    sample_size = as.numeric(analytical_sample_size)
  ) %>%
  filter(is.na(exposed_num) & is.na(exposed_perc))

missing_table <- tibble(
  violence_type = violence_types,
  n_missing_studies = sapply(
    violence_types,
    function(type) n_distinct(missing_exposure$study[missing_exposure$exposure_type == type])
  )
)

n_missing_total <- n_distinct(missing_exposure$study)

violence_summary_tables <- list(
  violence_counts = violence_counts_table,
  missing_counts = missing_table,
  n_missing_total = tibble(n_missing_total = n_missing_total)
)

# save
write_xlsx(violence_summary_tables, "Violence study and estimate counts.xlsx")

# unique studies reporting recent exposure
n_recent <- fsw_data_all %>% filter(exposure_tf_bin == "Recent") %>% summarise(n = n_distinct(study))
print(n_recent)

# unique studies reporting lifetime (ever) exposure
n_ever <- fsw_data_all %>% filter(exposure_tf_bin == "Ever") %>% summarise(n = n_distinct(study))
print(n_ever)

# relevant columns
violence_df <- fsw_data_all %>%
  filter(use_exposed == "yes") %>%
  select(study, exposure_tf_bin, outcome, exposure_type, exposure_definition_short, perpetrator, exposed_num, exposed_perc, analytical_sample_size, hiv_num)

# types of violence and their corresponding sheet names
violence_types <- c("Physical violence", "Sexual violence", "Physical or sexual", "Other violence")

all_violence_studies_list <- list()

# loop over each type of violence
for (violence in violence_types) {
  ever_stud <- violence_df %>%
    filter(exposure_type == violence, exposure_tf_bin == "Ever") %>%
    select(study)
  recent_stud <- violence_df %>%
    filter(exposure_type == violence, exposure_tf_bin == "Recent") %>%
    select(study)
  
  combined_studies <- bind_rows(ever_stud, recent_stud) %>%
    arrange(study) %>%
    distinct(study, .keep_all = TRUE)
  
  all_violence_studies_list[[violence]] <- combined_studies
}

# results for each type of violence
physical_violence_studies <- all_violence_studies_list[["Physical violence"]]
sexual_violence_studies <- all_violence_studies_list[["Sexual violence"]]
physical_or_sexual_violence_studies <- all_violence_studies_list[["Physical or sexual"]]

# recent and lifetime violence

# Physical violence studies
fsw_data_pv_ever_studies   <- violence_df %>% filter(exposure_tf_bin == "Ever", exposure_type == "Physical violence")
fsw_data_pv_recent_studies <- violence_df %>% filter(exposure_tf_bin == "Recent", exposure_type == "Physical violence")

# Sexual violence studies
fsw_data_sv_ever_studies   <- violence_df %>% filter(exposure_tf_bin == "Ever", exposure_type == "Sexual violence")
fsw_data_sv_recent_studies <- violence_df %>% filter(exposure_tf_bin == "Recent", exposure_type == "Sexual violence")

# Physical and/or sexual violence studies
fsw_data_psv_ever_studies   <- violence_df %>% filter(exposure_tf_bin == "Ever", exposure_type == "Physical and/or sexual violence")
fsw_data_psv_recent_studies <- violence_df %>% filter(exposure_tf_bin == "Recent", exposure_type == "Physical and/or sexual violence")

# Other violence studies
fsw_data_other_ever_studies   <- violence_df %>% filter(exposure_tf_bin == "Ever", exposure_type == "Other violence")
fsw_data_other_recent_studies <- violence_df %>% filter(exposure_tf_bin == "Recent", exposure_type == "Other violence")

# create list of dataframes
dfs_studies <- c(
  "fsw_data_pv_ever_studies", "fsw_data_pv_recent_studies",
  "fsw_data_sv_ever_studies", "fsw_data_sv_recent_studies",
  "fsw_data_psv_ever_studies", "fsw_data_psv_recent_studies",
  "fsw_data_other_ever_studies", "fsw_data_other_recent_studies"
)

# keep relevant columns
for (df_name in dfs_studies) {
  df <- get(df_name)
  df <- df %>%
    select(study, exposure_definition_short, perpetrator, exposed_num, exposed_perc, analytical_sample_size)
  assign(df_name, df, envir = .GlobalEnv)
}

# keep one type of violence per study
for (df_name in dfs_studies) {
  df <- get(df_name)
  # Ensure exposed_num is numeric for sequencing
  df$exposed_num <- as.numeric(df$exposed_num)
  df <- df %>%
    arrange(study, exposed_num) %>%
    group_by(study, exposed_num) %>%
    mutate(sequence = row_number()) %>%
    ungroup() %>%
    filter(sequence == 1) %>%
    select(-sequence)
  assign(df_name, df, envir = .GlobalEnv)
}

# create table for sheet
summary_table <- data.frame(
  dataframe = character(),
  num_studies = numeric(),
  num_estimates = numeric(),
  total_exposed_num = numeric(),
  total_sample_size = numeric(),
  exposed_percentage = numeric(),
  pooled_prev = numeric(),
  pooled_prev_ci_lb = numeric(),
  pooled_prev_ci_ub = numeric(),
  I2 = numeric(),
  stringsAsFactors = FALSE
)

# pooled proportions table
for (df_name in dfs_studies) {

  df <- get(df_name)
  df$exposed_num <- as.numeric(df$exposed_num)
  df$analytical_sample_size <- as.numeric(df$analytical_sample_size)

  df <- df[
    !is.na(df$exposed_num) &
    !is.na(df$analytical_sample_size) &
    df$analytical_sample_size > 0 &
    df$exposed_num >= 0,
  ]

  total_exposed_num <- sum(df$exposed_num)
  total_sample_size <- sum(df$analytical_sample_size)
  exposed_percentage <- 100 * total_exposed_num / total_sample_size
  num_studies <- length(unique(df$study))
  num_estimates <- nrow(df)

  if (nrow(df) >= 2) {

    escalc_res <- escalc(
      measure = "PLO",
      xi = df$exposed_num,
      ni = df$analytical_sample_size
    )

    escalc_res$effect_num <- seq_len(nrow(escalc_res))

    rma_res <- rma.mv(
      yi, vi,
      random = ~ 1 | study / effect_num,
      data = cbind(escalc_res, study = df$study)
    )

    pred <- predict(rma_res, transf = transf.ilogit)

    pooled_prev <- pred$pred
    pooled_prev_ci <- c(pred$ci.lb, pred$ci.ub)

    # I2 (total heterogeneity)
    sigma2 <- rma_res$sigma2
    mean_vi <- mean(rma_res$vi)
    total_var <- sum(sigma2) + mean_vi
    I2 <- 100 * sum(sigma2) / total_var

  } else {
    pooled_prev <- NA
    pooled_prev_ci <- c(NA, NA)
    I2 <- NA
  }

  summary_table <- rbind(
    summary_table,
    data.frame(
      dataframe = df_name,
      num_studies = num_studies,
      num_estimates = num_estimates,
      total_exposed_num = total_exposed_num,
      total_sample_size = total_sample_size,
      exposed_percentage = exposed_percentage,
      pooled_prev = pooled_prev,
      pooled_prev_ci_lb = pooled_prev_ci[1],
      pooled_prev_ci_ub = pooled_prev_ci[2],
      I2 = I2
    )
  )
}

print(summary_table)
write_xlsx(summary_table, "Summary Table.xlsx")

## study estimates tables
formatted_data <- fsw_data_all %>%
  select(
    author = author,
    year = year,
    violence_definition = exposure_definition_short,
    violence_time_frame = exposure_time_frame,
    exposed_num = exposed_num,
    exposed_perc = exposed_perc,
    perpetrator = perpetrator,
    outcome = outcome,
    adjusted_for = adjusted_for,
    unadj_effect = unadj_est,
    unadj_lower = un_lower,
    unadj_upper = un_upper,
    adj_effect = adj_est,
    adj_lower = adj_lower,
    adj_upper = adj_upper
  ) %>%
  mutate(
    exposed_num = as.numeric(exposed_num),
    exposed_perc = as.numeric(exposed_perc),
    author_year = paste0(author, " (", year, ")"),
    exposed = ifelse(!is.na(exposed_num) & !is.na(exposed_perc),
                     paste0(exposed_num, " (", round(exposed_perc * 100, 0), "%)"),
                     "NR"),
    perpetrator = ifelse(is.na(perpetrator), "Any perpetrator", perpetrator),
    adjusted_for = ifelse(is.na(adjusted_for), "NR", adjusted_for),
    unadj_ci = ifelse(!is.na(unadj_effect), 
                      paste0("OR: ", trimws(format(round(unadj_effect, 2), nsmall = 2)), " (", 
                             trimws(format(round(unadj_lower, 2), nsmall = 2)), "–", 
                             trimws(format(round(unadj_upper, 2), nsmall = 2)), ")"), 
                      NA),
    adj_ci = ifelse(!is.na(adj_effect), 
                    paste0("aOR: ", trimws(format(round(adj_effect, 2), nsmall = 2)), " (", 
                           trimws(format(round(adj_lower, 2), nsmall = 2)), "–", 
                           trimws(format(round(adj_upper, 2), nsmall = 2)), ")"), 
                    NA),
    effect_size = paste0(
      ifelse(!is.na(unadj_ci), unadj_ci, ""),
      ifelse(!is.na(unadj_ci) & !is.na(adj_ci), "; ", ""),
      ifelse(!is.na(adj_ci), adj_ci, "")
    )
  ) %>%
  select(
    `Author (year)` = author_year,
    `Violence definition` = violence_definition,
    `Violence time frame` = violence_time_frame,
    `Exposed, n (%)` = exposed,
    `Perpetrator` = perpetrator,
    `Outcome` = outcome,
    `Effect size (95% CI)` = effect_size,
    `Adjusted for` = adjusted_for
  )

violence_sheets <- split(formatted_data, fsw_data_all$exposure_type)
names(violence_sheets) <- names(violence_sheets) %>%
  str_replace_all("[\\[\\]:*?/\\\\]", "_")

write_xlsx(violence_sheets, "Formatted Violence Data by Type.xlsx")