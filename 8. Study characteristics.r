# Load required packages
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

## this is a mistake - need to filter by prevalence first
## also need to think about exp_prop 

# Create categories for hiv_perc2
hiv_prevalence_studies <- fsw_data_prev %>%
  mutate(
    hiv_perc_cat = case_when(
      hiv_perc2 >= 0 & hiv_perc2 < 0.1 ~ "0-0.1",
      hiv_perc2 >= 0.1 & hiv_perc2 < 0.25 ~ "0.1-0.25",
      hiv_perc2 >= 0.25 & hiv_perc2 < 0.50 ~ "0.25-0.50",
      hiv_perc2 >= 0.50 ~ "0.50+",
      TRUE ~ NA_character_
    ),
    hiv_perc_cat = factor(hiv_perc_cat, levels = c("0-0.1", "0.1-0.25", "0.25-0.50", "0.50+"))
  )

# Summarize hiv_perc_cat by sample_size_cat
summary_table <- hiv_prevalence_studies %>%
  group_by(sample_size_cat, hiv_perc_cat) %>%
  summarise(
    count = n(),  # Count the number of studies in each category
    percentage = round((n() / sum(n())) * 100, 2)  # Calculate percentage within each sample_size_cat
  ) %>%
  ungroup()

# Print the summary table
print(summary_table)

# Create sample_size_quartile (quartiles for sample_size)
hiv_prevalence_studies <- hiv_prevalence_studies %>%
  mutate(
    sample_size_quartile = ntile(sample_size, 4)  # Divide sample_size into 4 quartiles
  )

# Calculate the average hiv_perc2 within each sample_size_quartile
average_hiv_prev <- hiv_prevalence_studies %>%
  group_by(sample_size_quartile) %>%
  summarise(
    avg_hiv_prev = mean(hiv_perc2, na.rm = TRUE),  # Calculate the mean, ignoring NAs
    n = n()  # Count the number of studies in each quartile
  )

# Print the result
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

# Calculate the number of unique values in the country variable
num_unique_countries <- all_studies %>%
  summarise(unique_countries = n_distinct(country))

# Print the result
print(num_unique_countries)

# Update the all_studies dataframe
all_studies <- all_studies %>%
  # Select and reorder the columns
  select(
    study,
    country,
    who_region,
    design,
    pub_type,
    sample_size,
    exposed_num,
    exposed_perc,
    hiv_num,
    hiv_perc,
    recruitment,
    rob_score
  ) %>%

  # Rename the columns
  rename(
    "study" = study,
    "Country" = country,
    "WHO region" = who_region,
    "Study design" = design,
    "Publication type" = pub_type,
    "Sample size" = sample_size,
    "Exposed (n)" = exposed_num,
    "Exposed (%)" = exposed_perc,
    "HIV (n)" = hiv_num,
    "HIV (%)" = hiv_perc,
    "Recruitment strategy" = recruitment,
    "ROB score" = rob_score
  )

# save the all_studies dataframe to an Excel file
write_xlsx(all_studies, "All studies.xlsx")

# characteristics of included studies

# Descriptive table for WHO region
who_region_table <- all_studies %>%
  count(`WHO region`) %>%  # Use backticks for column names with spaces
  mutate(percentage = n / sum(n) * 100)  # Calculate percentages

# Descriptive table for study design
study_design_table <- all_studies %>%
  count(`Study design`) %>%  # Use backticks for column names with spaces
  mutate(percentage = n / sum(n) * 100)  # Calculate percentages

# Count the number of unique countries in all_studies
num_unique_countries <- all_studies %>% summarise(unique_countries = n_distinct(Country))
print(num_unique_countries)

# Descriptive table for "Publication type"
pub_type_table <- all_studies %>%
  count(`Publication type`) %>%  # Use backticks for column names with spaces
  mutate(percentage = n / sum(n) * 100)  # Calculate percentages

# Descriptive table for study quality
study_quality_table <- all_studies %>%
  count(`ROB score`) %>%
  mutate(percentage = n / sum(n) * 100)
 
# Print the tables
print(who_region_table)
print(pub_type_table)
print(study_design_table)
print(study_quality_table)

# Calculate total Sample size and total HIV (n)
totals_table <- all_studies %>%
  summarise(
    total_sample_size = sum(as.numeric(`Sample size`), na.rm = TRUE),  # Convert Sample size to numeric and sum
    total_hiv_n = sum(as.numeric(`HIV (n)`), na.rm = TRUE)             # Convert HIV (n) to numeric and sum
  )

# Print the totals table
print(totals_table)

## total studies for violence types

# Function to get study and estimate counts
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

# Create summary table for each violence type
violence_types <- c(
  "Physical violence",
  "Sexual violence",
  "Physical and/or sexual violence",
  "Other violence"
)

violence_counts_table <- bind_rows(lapply(violence_types, get_counts))

# Studies missing exposure data (no exposed_num and no exposed_perc)
missing_exposure <- violence_df %>%
  mutate(
    exposed_num = as.numeric(exposed_num),
    exposed_perc = as.numeric(exposed_perc),
    sample_size = as.numeric(sample_size)
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

# Combine all results into one list for easy export
violence_summary_tables <- list(
  violence_counts = violence_counts_table,
  missing_counts = missing_table,
  n_missing_total = tibble(n_missing_total = n_missing_total)
)

# Write to Excel (each table as a sheet)
write_xlsx(violence_summary_tables, "Violence study and estimate counts.xlsx")

# Filter to relevant columns
violence_df <- fsw_data_all %>%
  filter(use_exposed == "yes") %>%
  select(study, exposure_tf_bin, outcome, exposure_type, exposure_definition_short, perpetrator, exposed_num, exposed_perc, sample_size, hiv_num)

# Define the types of violence and their corresponding sheet names
violence_types <- c("Physical violence", "Sexual violence", "Physical or sexual", "Other violence")

# Initialize an empty list to store the results for each type of violence
all_violence_studies_list <- list()

# Loop over each type of violence
for (violence in violence_types) {
  # Filter for "Ever" and "Recent" exposures for the current type of violence
  ever_stud <- violence_df %>%
    filter(exposure_type == violence, exposure_tf_bin == "Ever") %>%
    select(study)
  recent_stud <- violence_df %>%
    filter(exposure_type == violence, exposure_tf_bin == "Recent") %>%
    select(study)
  
  # Combine and deduplicate studies
  combined_studies <- bind_rows(ever_stud, recent_stud) %>%
    arrange(study) %>%
    distinct(study, .keep_all = TRUE)
  
  # Store the result in the list
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

# keep relevant rows
for (df_name in dfs_studies) {
  df <- get(df_name)
  # Keep only the specified columns
  df <- df %>%
    select(study, exposure_definition_short, perpetrator, exposed_num, exposed_perc, sample_size)
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
  df$sample_size <- as.numeric(df$sample_size)
  df <- df[!is.na(df$exposed_num) & !is.na(df$sample_size) & df$sample_size > 0 & df$exposed_num >= 0, ]

  total_exposed_num <- sum(df$exposed_num)
  total_sample_size <- sum(df$sample_size)
  exposed_percentage <- (total_exposed_num / total_sample_size) * 100
  num_studies <- length(unique(df$study))
  num_estimates <- nrow(df)

  if (nrow(df) >= 2) {
    escalc_res <- escalc(
      measure = "PLO",
      xi = as.numeric(df$exposed_num),
      ni = as.numeric(df$sample_size)
    )
    escalc_res$effect_num <- seq_len(nrow(escalc_res))
    rma_res <- rma.mv(
      yi, vi,
      random = ~ 1 | study/effect_num,
      data = cbind(escalc_res, study = df$study)
    )
    pred <- predict(rma_res, transf=transf.ilogit)
    pooled_prev <- pred$pred
    pooled_prev_ci <- c(pred$ci.lb, pred$ci.ub)
    # I2 calculation
    var_comp <- as.numeric(rma_res$sigma2)
    total_var <- sum(var_comp) + mean(rma_res$vi)
    I2 <- 100 * var_comp[1] / total_var
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
    unadj_effect = unadj_est,
    unadj_lower = un_lower,
    unadj_upper = un_upper,
    adj_effect = adj_est,
    adj_lower = adj_lower,
    adj_upper = adj_upper,
    exposure_type = exposure_type
  ) %>%
  mutate(
    exposed_num = as.numeric(exposed_num),
    exposed_perc = as.numeric(exposed_perc),
    author_year = paste0(author, " (", year, ")"),
    exposed = ifelse(!is.na(exposed_num) & !is.na(exposed_perc),
                     paste0(exposed_num, " (", round(exposed_perc * 100, 0), "%)"),
                     "NR"),
    perpetrator = ifelse(is.na(perpetrator), "Any perpetrator", perpetrator),
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
    `Violence type` = exposure_type
  )

violence_sheets <- split(formatted_data, formatted_data$`Violence type`)
names(violence_sheets) <- names(violence_sheets) %>%
  str_replace_all("[\\[\\]:*?/\\\\]", "_") 

write_xlsx(violence_sheets, "Formatted Violence Data by Type.xlsx")