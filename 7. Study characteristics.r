# Load required packages
pacman::p_load("readxl", "writexl", "tidyverse")

# Combine all dataframes into one
combine_and_process_data <- function(dfs) {
  # Combine all dataframes
  all_studies <- bind_rows(lapply(dfs, function(var) {
    df <- get(var)
    
    # Ensure consistent columns and data types
    required_columns <- c(
      "study", "title", "pub_type", "design", "study_setting", "sample_size",
      "recruitment", "country", "cities", "who_region", "rob_score",
      "sw_time_frame", "hiv_num", "hiv_perc", "fsw_prop", "exposed_num", "exposed_perc"
    )
    
    df <- df %>%
      mutate(across(all_of(required_columns), ~ ifelse(is.null(.), NA, .))) %>%
      mutate(
        sample_size = as.character(sample_size),
        hiv_num = as.character(hiv_num),
        hiv_perc = as.character(hiv_perc),
        fsw_prop = as.character(fsw_prop),
        exposed_num = as.character(exposed_num),
        exposed_perc = as.character(exposed_perc)
      )
    
    return(df)
  }))
  
  # Deduplicate based on relevant columns
  all_studies <- all_studies %>%
    distinct(
      study, design, pub_type, sample_size, hiv_num, hiv_perc,
      recruitment, who_region, rob_score, cities, country, .keep_all = TRUE
    )
  
  # Create a 'location' column
  all_studies <- all_studies %>%
    mutate(location = paste0(cities, ", ", country))
  
  # Reorder and rename columns
  all_studies <- all_studies %>%
    select(
      study,
      location,
      design,
      pub_type,
      sample_size,
      hiv_num,
      hiv_perc,
      recruitment,
      who_region,
      rob_score
    ) %>%
    rename(
      " " = study,
      "Location(city, country)" = location,
      "Study design" = design,
      "Publication type" = pub_type,
      "Sample size" = sample_size,
      "HIV+" = hiv_num,
      "HIV %" = hiv_perc,
      "Recruitment strategy" = recruitment,
      "WHO region" = who_region,
      "ROB score" = rob_score
    )
  
  return(all_studies)
}

# Define the list of dataframes and output file
dfs <- c("fsw_data_test", "fsw_data_pv_ever", "fsw_data_pv_recent", 
         "fsw_data_sv_ever", "fsw_data_sv_recent", "fsw_data_psv_ever", 
         "fsw_data_psv_recent", "fsw_data_art", "fsw_data_vs")
output_file <- "All Studies.xlsx"

# Process the data
all_studies <- combine_and_process_data(dfs)

# Save the final dataframe to an Excel file
write_xlsx(list(all_studies = all_studies), path = output_file)

# Print the number of rows in the final dataframe
print(paste("Number of rows in the final dataframe:", nrow(all_studies)))

# Print the final dataframe
print(all_studies)