# load packages
pacman::p_load("readxl", "writexl", "tidyverse")

# Function to process and save study data
process_and_save_study_data <- function(dfs, output_file) {
  # Create a list to store processed dataframes
  processed_dataframes <- list()
  
  # Loop through each dataframe in the list
  for (var in dfs) {
    # Get the dataframe by name
    df <- get(var)
    
    # Process the dataframe
    df <- df %>%
      arrange(study) %>% # Sort by the column 'study'
      group_by(study) %>% # Group by 'study'
      mutate(sequence = row_number()) %>% # Create a sequence within each group
      ungroup() %>%
      filter(sequence == 1) %>% # Keep only rows where sequence equals 1
      select(study, title, pub_type, design, study_setting, sample_size, recruitment, 
             country, cities, who_region, rob_score, sw_time_frame, hiv_num, hiv_perc) # Keep specified columns
    
    # Sort the dataframe by 'study' before saving
    df <- df %>% arrange(study)

    # Add the processed dataframe to the list
    # Use the original name (without the "_characteristics" suffix) as the sheet name
    processed_dataframes[[var]] <- df
  }
  
  # Save all processed dataframes into an Excel file with sheet names
  write_xlsx(processed_dataframes, path = output_file)
}

# Define the output file name
output_file <- "Study characteristics.xlsx"

# Call the function
process_and_save_study_data(dfs, output_file)

# Combine all processed dataframes into a single dataframe
all_studies <- bind_rows(lapply(dfs, function(var) {
  get(var) %>%
    mutate(
      sample_size = as.character(sample_size), # Ensure sample_size is of the same type
      hiv_num = as.character(hiv_num),         # Ensure hiv_num is of the same type
      hiv_perc = as.character(hiv_perc)        # Ensure hiv_perc is of the same type
    ) %>%
    arrange(study) %>% # Sort by the column 'study'
    group_by(study) %>%
    mutate(sequence = row_number()) %>%
    ungroup() %>%
    filter(sequence == 1) %>%
    select(study, title, pub_type, design, study_setting, sample_size, recruitment, 
           country, cities, who_region, rob_score, sw_time_frame, hiv_num, hiv_perc)
}), .id = "source_dataframe")

# Sort the combined dataframe by the 'study' column
all_studies <- all_studies %>% arrange(study)

# View the combined dataframe
print(all_studies)

# Save the combined dataframe to a separate Excel file
write_xlsx(list(all_studies = all_studies), path = "All Studies.xlsx")