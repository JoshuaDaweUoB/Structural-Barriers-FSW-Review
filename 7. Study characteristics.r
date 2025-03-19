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

