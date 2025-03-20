# Load required packages
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
 
# Combine all dataframes in the list dfs into a single dataframe
all_studies <- bind_rows(lapply(dfs, function(var) {
  df <- get(var)
  
  # Standardize column types (e.g., convert all problematic columns to character)
  df <- df %>%
    mutate(
      fsw_prop = as.character(fsw_prop),        # Ensure fsw_prop is character
      sample_size = as.character(sample_size),  # Standardize sample_size
      hiv_num = as.character(hiv_num),          # Standardize hiv_num
      hiv_perc = as.character(hiv_perc),        # Standardize hiv_perc
      exposed_num = as.character(exposed_num),  # Standardize exposed_num
      exposed_perc = as.character(exposed_perc) # Standardize exposed_perc
    )
  
  return(df)
}))

# Sort, group by study, and keep the first row of each group
all_studies <- all_studies %>%
  arrange(study) %>%
  group_by(study) %>%
  mutate(sequence = row_number()) %>%
  ungroup() %>%
  filter(sequence == 1)

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

# Update the all_studies dataframe
all_studies <- all_studies %>%
  # Select and reorder the columns
  select(
    study,
    location,
    design,
    pub_type,
    sample_size,
    exposed_num,
    exposed_perc,
    hiv_num,
    hiv_perc,
    recruitment,
    who_region,
    rob_score
  ) %>%
  
  # Rename the columns
  rename(
    "study" = study,
    "Location (city, country)" = location,
    "Study design" = design,
    "Publication type" = pub_type,
    "Sample size" = sample_size,
    "Exposed (n)" = exposed_num,
    "Exposed (%)" = exposed_perc,
    "HIV (n)" = hiv_num,
    "HIV (%)" = hiv_perc,
    "Recruitment strategy" = recruitment,
    "WHO region" = who_region,
    "ROB score" = rob_score
  )
