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
    location,
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
    "Location (city, country)" = location,
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

View(all_studies)
# save the all_studies dataframe to an Excel file
write_xlsx(all_studies, "All studies.xlsx")

# characteristics of included studies

# Descriptive table for "WHO region"
who_region_table <- all_studies %>%
  count(`WHO region`) %>%  # Use backticks for column names with spaces
  mutate(percentage = n / sum(n) * 100)  # Calculate percentages

# Descriptive table for "WHO region"
study_design_table <- all_studies %>%
  count(`Study design`) %>%  # Use backticks for column names with spaces
  mutate(percentage = n / sum(n) * 100)  # Calculate percentages

# Descriptive table for "Publication type"
pub_type_table <- all_studies %>%
  count(`Publication type`) %>%  # Use backticks for column names with spaces
  mutate(percentage = n / sum(n) * 100)  # Calculate percentages

# Print the tables
print(who_region_table)
print(pub_type_table)
print(study_design_table)

# Calculate total Sample size and total HIV (n)
totals_table <- all_studies %>%
  summarise(
    total_sample_size = sum(as.numeric(`Sample size`), na.rm = TRUE),  # Convert Sample size to numeric and sum
    total_hiv_n = sum(as.numeric(`HIV (n)`), na.rm = TRUE)             # Convert HIV (n) to numeric and sum
  )

# Print the totals table
print(totals_table)

## total studies for violence types

# Define the types of violence and their corresponding sheet names
violence_types <- c("Physical violence", "Sexual violence", "Physical or sexual")

# Initialize an empty list to store the results for each type of violence
all_violence_studies_list <- list()

# Loop over each type of violence
for (violence in violence_types) {
  # Read the "Ever" and "Recent" sheets for the current type of violence
  ever_stud <- read_excel("All violence studies.xlsx", paste0(violence, " - Ever")) %>%
    select(study)
  recent_stud <- read_excel("All violence studies.xlsx", paste0(violence, " - Recent")) %>%
    select(study)
  
  # Combine the "Ever" and "Recent" data
  combined_studies <- bind_rows(ever_stud, recent_stud) %>%
    arrange(study) %>% # Sort by the column 'study'
    group_by(study) %>% # Group by 'study'
    mutate(sequence = row_number()) %>% # Create a sequence within each group
    ungroup() %>%
    filter(sequence == 1) # Keep only rows where sequence equals 1
  
  # Store the result in the list with the type of violence as the key
  all_violence_studies_list[[violence]] <- combined_studies
}

# Access the results for each type of violence
physical_violence_studies <- all_violence_studies_list[["Physical violence"]]
sexual_violence_studies <- all_violence_studies_list[["Sexual violence"]]
physical_or_sexual_violence_studies <- all_violence_studies_list[["Physical or sexual"]]
View(physical_violence_studies)
View(sexual_violence_studies)
View(physical_or_sexual_violence_studies)

# recent and lifetime violence

# Load physical violence data
fsw_data_pv_ever_studies <- read_excel("All violence studies.xlsx", "Physical violence - Ever") 
fsw_data_pv_recent_studies <- read_excel("All violence studies.xlsx", "Physical violence - Recent") 

# Load sexual violence data
fsw_data_sv_ever_studies<- read_excel("All violence studies.xlsx", "Sexual violence - Ever") 
fsw_data_sv_recent_studies <- read_excel("All violence studies.xlsx", "Sexual violence - Recent") 

# Load physical & sexual violence data
fsw_data_psv_ever_studies <- read_excel("All violence studies.xlsx", "Physical or sexual - Ever") 
fsw_data_psv_recent_studies <- read_excel("All violence studies.xlsx", "Physical or sexual - Recent") 

# create list of dataframes
dfs_studies <- c("fsw_data_pv_ever_studies", "fsw_data_pv_recent_studies", "fsw_data_sv_ever_studies", "fsw_data_sv_recent_studies", "fsw_data_psv_ever_studies", "fsw_data_psv_recent_studies")

# Loop over the list of dataframe names
for (df_name in dfs_studies) {
  # Access the dataframe by its name
  df <- get(df_name)
  
  # Convert exposed_num and sample_size to numeric, treating "NA" strings as NA
  df <- df %>%
    mutate(
      exposed_num = as.numeric(ifelse(exposed_num == "NA", NA, exposed_num)),
      sample_size = as.numeric(ifelse(sample_size == "NA", NA, sample_size))
    ) %>%
    # Filter rows with no missing values for exposed_num and sample_size
    filter(!is.na(exposed_num) & !is.na(sample_size)) %>%
    # Create a new variable exposed_perc2
    mutate(
      exposed_perc2 = exposed_num / sample_size * 100  # Calculate percentage
    ) %>%
    arrange(study, desc(exposed_perc)) %>% # Sort by 'study' and then by 'exposed_perc' in descending order
    group_by(study) %>% # Group by 'study'
    mutate(sequence = row_number()) %>% # Create a sequence within each group
    ungroup() %>%
    filter(sequence == 1) %>% # Keep only rows where sequence equals 1
    select(study, exposed_num, exposed_perc, sample_size, exposed_perc2) # Keep specified columns
  
  # Assign the processed dataframe back to its original name
  assign(df_name, df)
}

# Initialize an empty dataframe to store the summary results
summary_table <- data.frame(
  dataframe = character(),
  total_exposed_num = numeric(),
  total_sample_size = numeric(),
  exposed_percentage = numeric(),
  stringsAsFactors = FALSE
)

# Loop over the list of dataframe names
for (df_name in dfs_studies) {
  # Access the dataframe by its name
  df <- get(df_name)
  
  # Calculate the totals
  total_exposed_num <- sum(df$exposed_num, na.rm = TRUE)  # Total exposed_num
  total_sample_size <- sum(df$sample_size, na.rm = TRUE)  # Total sample_size
  exposed_percentage <- (total_exposed_num / total_sample_size) * 100  # Percentage
  
  # Add the results to the summary table
  summary_table <- rbind(
    summary_table,
    data.frame(
      dataframe = df_name,
      total_exposed_num = total_exposed_num,
      total_sample_size = total_sample_size,
      exposed_percentage = exposed_percentage
    )
  )
}

# Print the summary table
print(summary_table)
write_xlsx(summary_table, "Summary Table.xlsx")


