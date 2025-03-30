# load packages 
pacman::p_load("meta", "metafor", "readxl", "openxlsx", "tidyverse", "kableExtra", "robumeta", "clubSandwich") 

# set working directory
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence")

# Call the function for "recently exposed to physical violence"
process_and_plot(
  data = fsw_data_pv_recent,
  data_name = "fsw_data_pv_recent",
  output_plot_filename = "Plots/subgroups/recent_pv_subgroup.png"
)

# Call the function for "ever exposed to physical violence"
process_and_plot(
  data = fsw_data_pv_ever,
  data_name = "fsw_data_pv_ever",
  output_plot_filename = "Plots/subgroups/ever_pv_subgroup.png"
)

# Call the function for "recently exposed to sexual violence"
process_and_plot(
  data = fsw_data_sv_recent,
  data_name = "fsw_data_sv_recent",
  output_plot_filename = "Plots/subgroups/recent_sv_subgroup.png"
)

# Call the function for "ever exposed to sexual violence"
process_and_plot(
  data = fsw_data_sv_ever,
  data_name = "fsw_data_sv_ever",
  output_plot_filename = "Plots/subgroups/ever_sv_subgroup.png"
)

# Call the function for "recently exposed to physical and/or sexual violence"
process_and_plot(
  data = fsw_data_psv_recent,
  data_name = "fsw_data_psv_recent",
  output_plot_filename = "Plots/subgroups/recent_psv_subgroup.png"
)

# Call the function for "ever exposed to physical and/or sexual violence"
process_and_plot(
  data = fsw_data_psv_ever,
  data_name = "fsw_data_psv_ever",
  output_plot_filename = "Plots/subgroups/ever_psv_subgroup.png"
)


