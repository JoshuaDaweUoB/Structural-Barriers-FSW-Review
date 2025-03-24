
# load packages 
pacman::p_load("meta", "metafor", "readxl", "openxlsx", "tidyverse", "kableExtra", "robumeta", "clubSandwich") 

# set working directory
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence")

# load data
fsw_data_overall <- read_excel("Violence estimates.xlsx", "All studies") 

# settings
settings.meta(CIbracket = "(") 
settings.meta(CIseparator = "-") 

fsw_data_overall$i2 <- round(fsw_data_overall$i2 * 100, 0)  # Convert i2 to percentage and round to no decimals
fsw_data_overall$i2 <- paste0(fsw_data_overall$i2, "%")  # Add percentage sign to i2 values

# Perform meta-analysis
  summary_hiv_violence <- metagen(TE = effect_ln,
                                  lower = lower_ln,
                                  upper = upper_ln,
                                  data = fsw_data_overall,
                                  sm = "OR",
                                  method.tau = "DL",
                                  comb.fixed = FALSE,
                                  comb.random = FALSE, 
                                  backtransf = TRUE,
                                  byvar = outcome,
                                  text.random = "Overall")
  
  # Print summary
  print(summary(summary_hiv_violence))
  
  # Save forest plot
  png(filename = "Plots/overall plots/violence_all_analyses.png", width = 30, height = 28, units = "cm", res = 600)
  forest(summary_hiv_violence, 
         sortvar = outcome,
         xlim = c(0.2, 4),             
         leftcols = c("violence", "name", "adjust"), 
         leftlabs = c("Violence exposure", "Exposure time frame", "Adjustment"),
         rightcols = c("or_95", "i2", "studies", "estimates"), 
         rightlabs = c("OR (95% CI)", "I²", "Nb studies", "Nb estimates"), 
         pooled.totals = FALSE,
         xintercept = 1,
         addrow.overall = TRUE,
         test.subgroup = FALSE,
         overall.hetstat = FALSE,
         overall = FALSE,
         labeltext = TRUE,
         col.subgroup = "black",
         print.subgroup.name = FALSE)
  dev.off()

