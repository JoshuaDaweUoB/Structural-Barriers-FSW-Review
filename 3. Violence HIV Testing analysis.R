# load packages 
pacman::p_load("meta", "metafor", "readxl", "tidyverse", "kableExtra", "robumeta", "clubSandwich") 

# set working directory
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence")

# settings
settings.meta(CIbracket = "(") 
settings.meta(CIseparator = "-") 

# columns
leftcols_recent <- c("study", "study_num", "effect_num", "exposure_definition_short", "exposure_time_frame", "outcome_time_frame", "perpetrator", "country")
leftlabs_recent <- c("Study", "Study number", "Effect number", "Exposure definition", "Exposure time frame", "Testing timeframe", "Perpetrator", "Country")
leftcols_lifetime <- c("study", "study_num", "effect_num", "exposure_definition_short",  "outcome_time_frame", "perpetrator", "country")
leftlabs_lifetime <- c("Study", "Study number", "Effect number", "Exposure definition", "Testing timeframe", "Perpetrator", "Country")
rightcols <- c("effect", "ci")
rightlabs = c("Estimate", "95% CI")

#### Multilevel random effects model with constant sampling correlation ####

# Constant sampling correlation
rho <- 0.6

# Define the different analyses
analyses <- c("unadj", "adj", "best")

# Define the corresponding variable names for each analysis
var_names <- list(
  unadj = list(est = "unadj_est_ln", var = "unadj_var_ln", lower = "un_lower_ln", upper = "un_upper_ln"),
  adj = list(est = "adj_est_ln", var = "adj_var_ln", lower = "adj_lower_ln", upper = "adj_upper_ln"),
  best = list(est = "effect_best_ln", var = "effect_best_var_ln", lower = "effect_best_lower_ln", upper = "effect_best_upper_ln")
)

## recent violence data 

# Define the corresponding plot filenames for each analysis
plot_filenames <- list(
  unadj = "Plots/recent_unadj_test.png",
  adj = "Plots/recent_adj_test.png",
  best = "Plots/recent_best_test.png"
)

# Function to perform the analysis and create forest plots
perform_analysis <- function(df, analysis) {
 
  # Filter the dataframe
  filtered_df <- df %>% filter(exposure_tf_bin == "Recent", use == "yes", !is.na(.[[var_names[[analysis]]$est]]))
  
  # Create study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)

  # Create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # Fit a multilevel random effects model using `rma.mv` from metafor
  result <- rma.mv(filtered_df[[var_names[[analysis]]$est]], 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE)       
  
  print(result)
  print(exp(coef(result)))
  
  result2 <- metagen(TE = filtered_df[[var_names[[analysis]]$est]],
                     lower = filtered_df[[var_names[[analysis]]$lower]],
                     upper = filtered_df[[var_names[[analysis]]$upper]],
                     studlab = filtered_df$study,
                     data = filtered_df,
                     sm = "OR",
                     method.tau = "REML",
                     common = FALSE,
                     random = TRUE, 
                     backtransf = TRUE,
                     text.random = "Overall")
  
  print(summary(result2))
  
  result2$TE.random <- result$b
  result2$lower.random <- result$ci.lb
  result2$upper.random <- result$ci.ub
  
  filename <- plot_filenames[[analysis]]
  png(filename = filename, width = 45, height = 14, units = "cm", res = 600)
  
  forest(result2,
         sortvar = filtered_df$study,
         xlim = c(0.2, 4),             
         leftcols = leftcols_recent, 
         leftlabs = leftlabs_recent,
         rightcols = rightcols,
         rightlabs = rightlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  dev.off()
}

# Loop over each analysis to perform the analysis and create forest plots
for (analysis in analyses) {
  perform_analysis(fsw_data_test, analysis)
}

## lifetime violence data

# Define the corresponding plot filenames for each analysis
plot_filenames <- list(
  unadj = "Plots/ever_unadj_test.png",
  adj = "Plots/ever_adj_test.png",
  best = "Plots/ever_best_test.png"
)

# Function to perform the analysis and create forest plots
perform_analysis <- function(df, analysis, exposure) {
  # Filter the dataframe
  filtered_df <- df %>% filter(exposure_tf_bin == exposure, use == "yes", !is.na(.[[var_names[[analysis]]$est]]))
  
  # Create study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)

  # Create a covariance matrix assuming constant sampling correlation
  V_mat <- impute_covariance_matrix(filtered_df[[var_names[[analysis]]$var]],
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # Fit a multilevel random effects model using `rma.mv` from metafor
  result <- rma.mv(filtered_df[[var_names[[analysis]]$est]], 
                   V = V_mat, 
                   random = ~ 1 | study_num / effect_num,
                   data = filtered_df,   
                   sparse = TRUE)       
  
  print(result)
  print(exp(coef(result)))
  
  result2 <- metagen(TE = filtered_df[[var_names[[analysis]]$est]],
                     lower = filtered_df[[var_names[[analysis]]$lower]],
                     upper = filtered_df[[var_names[[analysis]]$upper]],
                     studlab = filtered_df$study,
                     data = filtered_df,
                     sm = "OR",
                     method.tau = "REML",
                     common = FALSE,
                     random = TRUE, 
                     backtransf = TRUE,
                     text.random = "Overall")
  
  print(summary(result2))
  
  result2$TE.random <- result$b
  result2$lower.random <- result$ci.lb
  result2$upper.random <- result$ci.ub
  
  filename <- plot_filenames[[analysis]]
  png(filename = filename, width = 45, height = 14, units = "cm", res = 600)
  
  forest(result2,
         sortvar = filtered_df$study,
         xlim = c(0.2, 4),
         leftcols = leftcols_lifetime,
         leftlabs = leftlabs_lifetime,
         rightcols = rightcols,
         rightlabs = rightlabs,
         pooled.totals = TRUE,
         xintercept = 1,
         addrow.overall = TRUE,
         overall.hetstat = TRUE,
         overall = TRUE,
         labeltext = TRUE,
         col.subgroup = "black")
  
  dev.off()
}

# Loop over each analysis to perform the analysis and create forest plots for lifetime violence
for (analysis in analyses) {
  perform_analysis(fsw_data_test, analysis, "Ever")
}

#### overall forest plot ####

## load dataframe

summary_violence_test <- read_excel("Violence estimates.xlsx", "Testing") 

## forest plot

summary_hiv_violence_test  <- metagen(TE = effect_ln,
                                      lower = lower_ln,
                                      upper = upper_ln,
                                      studlab = name,
                                      data = summary_violence_test,
                                      sm = "OR",
                                      method.tau = "DL",
                                      comb.fixed = FALSE,
                                      comb.random = FALSE, 
                                      backtransf = TRUE,
                                      byvar = name,
                                      text.random = "Overall")

summary(summary_hiv_violence_test) 

filename <- paste0("Plots/overall plots/violence_test_overall.png")
png(filename = filename, width = 25, height = 14, units = "cm", res = 600)

summary_hiv_violence_test <- forest(summary_hiv_violence_test, 
                                    sortvar = name,
                                    xlim=c(0.2, 4),             
                                    leftcols = c("Adjust", "studies", "estimates"), 
                                    leftlabs = c("Model type", "Studies", "Estimates"),
                                    rightcols = c("or_95", "i2"), 
                                    rightlabs = c("OR (95% CI)", "I²"), 
                                    pooled.totals = F,
                                    xintercept=1,
                                    addrow.overall = T,
                                    test.subgroup = F,
                                    overall.hetstat = F,
                                    overall = F,
                                    labeltext = TRUE,
                                    col.subgroup = "black",
                                    print.subgroup.name = FALSE) 
dev.off()

### subgroup analysis

## WHO region

# africa 

filtered_df <- fsw_data_test %>% filter(exposure_tf_bin == "Ever")
filtered_df <- filtered_df %>% filter(effect_best != "NA")
filtered_df <- filtered_df %>% filter(who_region == "African Region")

# create a covariance matrix assuming constant sampling correlation
V_mat <- impute_covariance_matrix(filtered_df$effect_best_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)

# fit a multilevel random effects model using `rma.mv` from metafor
result <- rma.mv(effect_best_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)       

result 
exp(coef(result))


result2 <- metagen(TE = effect_best_ln,
                   lower = effect_best_lower_ln,
                   upper = effect_best_upper_ln,
                   studlab = study,
                   data = filtered_df,
                   sm = "OR",
                   method.tau = "REML",
                   common = FALSE,
                   random = TRUE, 
                   backtransf = TRUE,
                   text.random = "Overall")

print(summary(result2))

result2$TE.random <- result$b
result2$lower.random <- result$ci.lb
result2$upper.random <- result$ci.ub

filename <- paste0("Plots/subgroups/ever_best_test_africa.png")
png(filename = filename, width = 38, height = 14, units = "cm", res = 600)

forest(result2,
       sortvar = study,
       xlim = c(0.2, 4),             
       leftcols = leftcols_recent, 
       leftlabs = leftlabs_recent,
       rightcols = rightcols,
       rightlabs = rightlabs,
       pooled.totals = TRUE,
       xintercept = 1,
       addrow.overall = TRUE,
       overall.hetstat = TRUE,
       overall = TRUE,
       labeltext = TRUE,
       col.subgroup = "black")

dev.off()  

# asia

filtered_df <- fsw_data_test %>% filter(exposure_tf_bin == "Ever")
filtered_df <- filtered_df %>% filter(effect_best != "NA")
filtered_df <- filtered_df %>% filter(who_region == "South-East Asian Region")

# create a covariance matrix assuming constant sampling correlation
V_mat <- impute_covariance_matrix(filtered_df$effect_best_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)

# fit a multilevel random effects model using `rma.mv` from metafor
result <- rma.mv(effect_best_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)       

result 
exp(coef(result))


result2 <- metagen(TE = effect_best_ln,
                   lower = effect_best_lower_ln,
                   upper = effect_best_upper_ln,
                   studlab = study,
                   data = filtered_df,
                   sm = "OR",
                   method.tau = "REML",
                   common = FALSE,
                   random = TRUE, 
                   backtransf = TRUE,
                   text.random = "Overall")

print(summary(result2))

result2$TE.random <- result$b
result2$lower.random <- result$ci.lb
result2$upper.random <- result$ci.ub

filename <- paste0("Plots/subgroups/ever_best_test_asia.png")
png(filename = filename, width = 38, height = 14, units = "cm", res = 600)

forest(result2,
       sortvar = study,
       xlim = c(0.2, 4),             
       leftcols = leftcols_recent, 
       leftlabs = leftlabs_recent,
       rightcols = rightcols,
       rightlabs = rightlabs,
       pooled.totals = TRUE,
       xintercept = 1,
       addrow.overall = TRUE,
       overall.hetstat = TRUE,
       overall = TRUE,
       labeltext = TRUE,
       col.subgroup = "black")

dev.off()  

## LIC


## recruitment method

##perpetrator

## year

## ROB
















