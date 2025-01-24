# load packages 
pacman::p_load("meta", "metafor", "readxl", "tidyverse", "kableExtra", "robumeta", "clubSandwich") 

# set working directory
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence")

# settings
settings.meta(CIbracket = "(") 
settings.meta(CIseparator = "-") 

# columns
leftcols_recent <- c("study", "exposure_definition_short", "exposure_time_frame", "perpetrator", "country")
leftlabs_recent <- c("Study", "Exposure definition", "Exposure time frame", "Perpetrator", "Country")
leftcols_lifetime <- c("study", "exposure_definition_short", "perpetrator", "country")
leftlabs_lifetime <- c("Study", "Exposure definition", "Perpetrator", "Country")
rightcols <- c("effect", "ci")
rightlabs = c("Estimate", "95% CI")

#### approach one: meta analysis multi-level model ####

# Function to perform the analysis and generate forest plots
perform_analysis <- function(data, outcome_filter, estimate_filter, var_column, plot_filename, adj = FALSE) {
  filtered_df <- data %>% filter(outcome == outcome_filter)
  filtered_df <- filtered_df %>% filter(!is.na(!!sym(estimate_filter)))
  
  if (nrow(filtered_df) == 0) {
    message(paste("No data available for", estimate_filter, "in", plot_filename))
    return(NULL)
  }
  
  if (adj) {
    # Fit the multilevel random effects model for adjusted analyses
    result <- rma.mv(
      yi = !!sym(estimate_filter),                  
      V = filtered_df[[var_column]],                   
      slab = study,                      
      data = filtered_df,              
      random = ~ 1 | study/effect_num,     
      test = "t",                        
      method = "REML"                 
    )
  } else {
    # Create a covariance matrix assuming constant sampling correlation
    V_mat <- impute_covariance_matrix(filtered_df[[var_column]],
                                      cluster = filtered_df$study_num,
                                      r = rho,
                                      smooth_vi = TRUE)
    
    # Fit the multilevel random effects model for unadjusted analyses
    result <- rma.mv(!!sym(estimate_filter), 
                     V = V_mat, 
                     random = ~ 1 | study_num / effect_num,
                     data = filtered_df,   
                     control = list(rel.tol = 1e-8),
                     sparse = TRUE)
  }
  
  # Use `result` to update the `metagen` object
  result2 <- metagen(TE = filtered_df[[estimate_filter]],
                     lower = filtered_df$effect_best_lower_ln,
                     upper = filtered_df$effect_best_upper_ln,
                     studlab = filtered_df$study,
                     data = filtered_df,
                     sm = "OR",
                     method.tau = "REML",
                     common = FALSE,
                     random = TRUE, 
                     backtransf = TRUE,
                     text.random = "Overall")
  
  # Update result2 with values from the multilevel model
  result2$TE.random <- result$b
  result2$lower.random <- result$ci.lb
  result2$upper.random <- result$ci.ub
  
  # Save the forest plot as a PNG file
  png(filename = plot_filename, width = 38, height = 18, units = "cm", res = 600)
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

# List of datasets and corresponding parameters
datasets <- list(fsw_data_pv_recent, fsw_data_sv_recent, fsw_data_psv_recent, fsw_data_pv_ever)
dataset_names <- c("fsw_data_pv_recent", "fsw_data_sv_recent", "fsw_data_psv_recent", "fsw_data_pv_ever")
outcome_filter <- "HIV prevalence"
estimate_filters <- c("unadj_est_ln", "effect_best_ln", "adj_est_ln")
var_columns <- c("unadj_var_ln", "effect_best_var_ln", "adj_var_ln")
plot_filenames <- c("Plots/psv_recent_unadj_csc.png", "Plots/psv_recent_best_csc.png", "Plots/psv_recent_adj_csc.png", "Plots/pv_ever_unadj_csc.png")

# Loop through each dataset and perform the analysis
for (i in 1:length(datasets)) {
  for (j in 1:length(estimate_filters)) {
    perform_analysis(datasets[[i]], outcome_filter, estimate_filters[j], var_columns[j], plot_filenames[j], adj = (estimate_filters[j] == "adj_est_ln"))
  }
}

#### approach two: multilevel random effects model with constant sampling correlation ####

# constant sampling correlation
rho <- 0.6

{ 
### physical violence ### 
  
## unadjusted analyses ##
  
{ 
# recent exposure #
    
filtered_df <- fsw_data_pv_recent %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(unadj_est != "NA")
    
# create a covariance matrix assuming constant sampling correlation
V_mat <- impute_covariance_matrix(filtered_df$unadj_var_ln,
                                      cluster = filtered_df$study_num,
                                      r = rho,
                                      smooth_vi = TRUE)
    
# fit a multilevel random effects model using `rma.mv` from metafor
result <- rma.mv(unadj_est_ln, 
                            V = V_mat, 
                            random = ~ 1 | study_num / effect_num,
                            data = filtered_df,   
                            sparse = TRUE)       
    
result 
exp(coef(result))

    
result2 <- metagen(TE = unadj_est_ln,
                       lower = un_lower_ln,
                       upper = un_upper_ln,
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
    
filename <- paste0("Plots/pv_recent_unadj_csc.png")
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
    
# lifetime exposure #
    
filtered_df <- fsw_data_pv_ever %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(unadj_est != "NA")
    
# create a covariance matrix assuming constant sampling correlation
V_mat <- impute_covariance_matrix(filtered_df$unadj_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)

# fit a multilevel random effects model using `rma.mv` from metafor
result <- rma.mv(unadj_est_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)       

result 
exp(coef(result))
    
result2 <- metagen(TE = unadj_est_ln,
                       lower = un_lower_ln,
                       upper = un_upper_ln,
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
    
filename <- paste0("Plots/pv_ever_unadj_csc.png")
    png(filename = filename, width = 35, height = 18, units = "cm", res = 600)
    
forest(result2,
           sortvar = study,
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
  
## adjusted analyses ##
  
{
# recent exposure #
    
filtered_df <- fsw_data_pv_recent %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(adj_est != "NA")
    
# create a covariance matrix assuming constant sampling correlation
V_mat <- impute_covariance_matrix(filtered_df$adj_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)
    
# fit a multilevel random effects model using `rma.mv` from metafor
result <- rma.mv(adj_est_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)       
    
result 
exp(coef(result))
    
result2 <- metagen(TE = adj_est_ln,
                   lower = adj_lower_ln,
                   upper = adj_upper_ln,
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
    
filename <- paste0("Plots/pv_recent_adj_csc.png")
png(filename = filename, width = 38, height = 18, units = "cm", res = 600)
    
forest(result2,
           sortvar = study,
           xlim = c(0.2, 4),             
           leftcols = leftcols_recent, 
           leftlabs = leftlabs_recent,
           rightcols = rightcols,
           rightlabs = rightlabs,
           xintercept = 1,
           pooled.totals = TRUE,
           addrow.overall = TRUE,
           overall.hetstat = TRUE,
           overall = TRUE,
           labeltext = TRUE)
    
dev.off()
  
# lifetime exposure #
    
filtered_df <- fsw_data_pv_ever %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(adj_est != "NA")
    
# create a covariance matrix assuming constant sampling correlation
V_mat <- impute_covariance_matrix(filtered_df$unadj_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)

# fit a multilevel random effects model using `rma.mv` from metafor
result <- rma.mv(unadj_est_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)       

result 
exp(coef(result))
    
result2 <- metagen(TE = adj_est_ln,
                   lower = adj_lower_ln,
                   upper = adj_upper_ln,
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
    
filename <- paste0("Plots/pv_ever_adj_csc.png")
png(filename = filename, width = 38, height = 18, units = "cm", res = 600)
    
forest(result2,
           sortvar = study,
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
  
## best analyses ##
  
{
# recent exposure #
    
filtered_df <- fsw_data_pv_recent %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(effect_best != "NA")
    
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
    
filename <- paste0("Plots/pv_recent_best_csc.png")
png(filename = filename, width = 38, height = 18, units = "cm", res = 600)
    
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
    
# lifetime exposure #
    
filtered_df <- fsw_data_pv_ever %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(effect_best != "NA")
    
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
    
filename <- paste0("Plots/pv_ever_best_csc.png")
png(filename = filename, width = 38, height = 18, units = "cm", res = 600)
    
forest(result2,
           sortvar = study,
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
  
### sexual violence ###
  
## unadjusted analyses ##
  
{
# recent exposure #
filtered_df <- fsw_data_sv_recent %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(unadj_est != "NA")
    
# create a covariance matrix assuming constant sampling correlation
V_mat <- impute_covariance_matrix(filtered_df$unadj_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)

# fit a multilevel random effects model using `rma.mv` from metafor
result <- rma.mv(unadj_est_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)       
    
result 
exp(coef(result))
    
result2 <- metagen(TE = unadj_est_ln,
                       lower = un_lower_ln,
                       upper = un_upper_ln,
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
    
filename <- paste0("Plots/sv_recent_unadj_csc.png")
png(filename = filename, width = 38, height = 18, units = "cm", res = 600)
    
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
    
# lifetime exposure#
    
filtered_df <- fsw_data_sv_ever %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(unadj_est != "NA")
    
# create a covariance matrix assuming constant sampling correlation
V_mat <- impute_covariance_matrix(filtered_df$unadj_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)

# fit a multilevel random effects model using `rma.mv` from metafor
result <- rma.mv(unadj_est_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)       
    
result 
exp(coef(result))
    
    result2 <- metagen(TE = unadj_est_ln,
                       lower = un_lower_ln,
                       upper = un_upper_ln,
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
    
filename <- paste0("Plots/sv_ever_unadj_csc.png")
png(filename = filename, width = 38, height = 18, units = "cm", res = 600)
    
forest(result2,
           sortvar = study,
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
  
## adjusted analyses ##
  
{
# recent exposure # 
    
filtered_df <- fsw_data_sv_recent %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(adj_est != "NA")
    
# create a covariance matrix assuming constant sampling correlation
V_mat <- impute_covariance_matrix(filtered_df$adj_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)

# fit a multilevel random effects model using `rma.mv` from metafor
result <- rma.mv(adj_est_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)       
    
result 
exp(coef(result))
    
result2 <- metagen(TE = adj_est_ln,
                   lower = adj_lower_ln,
                   upper = adj_upper_ln,
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
    
filename <- paste0("Plots/sv_recent_adj_csc.png")
png(filename = filename, width = 38, height = 18, units = "cm", res = 600)
    
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
    
# lifetime exposure #
    
filtered_df <- fsw_data_sv_ever %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(adj_est != "NA")
    
# create a covariance matrix assuming constant sampling correlation
V_mat <- impute_covariance_matrix(filtered_df$adj_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)

# fit a multilevel random effects model using `rma.mv` from metafor
result <- rma.mv(adj_est_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)       
    
result 
exp(coef(result))
    
result2 <- metagen(TE = adj_est_ln,
                   lower = adj_lower_ln,
                   upper = adj_upper_ln,
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
    
filename <- paste0("Plots/sv_ever_adj_csc.png")
png(filename = filename, width = 38, height = 18, units = "cm", res = 600)
    
forest(result2,
           sortvar = study,
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
  
## best analyses ##
  
{
    
# recent exposure #
    
filtered_df <- fsw_data_sv_recent %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(effect_best != "NA")
    
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
    
filename <- paste0("Plots/sv_recent_best_csc.png")
png(filename = filename, width = 38, height = 18, units = "cm", res = 600)
    
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
    
# lifetime exposure #
    
filtered_df <- fsw_data_sv_ever %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(effect_best != "NA")
    
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
    
filename <- paste0("Plots/sv_ever_best_csc.png")
png(filename = filename, width = 38, height = 18, units = "cm", res = 600)
    
forest(result2,
           sortvar = study,
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
  
### physical and/or sexual violence ###
  
## unadjusted analyses ##
  
{
# recent exposure #
filtered_df <- fsw_data_psv_recent %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(unadj_est != "NA")
    
# create a covariance matrix assuming constant sampling correlation
V_mat <- impute_covariance_matrix(filtered_df$unadj_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)

# fit a multilevel random effects model using `rma.mv` from metafor
result <- rma.mv(unadj_est_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)       
    
result 
exp(coef(result))
    
result2 <- metagen(TE = unadj_est_ln,
                   lower = un_lower_ln,
                   upper = un_upper_ln,
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
    
filename <- paste0("Plots/psv_recent_unadj_csc.png")
png(filename = filename, width = 35, height = 18, units = "cm", res = 600)
    
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
    
# lifetime exposure#
    
filtered_df <- fsw_data_psv_ever %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(unadj_est != "NA")
    
# create a covariance matrix assuming constant sampling correlation
V_mat <- impute_covariance_matrix(filtered_df$unadj_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)

# fit a multilevel random effects model using `rma.mv` from metafor
result <- rma.mv(unadj_est_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)       
    
result 
exp(coef(result))
    
result2 <- metagen(TE = unadj_est_ln,
                       lower = un_lower_ln,
                       upper = un_upper_ln,
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
    
filename <- paste0("Plots/psv_ever_unadj_csc.png")
png(filename = filename, width = 35, height = 18, units = "cm", res = 600)
    
forest(result2,
           sortvar = study,
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
  
## adjusted analyses ##
  
{
# recent exposure # 
  
filtered_df <- fsw_data_psv_recent %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(adj_est != "NA")
    
# create a covariance matrix assuming constant sampling correlation
V_mat <- impute_covariance_matrix(filtered_df$adj_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)

# fit a multilevel random effects model using `rma.mv` from metafor
result <- rma.mv(adj_est_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)       
    
result 
exp(coef(result))
    
result2 <- metagen(TE = adj_est_ln,
                       lower = adj_lower_ln,
                       upper = adj_upper_ln,
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
    
filename <- paste0("Plots/psv_recent_adj_csc.png")
png(filename = filename, width = 35, height = 18, units = "cm", res = 600)
    
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
    
# lifetime exposure #
  
filtered_df <- fsw_data_psv_ever %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(adj_est != "NA")
    
# create a covariance matrix assuming constant sampling correlation
V_mat <- impute_covariance_matrix(filtered_df$adj_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)

# fit a multilevel random effects model using `rma.mv` from metafor
result <- rma.mv(adj_est_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)       
    
result 
exp(coef(result))
    
result2 <- metagen(TE = adj_est_ln,
                   lower = adj_lower_ln,
                   upper = adj_upper_ln,
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
    
filename <- paste0("Plots/psv_ever_adj_csc.png")
png(filename = filename, width = 35, height = 18, units = "cm", res = 600)
    
forest(result2,
       sortvar = study,
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
  
## best analyses ##
  
{
# recent exposure #
    
filtered_df <- fsw_data_psv_recent %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(effect_best != "NA")
    
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
    
filename <- paste0("Plots/psv_recent_best_csc.png")
png(filename = filename, width = 35, height = 18, units = "cm", res = 600)
    
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
    
# lifetime exposure #
    
filtered_df <- fsw_data_psv_ever %>% filter(outcome == "HIV prevalence")
filtered_df <- filtered_df %>% filter(effect_best != "NA")
    
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
    
filename <- paste0("Plots/psv_ever_best_csc.png")
png(filename = filename, width = 35, height = 18, units = "cm", res = 600)
    
forest(result2,
       sortvar = study,
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
  
} 



#### overall forest plots ####

### ever exposed ###

#### overall forest plots ####

### ever exposed ###

## load dataframe

summary_violence_ever <- read_excel("Violence estimates.xlsx", "HIV infection - ever") 

## forest plot

summary_hiv_violence_ever  <- metagen(TE = effect_ln_2,
                                      lower = lower_ln_2,
                                      upper = upper_ln_2,
                                      studlab = name,
                                      data = summary_violence_ever,
                                      sm = "OR",
                                      method.tau = "DL",
                                      comb.fixed = FALSE,
                                      comb.random = FALSE, 
                                      backtransf = TRUE,
                                      byvar = Adjust,
                                      text.random = "Overall")

summary(summary_hiv_violence_ever) 

filename <- paste0("Plots/overall plots/violence_ever_overall.png")
png(filename = filename, width = 25, height = 14, units = "cm", res = 600)

summary_hiv_violence_ever <- forest(summary_hiv_violence_ever, 
                                         sortvar = name,
                                         xlim=c(0.2, 4),             
                                         leftcols = c("name", "studies", "estimates"), 
                                         leftlabs = c("Pooled exposure", "Studies", "Estimates"),
                                         rightcols = c("or_95_2", "i2"), 
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

### recent exposed ###

## load dataframe

summary_violence_recent <- read_excel("Violence estimates.xlsx", "HIV infection - recent") 

## forest plot

summary_hiv_violence_recent  <- metagen(TE = effect_ln_2,
                                      lower = lower_ln_2,
                                      upper = upper_ln_2,
                                      studlab = name,
                                      data = summary_violence_recent,
                                      sm = "OR",
                                      method.tau = "DL",
                                      comb.fixed = FALSE,
                                      comb.random = FALSE, 
                                      backtransf = TRUE,
                                      byvar = Adjust,
                                      text.random = "Overall")

summary(summary_hiv_violence_recent) 

filename <- paste0("Plots/overall plots/violence_recent_overall.png")
png(filename = filename, width = 25, height = 14, units = "cm", res = 600)

summary_hiv_violence_recent <- forest(summary_hiv_violence_recent, 
                                    sortvar = name,
                                    xlim=c(0.2, 4),             
                                    leftcols = c("name", "studies", "estimates"), 
                                    leftlabs = c("Pooled exposure", "Studies", "Estimates"),
                                    rightcols = c("or_95_2", "i2"), 
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

#### subgroup analyses

# Define rho (correlation coefficient)
rho <- 0.6

#### subgroup analyses

## recent violence

# List of dataframes to loop over
dataframes <- list(fsw_data_pv_recent, fsw_data_sv_recent, fsw_data_psv_recent)
dataframe_names <- c("fsw_data_pv_recent", "fsw_data_sv_recent", "fsw_data_psv_recent")

# Define the columns for subgroup analysis
subgroup_columns <- c("ldc_bin", "pre_2016", "recruitment", "perpetrator", "who_region")

# List to store results
results_list <- list()
skipped_results_list <- list()

# Loop through each dataframe
for (i in 1:length(dataframes)) {
  
  # Get the current dataframe and name
  current_df <- dataframes[[i]]
  current_df_name <- dataframe_names[i]
  
  # Filter the dataframe for outcome = "HIV prevalence" and effect_best != "NA"
  filtered_df <- current_df %>% 
    filter(outcome == "HIV prevalence") %>%
    filter(effect_best != "NA")
  
  # Loop through each subgroup column
  for (column in subgroup_columns) {
    # Get unique levels of the current column
    unique_levels <- unique(filtered_df[[column]])
    
    # Loop through each unique level in the current column
    for (level in unique_levels) {
      # Filter data for the current level
      subgroup_df <- filtered_df[filtered_df[[column]] == level, ]
      
      # Handle cases with only one study
      if (nrow(subgroup_df) == 1) {
        message(paste("Only one study for", current_df_name, column, "=", level))
        results_list[[length(results_list) + 1]] <- data.frame(
          dataframe = current_df_name,
          column = column,
          level = level,
          pooled_OR = exp(subgroup_df$effect_best_ln),
          lower_CI = exp(subgroup_df$effect_best_lower_ln),
          upper_CI = exp(subgroup_df$effect_best_upper_ln)
        )
        next
      }
      
      # Create a covariance matrix for the subgroup
      V_mat <- impute_covariance_matrix(subgroup_df$effect_best_var_ln,
                                        cluster = subgroup_df$study_num,
                                        r = rho,
                                        smooth_vi = TRUE)
      
      # Fit the multilevel random effects model for the subgroup
      result <- rma.mv(effect_best_ln, 
                       V = V_mat, 
                       random = ~ 1 | study_num / effect_num,
                       data = subgroup_df,   
                       control=list(rel.tol=1e-8),
                       sparse = TRUE)
      
      # Use `result` to update the `metagen` object
      result2 <- metagen(TE = subgroup_df$effect_best_ln,
                         lower = subgroup_df$effect_best_lower_ln,
                         upper = subgroup_df$effect_best_upper_ln,
                         studlab = subgroup_df$study,
                         data = subgroup_df,
                         sm = "OR",
                         method.tau = "REML",
                         common = FALSE,
                         random = TRUE, 
                         backtransf = TRUE,
                         text.random = "Overall")
      
      # Update result2 with values from the multilevel model
      result2$TE.random <- result$b
      result2$lower.random <- result$ci.lb
      result2$upper.random <- result$ci.ub
      
      # Store the results in the list
      results_list[[length(results_list) + 1]] <- data.frame(
        dataframe = current_df_name,
        column = column,
        level = level,
        pooled_OR = exp(result$b),
        lower_CI = exp(result$ci.lb),
        upper_CI = exp(result$ci.ub)
      )
      
      # Create a unique filename for the forest plot
      filename <- paste0("Plots/subgroups/", current_df_name, "_", column, "_", level, ".png")
      
      # Save the forest plot as a PNG file
      png(filename = filename, width = 38, height = 18, units = "cm", res = 600)
      forest(result2,
             sortvar = subgroup_df$study,
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
  }
}

# Combine the results list and skipped results list
all_results_list <- c(results_list, skipped_results_list)

# Convert the combined results list to a dataframe
results_df <- do.call(rbind, all_results_list)

# Split the results_df into three dataframes
fsw_data_pv_recent_results <- results_df[results_df$dataframe == "fsw_data_pv_recent", ]
fsw_data_sv_recent_results <- results_df[results_df$dataframe == "fsw_data_sv_recent", ]
fsw_data_psv_recent_results <- results_df[results_df$dataframe == "fsw_data_psv_recent", ]

