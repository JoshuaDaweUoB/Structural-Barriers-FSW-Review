# load packages 
pacman::p_load("meta", "metafor", "readxl", "openxlsx", "tidyverse", "kableExtra", "robumeta", "clubSandwich", "grid", "png", "gridExtra") 

# set working directory
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence")

# settings
settings.meta(CIbracket = "(") 
settings.meta(CIseparator = "-") 

# columns
leftcols_recent <- c("study", "study_num", "effect_num", "exposure_definition_short", "exposure_time_frame", "perpetrator", "country")
leftlabs_recent <- c("Study", "Study number", "Effect number", "Exposure definition", "Exposure time frame", "Perpetrator", "Country")
leftcols_lifetime <- c("study", "study_num", "effect_num", "exposure_definition_short", "perpetrator", "country")
leftlabs_lifetime <- c("Study", "Study number", "Effect number", "Exposure definition", "Perpetrator", "Country")
rightcols <- c("effect", "ci")
rightlabs = c("Estimate", "95% CI")

# constant sampling correlation
rho <- 0.6
View(fsw_data_incidence)
# Filter for recent exposure and non-missing RR effect
filtered_df <- fsw_data_incidence %>%
  filter(exposure_tf_bin == "Recent", !is.na(effect_best_ln))

# Add study_num and effect_num columns
filtered_df <- create_study_effect_nums(filtered_df)

# Covariance matrix
V_mat <- impute_covariance_matrix(filtered_df$effect_best_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)

# Multilevel random effects model
result <- rma.mv(filtered_df$effect_best_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)

print(result)
print(exp(coef(result))) 

# Meta-analysis summary
result2 <- metagen(
  TE = filtered_df$effect_best_ln,
  lower = filtered_df$effect_best_lower_ln,
  upper = filtered_df$effect_best_upper_ln,
  studlab = filtered_df$study,
  data = filtered_df,
  sm = "RR",
  method.tau = "REML",
  common = FALSE,
  random = TRUE, 
  backtransf = TRUE,
  text.random = "Overall"
)

print(summary(result2))

result2$TE.random <- result$b
result2$lower.random <- result$ci.lb
result2$upper.random <- result$ci.ub

# Save forest plot
filename <- "Plots/incidence/recent_best_incidence_rr.png"
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
       col.subgroup = "black"
)
dev.off()

# Unadjusted analysis
filtered_df_unadj <- fsw_data_incidence %>%
  filter(exposure_tf_bin == "Recent", !is.na(unadj_est_ln))

filtered_df_unadj <- create_study_effect_nums(filtered_df_unadj)

V_mat_unadj <- impute_covariance_matrix(filtered_df_unadj$unadj_var_ln,
                                        cluster = filtered_df_unadj$study_num,
                                        r = rho,
                                        smooth_vi = TRUE)

result_unadj <- rma.mv(
  filtered_df_unadj$unadj_est_ln, 
  V = V_mat_unadj, 
  random = ~ 1 | study_num / effect_num,
  data = filtered_df_unadj,   
  sparse = TRUE,
  control = list(optimizer = "optim", method = "BFGS")
)

print(result_unadj)
print(exp(coef(result_unadj))) 

result2_unadj <- metagen(
  TE = filtered_df_unadj$unadj_est_ln,
  lower = filtered_df_unadj$un_lower_ln,
  upper = filtered_df_unadj$un_upper_ln,
  studlab = filtered_df_unadj$study,
  data = filtered_df_unadj,
  sm = "RR",
  method.tau = "REML",
  common = FALSE,
  random = TRUE, 
  backtransf = TRUE,
  text.random = "Overall"
)

print(summary(result2_unadj))

result2_unadj$TE.random <- result_unadj$b
result2_unadj$lower.random <- result_unadj$ci.lb
result2_unadj$upper.random <- result_unadj$ci.ub

filename_unadj <- "Plots/incidence/recent_unadj_incidence_rr.png"
png(filename = filename_unadj, width = 45, height = 14, units = "cm", res = 600)
forest(result2_unadj,
       sortvar = filtered_df_unadj$study,
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
       col.subgroup = "black"
)
dev.off()

# Adjusted analysis
filtered_df_adj <- fsw_data_incidence %>%
  filter(exposure_tf_bin == "Recent", !is.na(adj_est_ln))

filtered_df_adj <- create_study_effect_nums(filtered_df_adj)

V_mat_adj <- impute_covariance_matrix(filtered_df_adj$adj_var_ln,
                                      cluster = filtered_df_adj$study_num,
                                      r = rho,
                                      smooth_vi = TRUE)

result_adj <- rma.mv(filtered_df_adj$adj_est_ln, 
                     V = V_mat_adj, 
                     random = ~ 1 | study_num / effect_num,
                     data = filtered_df_adj,   
                     sparse = TRUE)

print(result_adj)
print(exp(coef(result_adj))) 

result2_adj <- metagen(
  TE = filtered_df_adj$adj_est_ln,
  lower = filtered_df_adj$adj_lower_ln,
  upper = filtered_df_adj$adj_upper_ln,
  studlab = filtered_df_adj$study,
  data = filtered_df_adj,
  sm = "RR",
  method.tau = "REML",
  common = FALSE,
  random = TRUE, 
  backtransf = TRUE,
  text.random = "Overall"
)

print(summary(result2_adj))

result2_adj$TE.random <- result_adj$b
result2_adj$lower.random <- result_adj$ci.lb
result2_adj$upper.random <- result_adj$ci.ub

filename_adj <- "Plots/incidence/recent_adj_incidence_rr.png"
png(filename = filename_adj, width = 45, height = 14, units = "cm", res = 600)
forest(result2_adj,
       sortvar = filtered_df_adj$study,
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
       col.subgroup = "black"
)
dev.off()