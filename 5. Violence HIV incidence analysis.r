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
rightcols <- c("effect", "ci")
rightlabs = c("Estimate", "95% CI")

# constant sampling correlation
rho <- 0.6

# recent exposure and non-missing RR effect
filtered_df <- fsw_data_incidence %>%
  filter(exposure_tf_bin == "Recent", !is.na(effect_best_ln))

# study_num and effect_num columns
filtered_df <- create_study_effect_nums(filtered_df)

# Covariance matrix
V_mat <- impute_covariance_matrix(filtered_df$effect_best_var_ln,
                                  cluster = filtered_df$study_num,
                                  r = rho,
                                  smooth_vi = TRUE)

# multilevel random effects model
result <- rma.mv(filtered_df$effect_best_ln, 
                 V = V_mat, 
                 random = ~ 1 | study_num / effect_num,
                 data = filtered_df,   
                 sparse = TRUE)

print(result)
print(exp(coef(result))) 

# meta-analysis summary
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

# save forest plot
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

# unadjusted analysis
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

# adjusted analysis
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

## Egger's test and funnel plot

# egger's test
eggers <- metabias(result2, method.bias = "linreg")
eggers_p <- if (!is.null(eggers$p.value)) eggers$p.value else NA
eggers_p_str <- if (!is.na(eggers_p)) sprintf("p = %.3f", eggers_p) else ""

# Funnel plot
funnel_label <- paste0(analysis_labels[[analysis]], " - ", exposure_labels[["recent"]])
sanitized_label <- gsub("[/\\?<>\\:*|\"]", "-", funnel_label)
funnel_filename <- paste0("Plots/incidence/funnel plots/", sanitized_label, ".png")

png(filename = funnel_filename, width = 15, height = 15, units = "cm", res = 300)
funnel(result2, main = paste0(funnel_label, "\nEgger's test ", eggers_p_str))
dev.off()

## loop over each analysis for incidence
for (analysis in analyses) {
  V_mat <- impute_covariance_matrix(filtered_df[[paste0(analysis, "_var_ln")]],
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)

  result <- rma.mv(
    yi = filtered_df[[paste0(analysis, "_est_ln")]],
    V = V_mat,
    random = ~ 1 | study_num / effect_num,
    data = filtered_df,
    sparse = TRUE,
    control = list(optimizer = "optim", method = "BFGS", maxiter = 1000, reltol = 1e-8)
  )

  print(result)
  print(exp(coef(result)))
}

## combine funnel plots
get_incidence_funnel_path <- function(analysis, exposure) {
  paste0(
    "Plots/incidence/funnel plots/",
    analysis_labels[[analysis]], " - ", exposure_labels[[exposure]], ".png"
  )
}

all_funnel_imgs <- vector("list", length = length(exposures) * length(analyses))
idx <- 1
for (i in seq_along(exposures)) {
  for (j in seq_along(analyses)) {
    file <- get_incidence_funnel_path(analyses[j], exposures[i])
    if (file.exists(file)) {
      all_funnel_imgs[[idx]] <- rasterGrob(readPNG(file), interpolate = TRUE)
    } else {
      all_funnel_imgs[[idx]] <- nullGrob()
    }
    idx <- idx + 1
  }
}

combined_filename <- "Plots/incidence/funnel plots/incidence_funnel_grid.png"
png(combined_filename, width = 1800, height = 1200, res = 150)
grid.arrange(
  grobs = all_funnel_imgs,
  nrow = length(exposures),
  ncol = length(analyses),
  top = "Incidence: Funnel plots"
)
dev.off()

## subgroup analysis 

# function for "recently exposed to any violence"
process_and_plot(
  data = filtered_df,
  data_name = "filtered_df_recent_incidence",
  output_plot_filename = "Plots/subgroups/recent_any_violence_incidence_subgroup.png"
)

## sensitivity analysis

# Define rho values for sensitivity analysis
rho_values <- c(0.4, 0.6)

# Create folder for sensitivity plots
dir.create("Plots/sensitivity/", recursive = TRUE, showWarnings = FALSE)

# Loop through each rho value
for (rho in rho_values) {
  # Filter for recent exposure and non-missing RR effect
  filtered_df <- fsw_data_incidence %>%
    filter(exposure_tf_bin == "Recent", !is.na(effect_best_ln))
  
  # Create study_num and effect_num columns
  filtered_df <- create_study_effect_nums(filtered_df)
  
  # Covariance matrix
  V_mat <- impute_covariance_matrix(filtered_df$effect_best_var_ln,
                                    cluster = filtered_df$study_num,
                                    r = rho,
                                    smooth_vi = TRUE)
  
  # Multilevel random effects model with adjusted optimizer
  result <- rma.mv(
    yi = filtered_df$effect_best_ln, 
    V = V_mat, 
    random = ~ 1 | study_num / effect_num,
    data = filtered_df,   
    sparse = TRUE,
    control = list(optimizer = "optim", method = "BFGS", maxiter = 1000, reltol = 1e-8)
  )
  
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
  filename <- paste0("Plots/sensitivity/recent_best_incidence_rho_", gsub("\\.", "_", rho), ".png")
  png(filename = filename, width = 50, height = 14, units = "cm", res = 600)
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
}