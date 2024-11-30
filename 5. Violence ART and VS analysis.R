# load packages 
pacman::p_load("meta", "metafor", "readxl", "tidyverse") 

# settings
settings.meta(CIbracket = "(") 
settings.meta(CIseparator = "-") 

# forest plot for unadjusted ART estimates
# ART and VS dataframes
dfs_art_vs <- c("fsw_data_art_pool", "fsw_data_vs_pool")

# unadjusted recent
for (var in dfs_art_vs) {
  df <- get(var)
  meta_analysis_art_vs_unadj_rec(df, var)
}

# adjusted recent
for (var in dfs_art_vs) {
  df <- get(var)
  meta_analysis_art_vs_adj_rec(df, var)
}

# best recent
for (var in dfs_art_vs) {
  df <- get(var)
  meta_analysis_art_vs_best_rec(df, var)
}

# unadjusted ever
for (var in dfs_art_vs) {
  df <- get(var)
  meta_analysis_art_vs_unadj_ever(df, var)
}

# adjusted ever
for (var in dfs_art_vs) {
  df <- get(var)
  meta_analysis_art_vs_adj_ever(df, var)
}

# best ever
for (var in dfs_art_vs) {
  df <- get(var)
  meta_analysis_art_vs_best_ever(df, var)
}











# Unused ------------------------------------------------------------------

filtered_df <- fsw_data_art %>% 
  filter(una_effect != "NR")

art_unadj <- metagen(TE = unadj_est_ln,
                     lower = un_lower_ln,
                     upper = un_upper_ln,
                     studlab = study,
                     data = filtered_df,
                     sm = "OR",
                     method.tau = "DL",
                     comb.fixed = FALSE,
                     comb.random = FALSE, 
                     backtransf = TRUE,
                     byvar = outcome,
                     text.random = "Overall")

# Print summary
print(summary(art_unadj)) 

# Create PNG for forest plot
png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/art_unadj_0511.png", width = 45, height = 28, units = "cm", res = 600)

# Generate forest plot
forest(art_unadj, 
            sortvar = TE,
            xlim = c(0.2, 4),             
            leftcols = c("study", "exposure_definition_short", "exposed_perc_string", "exposure_time_frame"), 
            leftlabs = c("Study", "Exposure definition", "Exposed (%)", "Time frame"),
            rightcols = c("unadj_or_95_2" , "who_region", "perpetrator"), 
            rightlabs = c("OR (95% CI)","WHO region", "Perpetrator"), 
            pooled.totals = FALSE,
            xintercept = 1,
            addrow.overall = TRUE,
            test.subgroup = FALSE,
            overall.hetstat = FALSE,
            overall = FALSE,
            labeltext = TRUE,
            col.subgroup = "black",
            print.subgroup.name = F)

# Close the PNG device
dev.off()

# forest plot for adjusted ART estimates

filtered_df <- fsw_data_art %>% 
  filter(adj_effect != "NR")

art_adj <- metagen(TE = adj_est_ln,
                   lower = adj_lower_ln,
                   upper = adj_upper_ln,
                   studlab = study,
                   data = filtered_df,
                   sm = "OR",
                   method.tau = "DL",
                   comb.fixed = FALSE,
                   comb.random = FALSE, 
                   backtransf = TRUE,
                   byvar = outcome,
                   text.random = "Overall")

# Print summary
print(summary(art_adj)) 

# Create PNG for forest plot
png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/art_adj.png", width = 45, height = 28, units = "cm", res = 600)

# Generate forest plot
forest.meta(art_adj, 
            sortvar = TE,
            xlim = c(0.2, 4),             
            leftcols = c("study", "exposure_definition_short", "exposed_perc_string", "exposure_time_frame"), 
            leftlabs = c("Study", "Exposure definition", "Exposed (%)", "Time frame"),
            rightcols = c("or_95_2" , "who_region", "perpetrator"), 
            rightlabs = c("OR (95% CI)","WHO region", "Perpetrator"), 
            pooled.totals = FALSE,
            xintercept = 1,
            addrow.overall = TRUE,
            test.subgroup = FALSE,
            overall.hetstat = FALSE,
            overall = FALSE,
            labeltext = TRUE,
            col.subgroup = "black",
            print.subgroup.name = F)

# Close the PNG device
dev.off()


# forest plot for combined ART estimates

art_best <- metagen(TE = effect_best_ln,
                    lower = effect_best_lower_ln,
                    upper = effect_best_upper_ln,
                    studlab = study,
                    data = fsw_data_art,
                    sm = "OR",
                    method.tau = "DL",
                    comb.fixed = FALSE,
                    comb.random = FALSE, 
                    backtransf = TRUE,
                    byvar = outcome,
                    text.random = "Overall")

# Print summary
print(summary(art_best)) 

# Create PNG for forest plot
png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/art_best.png", width = 45, height = 28, units = "cm", res = 600)

# Generate forest plot
forest.meta(art_best, 
            sortvar = TE,
            xlim = c(0.2, 4),             
            leftcols = c("study", "exposure_definition_short", "exposed_perc_string", "exposure_time_frame"), 
            leftlabs = c("Study", "Exposure definition", "Exposed (%)", "Time frame"),
            rightcols = c("or_95_2" , "who_region", "perpetrator"), 
            rightlabs = c("OR (95% CI)","WHO region", "Perpetrator"), 
            pooled.totals = FALSE,
            xintercept = 1,
            addrow.overall = TRUE,
            test.subgroup = FALSE,
            overall.hetstat = FALSE,
            overall = FALSE,
            labeltext = TRUE,
            col.subgroup = "black",
            print.subgroup.name = F)

# Close the PNG device
dev.off()

# forest plot for unadjusted ART estimates

filtered_df <- fsw_data_art %>% 
  filter(una_effect != "NR")

art_unadj <- metagen(TE = unadj_est_ln,
                     lower = un_lower_ln,
                     upper = un_upper_ln,
                     studlab = study,
                     data = filtered_df,
                     sm = "OR",
                     method.tau = "DL",
                     comb.fixed = FALSE,
                     comb.random = FALSE, 
                     backtransf = TRUE,
                     byvar = outcome,
                     text.random = "Overall")

# Print summary
print(summary(art_unadj)) 

# Create PNG for forest plot
png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/art_unadj.png", width = 45, height = 28, units = "cm", res = 600)

# Generate forest plot
forest.meta(art_unadj, 
            sortvar = TE,
            xlim = c(0.2, 4),             
            leftcols = c("study", "exposure_definition_short", "exposed_perc_string", "exposure_time_frame"), 
            leftlabs = c("Study", "Exposure definition", "Exposed (%)", "Time frame"),
            rightcols = c("or_95_2" , "who_region", "perpetrator"), 
            rightlabs = c("OR (95% CI)","WHO region", "Perpetrator"), 
            pooled.totals = FALSE,
            xintercept = 1,
            addrow.overall = TRUE,
            test.subgroup = FALSE,
            overall.hetstat = FALSE,
            overall = FALSE,
            labeltext = TRUE,
            col.subgroup = "black",
            print.subgroup.name = F)

# Close the PNG device
dev.off()

# forest plot for adjusted ART estimates

filtered_df <- fsw_data_art %>% 
  filter(adj_effect != "NR")

art_adj <- metagen(TE = adj_est_ln,
                   lower = adj_lower_ln,
                   upper = adj_upper_ln,
                   studlab = study,
                   data = filtered_df,
                   sm = "OR",
                   method.tau = "DL",
                   comb.fixed = FALSE,
                   comb.random = FALSE, 
                   backtransf = TRUE,
                   byvar = outcome,
                   text.random = "Overall")

# Print summary
print(summary(art_adj)) 

# Create PNG for forest plot
png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/art_adj.png", width = 45, height = 28, units = "cm", res = 600)

# Generate forest plot
forest.meta(art_adj, 
            sortvar = TE,
            xlim = c(0.2, 4),             
            leftcols = c("study", "exposure_definition_short", "exposed_perc_string", "exposure_time_frame"), 
            leftlabs = c("Study", "Exposure definition", "Exposed (%)", "Time frame"),
            rightcols = c("or_95_2" , "who_region", "perpetrator"), 
            rightlabs = c("OR (95% CI)","WHO region", "Perpetrator"), 
            pooled.totals = FALSE,
            xintercept = 1,
            addrow.overall = TRUE,
            test.subgroup = FALSE,
            overall.hetstat = FALSE,
            overall = FALSE,
            labeltext = TRUE,
            col.subgroup = "black",
            print.subgroup.name = F)

# Close the PNG device
dev.off()


# forest plot for combined ART estimates

art_best <- metagen(TE = effect_best_ln,
                    lower = effect_best_lower_ln,
                    upper = effect_best_upper_ln,
                    studlab = study,
                    data = fsw_data_art,
                    sm = "OR",
                    method.tau = "DL",
                    comb.fixed = FALSE,
                    comb.random = FALSE, 
                    backtransf = TRUE,
                    byvar = outcome,
                    text.random = "Overall")

# Print summary
print(summary(art_best)) 

# Create PNG for forest plot
png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/art_adj.png", width = 45, height = 28, units = "cm", res = 600)

# Generate forest plot
forest.meta(art_adj, 
            sortvar = TE,
            xlim = c(0.2, 4),             
            leftcols = c("study", "exposure_definition_short", "exposed_perc_string", "exposure_time_frame"), 
            leftlabs = c("Study", "Exposure definition", "Exposed (%)", "Time frame"),
            rightcols = c("or_95_2" , "who_region", "perpetrator"), 
            rightlabs = c("OR (95% CI)","WHO region", "Perpetrator"), 
            pooled.totals = FALSE,
            xintercept = 1,
            addrow.overall = TRUE,
            test.subgroup = FALSE,
            overall.hetstat = FALSE,
            overall = FALSE,
            labeltext = TRUE,
            col.subgroup = "black",
            print.subgroup.name = F)

# Close the PNG device
dev.off()

# viral suppression

# forest plot for unadjusted viral load estimates

filtered_df <- fsw_data_vs %>% 
  filter(una_effect != "NR")

vs_unadj <- metagen(TE = unadj_est_ln,
                    lower = un_lower_ln,
                    upper = un_upper_ln,
                    studlab = study,
                    data = filtered_df,
                    sm = "OR",
                    method.tau = "DL",
                    comb.fixed = FALSE,
                    comb.random = FALSE, 
                    backtransf = TRUE,
                    byvar = outcome,
                    text.random = "Overall")

# Print summary
print(summary(vs_unadj)) 

# Create PNG for forest plot
png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/vs_unadj.png", width = 45, height = 28, units = "cm", res = 600)

# Generate forest plot
forest.meta(vs_unadj, 
            sortvar = TE,
            xlim = c(0.2, 4),             
            leftcols = c("study", "exposure_definition_short", "exposed_perc_string", "exposure_time_frame"), 
            leftlabs = c("Study", "Exposure definition", "Exposed (%)", "Time frame"),
            rightcols = c("or_95_2" , "outcome_definition", "who_region", "perpetrator"), 
            rightlabs = c("OR (95% CI)","VS definition",  "WHO region", "Perpetrator"), 
            pooled.totals = FALSE,
            xintercept = 1,
            addrow.overall = TRUE,
            test.subgroup = FALSE,
            overall.hetstat = FALSE,
            overall = FALSE,
            labeltext = TRUE,
            col.subgroup = "black",
            print.subgroup.name = F)

# Close the PNG device
dev.off()

# forest plot for adjusted ART estimates

filtered_df <- fsw_data_vs %>% 
  filter(adj_effect != "NR")

vs_adj <- metagen(TE = adj_est_ln,
                  lower = adj_lower_ln,
                  upper = adj_upper_ln,
                  studlab = study,
                  data = filtered_df,
                  sm = "OR",
                  method.tau = "DL",
                  comb.fixed = FALSE,
                  comb.random = FALSE, 
                  backtransf = TRUE,
                  byvar = outcome,
                  text.random = "Overall")

# Print summary
print(summary(vs_adj)) 

# Create PNG for forest plot
png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/vs_adj.png", width = 45, height = 28, units = "cm", res = 600)

# Generate forest plot
forest(vs_adj, 
            sortvar = TE,
            xlim = c(0.2, 4),             
            leftcols = c("study", "exposure_definition_short", "exposed_perc_string", "exposure_time_frame"), 
            leftlabs = c("Study", "Exposure definition", "Exposed (%)", "Time frame"),
            rightcols = c("or_95_2" , "outcome_definition", "who_region", "perpetrator"), 
            rightlabs = c("OR (95% CI)","VS definition",  "WHO region", "Perpetrator"), 
            pooled.totals = FALSE,
            xintercept = 1,
            addrow.overall = TRUE,
            test.subgroup = FALSE,
            overall.hetstat = FALSE,
            overall = FALSE,
            labeltext = TRUE,
            col.subgroup = "black",
            print.subgroup.name = F)

# Close the PNG device
dev.off()


# forest plot for combined vs estimates

vs_best <- metagen(TE = effect_best_ln,
                   lower = effect_best_lower_ln,
                   upper = effect_best_upper_ln,
                   studlab = study,
                   data = fsw_data_vs,
                   sm = "OR",
                   method.tau = "DL",
                   comb.fixed = FALSE,
                   comb.random = FALSE, 
                   backtransf = TRUE,
                   byvar = outcome,
                   text.random = "Overall")

# Print summary
print(summary(vs_best)) 

# Create PNG for forest plot
png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/vs_best.png", width = 45, height = 28, units = "cm", res = 600)

# Generate forest plot
forest.meta(vs_best, 
            sortvar = TE,
            xlim = c(0.2, 4),             
            leftcols = c("study", "exposure_definition_short", "exposed_perc_string", "exposure_time_frame"), 
            leftlabs = c("Study", "Exposure definition", "Exposed (%)", "Time frame"),
            rightcols = c("or_95_2" , "outcome_definition", "who_region", "perpetrator"), 
            rightlabs = c("OR (95% CI)","VS definition",  "WHO region", "Perpetrator"), 
            pooled.totals = FALSE,
            xintercept = 1,
            addrow.overall = TRUE,
            test.subgroup = FALSE,
            overall.hetstat = FALSE,
            overall = FALSE,
            labeltext = TRUE,
            col.subgroup = "black",
            print.subgroup.name = F)

# Close the PNG device
dev.off()

# Print summary
print(summary(art_best)) 

# Create PNG for forest plot
png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/art_best.png", width = 45, height = 28, units = "cm", res = 600)

# Generate forest plot
forest.meta(art_best, 
            sortvar = TE,
            xlim = c(0.2, 4),             
            leftcols = c("study", "exposure_definition_short", "exposed_perc_string", "exposure_time_frame"), 
            leftlabs = c("Study", "Exposure definition", "Exposed (%)", "Time frame"),
            rightcols = c("or_95_2" , "who_region", "perpetrator"), 
            rightlabs = c("OR (95% CI)","WHO region", "Perpetrator"), 
            pooled.totals = FALSE,
            xintercept = 1,
            addrow.overall = TRUE,
            test.subgroup = FALSE,
            overall.hetstat = FALSE,
            overall = FALSE,
            labeltext = TRUE,
            col.subgroup = "black",
            print.subgroup.name = F)

# Close the PNG device
dev.off()







