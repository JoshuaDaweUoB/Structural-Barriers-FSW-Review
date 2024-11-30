# load packages 
pacman::p_load("meta", "metafor", "readxl", "tidyverse") 

# settings
settings.meta(CIbracket = "(") 
settings.meta(CIseparator = "-") 

# ever summary plot
summary_hiv_violence_ever  <- metagen(TE = effect_ln,
                                      lower = lower_ln,
                                      upper = upper_ln,
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

png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/summary_hiv_violence_ever.png", width = 24, height = 12, units = "cm",
    res = 600)

summary_hiv_violence_ever <- forest.meta(summary_hiv_violence_ever, 
                                         sortvar = name,
                                         xlim=c(0.2, 4),             
                                         leftcols = c("name", "studies"), 
                                         leftlabs = c("Pooled exposure", "Studies"),
                                         rightcols = c("or_95_2", "fuck"), 
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

## recent summary plot
summary_hiv_violence_rec  <- metagen(TE = effect_ln,
                                     lower = lower_ln,
                                     upper = upper_ln,
                                     studlab = name,
                                     data = summary_violence_rec,
                                     sm = "OR",
                                     method.tau = "DL",
                                     comb.fixed = FALSE,
                                     comb.random = FALSE, 
                                     backtransf = TRUE,
                                     byvar = Adjust,
                                     text.random = "Overall")

summary(summary_hiv_violence_rec) 

png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/Plots/summary_hiv_violence_recent.png", width = 24, height = 12, units = "cm",
    res = 600)

summary_hiv_violence_rec <- forest.meta(summary_hiv_violence_rec, 
                                        sortvar = name,
                                        xlim=c(0.2, 4),             
                                        leftcols = c("name", "studies"), 
                                        leftlabs = c("Pooled exposure", "Studies"),
                                        rightcols = c("or_95_2", "fuck"), 
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

# all physical violence and hiv infection data


# all ART data
summary_hiv_art  <- metagen(TE = effect_best_ln,
                            lower = effect_best_lower_ln,
                            upper = effect_best_upper_ln,
                            studlab = study,
                            data = fsw_data_art,
                            sm = "OR",
                            method.tau = "DL",
                            comb.fixed = FALSE,
                            comb.random = FALSE, 
                            backtransf = TRUE,
                            byvar = exposure_type,
                            text.random = "Overall")

summary(summary_hiv_art) 

png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/hiv_art.png", width = 20, height = 18, units = "cm",
    res = 500)

summary_hiv_art <- forest.meta(summary_hiv_art, 
                               sortvar = study,
                               xlim=c(0.2, 4),             
                               leftcols = c("study", "outcome"), 
                               leftlabs = c("Study", "Outcome"),
                               rightcols = c("exposure_time_frame", "perpetrator"), 
                               rightlabs = c("Time frame", "Perpetrator"), 
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

# all VS data
summary_hiv_vs  <- metagen(TE = effect_best_ln,
                           lower = effect_best_lower_ln,
                           upper = effect_best_upper_ln,
                           studlab = study,
                           data = fsw_data_vs,
                           sm = "OR",
                           method.tau = "DL",
                           comb.fixed = FALSE,
                           comb.random = FALSE, 
                           backtransf = TRUE,
                           byvar = exposure_type,
                           text.random = "Overall")

summary(summary_hiv_vs) 

png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/hiv_vs.png", width = 22, height = 10, units = "cm",
    res = 500)

summary_hiv_vs <- forest.meta(summary_hiv_vs, 
                              sortvar = study,
                              xlim=c(0.2, 4),             
                              leftcols = c("study", "outcome_definition"), 
                              leftlabs = c("Study", "VS definition"),
                              rightcols = c("exposure_time_frame", "perpetrator"), 
                              rightlabs = c("Time frame", "Perpetrator"), 
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

# SV and HIV data
summary_hiv_sv  <- metagen(TE = effect_best_ln,
                           lower = effect_best_lower_ln,
                           upper = effect_best_upper_ln,
                           studlab = study,
                           data = fsw_data_sv_hiv,
                           sm = "OR",
                           method.tau = "DL",
                           comb.fixed = FALSE,
                           comb.random = FALSE, 
                           backtransf = TRUE,
                           text.random = "Overall")

summary(summary_hiv_sv) 

png(filename = "C:/Users/vl22683/OneDrive - University of Bristol/Documents/Misc/UNAIDS/FSW/Analysis/Violence/hiv_sv.png", width = 25, height = 25, units = "cm",
    res = 500)

summary_hiv_sv <- forest.meta(summary_hiv_sv, 
                              sortvar = study,
                              xlim=c(0.2, 4),             
                              leftcols = c("study"), 
                              leftlabs = c("Study"),
                              rightcols = c("exposure_time_frame", "perpetrator"), 
                              rightlabs = c("Time frame", "Perpetrator"), 
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
