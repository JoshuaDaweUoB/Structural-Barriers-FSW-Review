#  load packages
pacman::p_load("tidyverse")

# Cleaning functions ------------------------------------------------------

# convert to numeric
convert_to_numeric <- function(df) {
  df <- transform(df, 
                  unadj_est = as.numeric(unadj_est), 
                  un_lower = as.numeric(un_lower),
                  un_upper = as.numeric(un_upper),
                  adj_est = as.numeric(adj_est), 
                  adj_lower = as.numeric(adj_lower),
                  adj_upper = as.numeric(adj_upper),
                  effect_best = as.numeric(effect_best), 
                  effect_best_lower = as.numeric(effect_best_lower),
                  effect_best_upper = as.numeric(effect_best_upper))
  return(df)
}

# log transform
log_transform <- function(df) {
  df <- transform(df, 
                  unadj_est_ln = log(unadj_est),
                  un_lower_ln = log(un_lower),
                  un_upper_ln = log(un_upper),
                  adj_est_ln = log(adj_est),
                  adj_lower_ln = log(adj_lower),
                  adj_upper_ln = log(adj_upper),
                  effect_best_ln = log(effect_best),
                  effect_best_lower_ln = log(effect_best_lower),
                  effect_best_upper_ln = log(effect_best_upper))
  return(df)
}

# log transform variance ORs
log_transform_var <- function(df) {
  df <- transform(df, 
                  unadj_var_ln = ((un_upper_ln - un_lower_ln) / (2 * 1.96))^2,
                  adj_var_ln = ((adj_upper_ln - adj_lower_ln) / (2 * 1.96))^2,
                  effect_best_var_ln = ((effect_best_upper_ln - effect_best_lower_ln) / (2 * 1.96))^2)
  return(df)
}

# concatenate unadjusted estimates
format_violence_data_unadj <- function(df) {
  df$effect_unadj_str <- sprintf("%.2f", df$unadj_est)
  df$unaj_lower_str <- sprintf("%.2f", df$un_lower)
  df$unaj_upper_str <- sprintf("%.2f", df$un_upper)
  df$unadj_or_95 <- paste0(df$unaj_lower_str, "–", df$unaj_upper_str)
  df$unadj_or_95_2 <- paste0(df$effect_unadj_str, " (", df$unadj_or_95, ")")
  return(df)
}

# concatenate adjusted estimates
format_violence_data_adj <- function(df) {
  df$effect_adj_str <- sprintf("%.2f", df$adj_est)
  df$adj_lower_str <- sprintf("%.2f", df$adj_lower)
  df$adj_upper_str <- sprintf("%.2f", df$adj_lower)
  df$adj_or_95 <- paste0(df$adj_lower_str, "–", df$adj_upper_str)
  df$adj_or_95_2 <- paste0(df$effect_adj__str, " (", df$adj_or_95, ")")
  return(df)
}

# concatenate best estimates
format_violence_data <- function(df) {
  df$effect_best_str <- sprintf("%.2f", df$effect_best)
  df$lower_str <- sprintf("%.2f", df$effect_best_lower)
  df$upper_str <- sprintf("%.2f", df$effect_best_upper)
  df$or_95 <- paste0(df$lower_str, "–", df$upper_str)
  df$or_95_2 <- paste0(df$effect_best_str, " (", df$or_95, ")")
  return(df)
}
