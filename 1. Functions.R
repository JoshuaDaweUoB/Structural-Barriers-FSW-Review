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

# create study number and effect number
create_study_effect_nums <- function(df) {
  # Sort the data by author and title
  df <- df %>% arrange(author, title)
  
  # Create study_num and effect_num columns
  df <- df %>%
    mutate(effect_num = row_number()) %>%
    mutate(study_num = cumsum(!duplicated(title)))
  
  return(df)
}

## subgroup function

# lists for loops and functions
subgroup_columns <- c("ldc_bin", "lmic_bin", "who_region", "hiv_decrim", "sw_decrim", "gbv_law", "pre_2016", "recruitment", "perpetrator", "rob_score_3cat")

# Define a reusable function for processing, modeling, and plotting subgroups
process_and_plot <- function(data, data_name, output_plot_filename) {
  # Create a list to store the dataframes for all subgroup columns
  all_subgroup_dataframes <- list()
  
  # Loop through each column in subgroup_columns
  for (column in subgroup_columns) {
    unique_levels <- unique(data[[column]])
    message(paste("Processing column:", column))
    print(unique_levels)
    
    # Loop through each level of the current column
    for (level in unique_levels) {
      filtered_df <- data %>% filter(.data[[column]] == level)
      dataframe_name <- paste0(data_name, "_", column, "_", level)
      all_subgroup_dataframes[[dataframe_name]] <- filtered_df
    }
  }
  
  # Create an empty dataframe to store subgroup results
  subgroup_results <- data.frame(
    subgroup_level = character(),
    subgroup = character(),
    model_coef = numeric(),
    model_ci_lower = numeric(),
    model_ci_upper = numeric(),
    studies = numeric(),
    estimates = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each subgroup dataframe to fit models
  for (name in names(all_subgroup_dataframes)) {
    filtered_df <- all_subgroup_dataframes[[name]]
    if (nrow(filtered_df) == 0) next
    
    if (!"study_num" %in% colnames(filtered_df) || !"effect_num" %in% colnames(filtered_df)) {
      filtered_df <- create_study_effect_nums(filtered_df)
      filtered_df <- filtered_df %>% filter(use == "yes")
    }
    
    if (!"effect_best_var_ln" %in% colnames(filtered_df) || !"effect_best_ln" %in% colnames(filtered_df)) next
    
    rho <- 0.6
   
    tryCatch({
      if (nrow(filtered_df) > 1) {
        V_mat <- impute_covariance_matrix(filtered_df$effect_best_var_ln,
                                          cluster = filtered_df$study_num,
                                          r = rho,
                                          smooth_vi = TRUE)

       result <- rma.mv(
  yi = filtered_df$effect_best_ln,  # Specify the column for effect sizes
  V = V_mat, 
  random = ~ 1 | study_num / effect_num,
  data = filtered_df,
  control = list(optimizer = "optim", maxiter = 100000, method = "BFGS", stepadj = 0.1),
  sparse = TRUE
)
                         
        coef <- exp(coef(result))
        ci_lower <- exp(result$ci.lb)
        ci_upper <- exp(result$ci.ub)
      } else {
        coef <- exp(filtered_df$effect_best_ln)
        ci_lower <- exp(filtered_df$effect_best_ln - 1.96 * sqrt(filtered_df$effect_best_var_ln))
        ci_upper <- exp(filtered_df$effect_best_ln + 1.96 * sqrt(filtered_df$effect_best_var_ln))
      }
      
      subgroup_name <- sub(paste0(data_name, "_"), "", name)
      subgroup <- sub("^(.*?)_.*$", "\\1", subgroup_name)
      subgroup_level <- sub("^[^_]*_", "", subgroup_name)
      
      num_studies <- length(unique(filtered_df$study_num))
      num_estimates <- nrow(filtered_df)
      
      subgroup_results <- rbind(
        subgroup_results,
        data.frame(
          subgroup_level = subgroup_level,
          subgroup = subgroup,
          model_coef = coef,
          model_ci_lower = ci_lower,
          model_ci_upper = ci_upper,
          studies = num_studies,
          estimates = num_estimates,
          stringsAsFactors = FALSE
        )
      )
    }, error = function(e) {
      message(paste("Error processing dataframe:", name))
      message(e)
    })
  }
  
  # Drop rows with missing values
  subgroup_results <- subgroup_results %>% drop_na()
  
  # Log-transform the results
  subgroup_results <- subgroup_results %>%
    mutate(
      log_model_coef = log(model_coef),
      log_model_ci_lower = log(model_ci_lower),
      log_model_ci_upper = log(model_ci_upper)
    )
  
  # Rename rows
  subgroup_results <- subgroup_results %>%
    mutate(subgroup_level = case_when(
    subgroup_level == "bin_no" ~ "No",
    subgroup_level == "bin_yes" ~ "Yes",
    subgroup_level == "decrim_partial" ~ "Partial",
    subgroup_level == "decrim_yes" ~ "Yes",
    subgroup_level == "decrim_no" ~ "No",
    subgroup_level == "law_yes" ~ "Yes",
    subgroup_level == "law_no" ~ "No",
    subgroup_level == "law_NA" ~ "NA",
    subgroup_level == "2016_FALSE" ~ "No",
    subgroup_level == "2016_TRUE" ~ "Yes",
    subgroup_level == "decrim_NA" ~ "NA",
    subgroup_level == "region_African Region" ~ "African Region",
    subgroup_level == "region_Region of the Americas" ~ "Region of the Americas",
    subgroup_level == "region_South-East Asia Region" ~ "South-East Asia Region",
    subgroup_level == "region_European Region" ~ "European Region",
    subgroup_level == "region_Eastern Mediterranean Region" ~ "Eastern Mediterranean Region",
    subgroup_level == "region_Western Pacific Region" ~ "Western Pacific Region",
    subgroup_level == "score_3cat_Good/ very good" ~ "Good/very good",
    subgroup_level == "score_3cat_Satisfactory" ~ "Satisfactory",
    subgroup_level == "score_3cat_Unsatisfactory" ~ "Unsatisfactory",
    TRUE ~ subgroup_level
    )) %>%
    mutate(subgroup = case_when(
      subgroup == "ldc" ~ "Least developed country",
      subgroup == "lmic" ~ "Lower-middle income country",
      subgroup == "hiv" ~ "HIV decriminalisation",
      subgroup == "gbv" ~ "Gender-based violence laws",
      subgroup == "sw" ~ "Decriminalised sex work",
      subgroup == "pre" ~ "Published before 2016",
      subgroup == "recruitment" ~ "Recruitment",
      subgroup == "perpetrator" ~ "Perpetrator",
      subgroup == "who" ~ "WHO region",
      subgroup == "rob" ~ "Risk of bias score",
      TRUE ~ subgroup
    ))
  
  # Perform meta-analysis
  forest_plot <- metagen(
    TE = log_model_coef,
    lower = log_model_ci_lower,
    upper = log_model_ci_upper,
    data = subgroup_results,
    sm = "OR",
    method.tau = "DL",
    comb.fixed = FALSE,
    comb.random = FALSE, 
    backtransf = TRUE,
    byvar = subgroup, 
    text.random = "Overall"
  )
  
  # Save forest plot
  png(filename = output_plot_filename, width = 30, height = 30, units = "cm", res = 600)
  forest(
    forest_plot,
    sortvar = subgroup,
    xlim = c(0.2, 4),             
    leftcols = c("subgroup_level", "studies", "estimates"), 
    leftlabs = c("Category", "Nb studies", "Nb estimates"),
    pooled.totals = FALSE,
    xintercept = 1,
    addrow.overall = TRUE,
    test.subgroup = FALSE,
    overall.hetstat = FALSE,
    overall = FALSE,
    labeltext = TRUE,
    col.subgroup = "black",
    print.subgroup.name = FALSE
  )
  dev.off()
  
  # Print summary
  print(summary(forest_plot))
}
