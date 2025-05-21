# Create power curves for case 0 data
create_power_analysis <- function(results) {
  power_data <- data.frame()
  
  for (key in names(results)) {
    params <- results[[key]]$params
    
    if (params$case_type == 0) {
      power_data <- rbind(power_data, data.frame(
        beta = params$beta,
        sampling_type = params$sampling_label,
        raw_data = results[[key]]$unweighted$reject_rate,
        survey_weighted = results[[key]]$weighted$reject_rate
      ))
    }
  }
  
  # Compute average across all iterations
  power_data <- power_data %>%
    group_by(beta, sampling_type) %>%
    dplyr::summarize(
      raw_data = mean(raw_data, na.rm = TRUE),
      survey_weighted = mean(survey_weighted, na.rm = TRUE),
      .groups = 'drop'
    )
  
  power_table <- power_data %>%
    mutate(sampling_type = ifelse(sampling_type == "constant", 
                                  "Sampling Probability = 0.5", 
                                  "Covariate-Dependent Sampling")) %>%
    select(sampling_type, beta, raw_data, survey_weighted) %>%
    arrange(sampling_type, beta)
  
  power_long <- power_data %>%
    pivot_longer(
      cols = c(raw_data, survey_weighted),
      names_to = "method",
      values_to = "reject_rate"
    ) %>%
    mutate(
      method = ifelse(method == "survey_weighted", "Survey-Weighted", "Raw Data"),
      sampling_label = ifelse(sampling_type == "constant", 
                              "Sampling Probability = 0.5", 
                              "Covariate-Dependent Sampling Probability")
    )
  
  power_plot <- ggplot(power_long, aes(x = beta, y = reject_rate, color = method, group = method)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    facet_wrap(~ sampling_label, ncol = 2) +
    labs(
      title = "Power Curve",
      x = expression(beta[ijt]),
      y = "Rejection Rate",
      color = "Type"
    ) +
    theme_minimal(base_size = 12) +
    scale_y_continuous(limits = c(0, 1)) +
    theme(
      strip.text = element_text(size = 14, face = "bold"),
      strip.background = element_rect(fill = "grey", color = NA),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1, "lines")
    )
  
  return(list(
    plot = power_plot,
    table = power_table
  ))
}

# Create weight plots to analyze the impact of dissimilarity
create_weight_analysis <- function(results) {
  weight_data <- data.frame()
  
  for (key in names(results)) {
    params <- results[[key]]$params
    
    # Only include case 1 and 2 with no treatment effect
    if (params$case_type %in% c(1, 2) && params$beta == 0) {
      # Compute average across all iterations
      mean_weighted <- mean(results[[key]]$weighted$target_ctrl_weights, na.rm = TRUE) 
      mean_unweighted <- mean(results[[key]]$unweighted$target_ctrl_weights, na.rm = TRUE)
      
      case_label <- ifelse(params$case_type == 1, 
                           "Case 1: One Ctrl Similar to Trt", 
                           "Case 2: One Ctrl Diff from Others")
      
      sampling_label <- ifelse(params$sampling_label == "constant", 
                               "Constant Sampling (p=0.5)", 
                               "Covariate-Dependent Sampling")
      
      weight_data <- rbind(weight_data, data.frame(
        case_label = case_label,
        dissimilarity_factor = params$dissimilarity_factor,
        sampling_label = sampling_label,
        raw_data = mean_unweighted,
        survey_weighted = mean_weighted
      ))
    }
  }
  
  weight_data <- weight_data %>%
    group_by(case_label, sampling_label, dissimilarity_factor) %>%
    dplyr::summarize(
      raw_data = mean(raw_data, na.rm = TRUE),
      survey_weighted = mean(survey_weighted, na.rm = TRUE),
      .groups = 'drop'
    )
  
  weight_table <- weight_data %>%
    select(sampling_label, case_label, dissimilarity_factor, raw_data, survey_weighted) %>%
    arrange(sampling_label, case_label, dissimilarity_factor)
  
  weight_long <- weight_data %>%
    pivot_longer(
      cols = c(raw_data, survey_weighted),
      names_to = "method",
      values_to = "weight"
    ) %>%
    mutate(method = ifelse(method == "survey_weighted", "Survey-Weighted", "Raw Data"))
  
  weight_plot <- ggplot(weight_long, aes(x = dissimilarity_factor, y = weight, color = method, group = method)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 4) +
    facet_grid(rows = vars(case_label), cols = vars(sampling_label)) +
    labs(
      title = "Impact of Data Generation on Synthetic Control Weights",
      x = "Dissimilarity Factor",
      y = "Weight of Target Control",
      color = "Method"
    ) +
    theme_minimal(base_size = 12) +
    scale_y_continuous(limits = c(0, 1)) +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "grey", color = NA),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      axis.title = element_text(size = 12, face = "bold"),
      panel.spacing = unit(1, "lines")
    )
  
  return(list(
    plot = weight_plot,
    table = weight_table
  ))
}