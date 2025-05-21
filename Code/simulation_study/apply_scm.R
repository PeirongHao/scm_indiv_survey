# Compute control region averages with/without survey weights
ctrl_avg_simul <- function(df, use_weights = FALSE) {
  df %>% 
    filter(region != 1 & S == 1) %>%
    group_by(region, year) %>%
    dplyr::summarize(
      Z1 = mean(Z1, na.rm = TRUE), 
      Z2 = mean(Z2, na.rm = TRUE),
      X1 = ifelse(use_weights, weighted.mean(X1, W, na.rm = TRUE), mean(X1, na.rm = TRUE)),
      X2 = ifelse(use_weights, weighted.mean(X2, W, na.rm = TRUE), mean(X2, na.rm = TRUE)),
      Y = ifelse(use_weights, weighted.mean(Y, W, na.rm = TRUE), mean(Y, na.rm = TRUE)),
      .groups = 'drop'
    )
}

# Add treated region average to control average
trt_avg_simul <- function(df, df.ctrl) {
  df.trt <- df %>% 
    filter(region == 1 & S == 1) %>%
    group_by(year) %>%
    dplyr::summarize(
      region = 1,
      Z1 = mean(Z1, na.rm = TRUE),
      Z2 = mean(Z2, na.rm = TRUE),
      X1 = weighted.mean(X1, W, na.rm = TRUE),
      X2 = weighted.mean(X2, W, na.rm = TRUE),
      Y = weighted.mean(Y, W, na.rm = TRUE),
      .groups = 'drop'
    )
  
  return(bind_rows(df.ctrl, df.trt))
}

# Run synthetic control and extract results
run_synthetic_control_simul <- function(df.simul, T0 = 40, use_weights = FALSE, target_district = 2) {
  df.avg <- trt_avg_simul(df.simul, ctrl_avg_simul(df.simul, use_weights))
  
  df_out <- df.avg %>%
    synthetic_control(
      outcome = Y, 
      unit = region, 
      time = year, 
      i_unit = "1", 
      i_time = T0, 
      generate_placebos = TRUE
    ) %>%
    generate_predictor(
      time_window = 1:T0,
      Z1 = mean(Z1, na.rm = TRUE),
      Z2 = mean(Z2, na.rm = TRUE),
      X1 = mean(X1, na.rm = TRUE),
      X2 = mean(X2, na.rm = TRUE)
    )
  
  # Add outcome values for each pre-treatment period
  for (t in 1:T0) {
    df_out <- df_out %>% 
      generate_predictor(time_window = t, !!paste0("Y_", t) := Y)
  }
  
  df_out <- df_out %>%
    generate_weights(
      optimization_window = 1:T0,
      optimization_method = "L-BFGS-B",
      margin_ipop = .01, sigf_ipop = 7, bound_ipop = 6
    ) %>%
    generate_control()
  
  # Extract p-value
  p_value_tbl <- df_out %>% grab_significance()
  p_value_trt <- p_value_tbl %>% filter(type == 'Treated') %>% pull(fishers_exact_pvalue)
  
  # Extract treatment effect
  effect <- df_out %>% grab_synthetic_control(placebo = TRUE) %>%
    filter(.placebo == 0, time_unit == (T0+1)) %>%
    dplyr::summarize(effect = real_y - synth_y) %>%
    pull(effect)
  
  # Extract weight for target control
  if (!is.null(target_district)) {
    target_ctrl_weight <- df_out %>% grab_unit_weights() %>% 
      filter(unit == target_district) %>% 
      pull(weight)
  } else {
    target_ctrl_weight <- NULL
  }
  
  return(list(
    p_value = p_value_trt, 
    effect = effect,
    target_ctrl_weight = target_ctrl_weight
  ))
}