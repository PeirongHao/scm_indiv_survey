# Run across two SCM approaches across multiple iterations
run_simulation <- function(max_iter = 2, T0 = 40, num_regions = 21, beta_ijt = 0,
                           case_type = 0, dissimilarity_factor = 0,
                           target_district = 2, save_data = TRUE,
                           use_covariates_for_sampling = FALSE) {
  weighted_results <- list()
  unweighted_results <- list()
  
  for (iter in 1:max_iter) {
    # Generate data
    df.simul <- generate_data_simul(
      num_regions = num_regions,
      T0 = T0,
      beta_ijt = beta_ijt,
      case_type = case_type,
      dissimilarity_factor = dissimilarity_factor,
      target_district = target_district,
      use_covariates_for_sampling = use_covariates_for_sampling,
      iter = iter,
      save_data = save_data
    )
    
    # Run survey-weighted synthetic control
    weighted_results[[iter]] <- run_synthetic_control_simul(
      df.simul = df.simul, T0 = T0, use_weights = TRUE, target_district = target_district
    )
    
    # Run unweighted (raw-data) synthetic control
    unweighted_results[[iter]] <- run_synthetic_control_simul(
      df.simul, T0, use_weights = FALSE, target_district = target_district
    )
  }
  
  # Extract results
  weighted_p_values <- sapply(weighted_results, function(x) x$p_value)
  unweighted_p_values <- sapply(unweighted_results, function(x) x$p_value)
  
  weighted_effects <- sapply(weighted_results, function(x) x$effect)
  unweighted_effects <- sapply(unweighted_results, function(x) x$effect)
  
  weighted_target_ctrl_weights <- sapply(weighted_results, function(x) x$target_ctrl_weight)
  unweighted_target_ctrl_weights <- sapply(unweighted_results, function(x) x$target_ctrl_weight)
  
  # Calculate rejection rates
  type1_error <- 1 / num_regions
  weighted_reject_rate <- mean(weighted_p_values <= type1_error, na.rm = TRUE)
  unweighted_reject_rate <- mean(unweighted_p_values <= type1_error, na.rm = TRUE)
  
  return(list(
    weighted = list(
      p_values = weighted_p_values,
      estimated_effects = weighted_effects,
      reject_rate = weighted_reject_rate,
      target_ctrl_weights = weighted_target_ctrl_weights
    ),
    unweighted = list(
      p_values = unweighted_p_values,
      estimated_effects = unweighted_effects,
      reject_rate = unweighted_reject_rate,
      target_ctrl_weights = unweighted_target_ctrl_weights
    )
  ))
}

# Run simulation across all parameter combinations
run_full_simulation <- function(max_iter = 2,
                                beta_values = c(0, 2),
                                dissimilarity_factors = c(0, 2),
                                target_district = 2,
                                save_data = TRUE,
                                T0 = 40, num_regions = 21) {
  
  results <- list()

  for(sampling_type in c(FALSE, TRUE)) { # FALSE = constant, TRUE = covariate dependent
    sampling_label <- ifelse(sampling_type, "covariate", "constant")
    
    for(case_type in c(0, 1, 2)) { # 0 = baseline, 1 = similar, 2 = different
      for(beta in beta_values) {
        for(factor in dissimilarity_factors) {
          
          sim_results <- run_simulation(
            max_iter = max_iter,
            beta_ijt = beta,
            case_type = case_type,
            dissimilarity_factor = factor,
            target_district = target_district,
            use_covariates_for_sampling = sampling_type,
            save_data = save_data,
            num_regions = num_regions,
            T0 = T0
          )
          
          # Create a unique key
          config_key <- paste0(
            "case", case_type, "_",
            "beta_", beta, "_",
            "dissim_", factor, "_",
            sampling_label
          )
          
          results[[config_key]] <- list(
            params = list(
              case_type = case_type,
              beta = beta,
              dissimilarity_factor = factor,
              sampling_label = sampling_label
            ),
            weighted = list(
              p_values = sim_results$weighted$p_values,
              estimated_effects = sim_results$weighted$estimated_effects,
              reject_rate = sim_results$weighted$reject_rate,
              target_ctrl_weights = sim_results$weighted$target_ctrl_weights
            ),
            unweighted = list(
              p_values = sim_results$unweighted$p_values,
              estimated_effects = sim_results$unweighted$estimated_effects,
              reject_rate = sim_results$unweighted$reject_rate,
              target_ctrl_weights = sim_results$unweighted$target_ctrl_weights
            )
          )
        }
      }
    }
  }
  
  return(results)
}