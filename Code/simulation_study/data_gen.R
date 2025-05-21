# Case types: 
# 1: Treated region and target control are similar, different from others
# 2: Target control is different from treated and all other controls
# 0: No systematic differences 

generate_data_simul <- function(case_type = 0, dissimilarity_factor = 0, target_district = 2, 
                                use_covariates_for_sampling = FALSE, 
                                iter = 1, save_data = TRUE, beta_ijt = 0,
                                num_regions = 21, T0 = 40) { 
  # Validate case_type
  if (!case_type %in% c(0, 1, 2)) {
    stop("case_type must be either 0 (baseline), 1 (similar), or 2 (different)")
  }
  
  T_end <- T0 + 1         # Total periods (pre + post-intervention)
  
  # Generate region-level variables
  Z1 <- rnorm(num_regions, mean = 7, sd = 1)
  Z2 <- rnorm(num_regions, mean = 3, sd = 0.5)
  
  # Add dissimilarity based on case type
  if (case_type == 1) {
    Z1[c(1, target_district)] <- rnorm(2, mean = 7 + dissimilarity_factor, sd = 1)
    Z2[c(1, target_district)] <- rnorm(2, mean = 3 + dissimilarity_factor, sd = 0.5)
  } else if (case_type == 2) {
    Z1[target_district] <- rnorm(1, mean = 7 + dissimilarity_factor, sd = 1)
    Z2[target_district] <- rnorm(1, mean = 3 + dissimilarity_factor, sd = 0.5)
  }
  
  nu <- rnorm(num_regions, mean = 0, sd = 0.2)
  eta <- matrix(rnorm(num_regions * T_end, mean = 0, sd = 0.1), 
                nrow = num_regions, ncol = T_end)
  
  N <- matrix(runif(num_regions * T_end, min = 60, max = 100), 
              nrow = num_regions, ncol = T_end)
  
  # Time-varying parameters
  delta <- log(seq(1, T_end, length.out = T_end))
  theta1 <- 0.003 * seq(1, T_end)
  theta2 <- 0.01 * seq(T_end, 1)
  phi1 <- 0.05 * log(seq(1, T_end))
  phi2 <- 0.02 * seq(1, T_end)
  lambda <- seq(0, 1, length.out = T_end)
  psi <- seq(1, 0, length.out = T_end)
  
  df.simul <- list()
  
  for (j in 1:num_regions) {
    for (t in 1:T_end) {
      N_jt <- round(N[j, t])

      # Initialize individual-level variables with standard values
      X1.j <- rnorm(N_jt, mean = 5, sd = 1)
      X2.j <- rnorm(N_jt, mean = 6, sd = 0.5)
      
      # Modify X values based on case type and region
      if (case_type == 1 && (j == 1 || j == target_district)) {
        X1.j <- rnorm(N_jt, mean = 5 + dissimilarity_factor, sd = 1)
        X2.j <- rnorm(N_jt, mean = 6 + dissimilarity_factor, sd = 0.5)
      } else if (case_type == 2 && j == target_district) {
        X1.j <- rnorm(N_jt, mean = 5 + dissimilarity_factor, sd = 1)
        X2.j <- rnorm(N_jt, mean = 6 + dissimilarity_factor, sd = 0.5)
      }
      
      mu.j <- rnorm(N_jt, mean = 1, sd = 0.3)
      eps.j <- matrix(rnorm(N_jt * T_end, mean = 0, sd = 0.1), 
                      nrow = N_jt, ncol = T_end)
      
      for (i in 1:N_jt) {
        if (use_covariates_for_sampling) { # Covariate-dependent sampling probability
          alpha_level <- 0.01 * (1 + X1.j[i] + Z1[j] + X2.j[i] + Z2[j])
          prob <- 1 / (1 + exp(-alpha_level))
        } else { # Constant sampling probability
          prob <- 0.5
        }
        
        S.ij <- rbern(1, prob = prob) # Sampling indicator
        W.ij <- ifelse(S.ij == 1, 1 / prob, NA) # Sampling weight
        
        if (S.ij == 1) { # Generate outcome with possible treatment effect
          Y.ij <- delta[t] + theta1[t] * X1.j[i] + theta2[t] * X2.j[i] +
            phi1[t] * Z1[j] + phi2[t] * Z2[j] +
            lambda[t] * mu.j[i] + psi[t] * nu[j] +
            eps.j[i, t] + eta[j, t] + 
            beta_ijt * (j == 1 & t > T0)
          
          df.simul[[length(df.simul) + 1]] <- data.frame(
            region = j, year = t, id = i,
            Z1 = Z1[j], Z2 = Z2[j], 
            X1 = X1.j[i], X2 = X2.j[i],
            Y = Y.ij, S = S.ij, W = W.ij
          )
        }
      }
    }
  }
  
  df.simul <- rbindlist(df.simul)
  
  if (save_data) {
    sampling_type <- ifelse(use_covariates_for_sampling, "covariate", "constant")
    
    save_simulation_data(
      df.simul, 
      iter = iter, 
      dissimilarity_factor = dissimilarity_factor, 
      case_type = case_type, 
      sampling_type = sampling_type,
      beta_ijt = beta_ijt
    )
  }
  
  return(df.simul)
}

save_simulation_data <- function(df.simul, iter=1, dissimilarity_factor = 0, 
                                 case_type = 0, sampling_type = "constant",
                                 beta_ijt = 0) {
  # Create directory if it doesn't exist
  dir_name <- paste0("Data/simulation_study/", sampling_type)
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, recursive = TRUE)
  }
  
  filename <- paste0(dir_name, "/case_", case_type, "_", 
                     "beta_", beta_ijt, "_",
                     "dissim_", dissimilarity_factor, "_",
                     "iter_", iter, ".csv")
  
  # Save data
  fwrite(df.simul, filename)
  return(filename)
}
