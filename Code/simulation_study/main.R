library(Rlab)
library(dplyr)
library(data.table)
library(tidysynth)
library(stats)
library(ggplot2)
library(tidyr)
library(knitr)
library(grid)

set.seed(123456)

#setwd("~/Desktop/scm_indiv_survey")
source("Code/simulation_study/data_gen.R")
source("Code/simulation_study/apply_scm.R")
source("Code/simulation_study/sim.R")
source("Code/simulation_study/visual.R")

# Define simulation parameters
max_iter <- 10          # Number of Monte Carlo iterations (use 300 for real study)
beta_values <- c(0, 1, 2) # Individual treatment effect values (use seq(0, 2, by = 0.5) for real study)
dissimilarity_factors <- c(0, 1, 2)  # Dissimilarity factors (modify for real study)
target_district <- 2   # Target control for dissimilarity case studies (see data_generation.R for details)
num_regions <- 21      # Number of regions (1 treated, 20 control)
T0 <- 40               # Pre-intervention periods

simulation_results <- run_full_simulation(
  num_regions = num_regions,
  T0 = T0,
  beta_values = beta_values,               
  dissimilarity_factors = dissimilarity_factors,  
  target_district = target_district,   
  max_iter = max_iter, 
  save_data = TRUE                         
)

# Create plots
power_analysis <- create_power_analysis(simulation_results)
weight_analysis <- create_weight_analysis(simulation_results)
ggsave("Results/simulation_study/power_curves.png", power_analysis$plot, width = 10, height = 6)
ggsave("Results/simulation_study/weight_plots.png", weight_analysis$plot, width = 12, height = 8)

# Print tables
power_table <- power_analysis$table %>%
  kable(col.names = c("Individual Treatment Effect", "Sampling Type", 
                      "Raw Data Method", "Survey-Weighted Method"),
        format = "markdown", digits = 2,
        caption = "Rejection Rates")
print(power_table)

weight_table <- weight_analysis$table %>%
  kable(col.names = c("Case Type", "Dissimilarity Factor", "Sampling Type", 
                      "Raw Data Method", "Survey-Weighted Method"),
        format = "markdown", digits = 2,
        caption = "Target Control Weights")
print(weight_table)
