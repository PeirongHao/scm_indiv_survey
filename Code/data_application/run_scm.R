library(dplyr)
library(data.table)
library(tidysynth)
library(stats)
library(ggplot2)
library(Hmisc)
set.seed(123456)

load("Data/data_application/processed_data.RData")

#-------------------------------------------------
# Compute region averages with/without survey weights
#-------------------------------------------------

# Compute control region averages with or without survey weights
ctrl_avg <- function(df, use_weights = FALSE) {
  df.ctrl <- df %>% filter(sitename != "Philadelphia, PA (PH)") %>%
    group_by(sitename, year) %>%
    dplyr::summarize(
      age = ifelse(use_weights, weighted.mean(age, weight, na.rm = TRUE),
                   mean(age, na.rm = TRUE)),
      bmi = ifelse(use_weights, weighted.mean(bmi, weight, na.rm = TRUE),
                   mean(bmi, na.rm = TRUE)),
      overweight = ifelse(use_weights, weighted.mean(overweight, weight, na.rm = TRUE), 
                          mean(overweight, na.rm = TRUE)),
      obesity = ifelse(use_weights, weighted.mean(obesity, weight, na.rm = TRUE), 
                       mean(obesity, na.rm = TRUE)),
      female = ifelse(use_weights, weighted.mean(female, weight, na.rm = TRUE),
                      mean(female, na.rm = TRUE)),
      black = ifelse(use_weights, weighted.mean(black, weight, na.rm = TRUE),
                     mean(black, na.rm = TRUE)),
      hispanic = ifelse(use_weights, weighted.mean(hispanic, weight, na.rm = TRUE),
                        mean(hispanic, na.rm = TRUE)),
      white = ifelse(use_weights, weighted.mean(white, weight, na.rm = TRUE),
                     mean(white, na.rm = TRUE)),
      other_race = ifelse(use_weights, weighted.mean(other_race, weight, na.rm = TRUE),
                          mean(other_race, na.rm = TRUE)),
      soda = ifelse(use_weights, weighted.mean(soda, weight, na.rm = TRUE), 
                    mean(soda, na.rm = TRUE)),
      .groups = 'drop'
    )
  
  return(df.ctrl)
}

# Compute treated region averages with survey weights
trt_avg <- function(df, df.ctrl) {
  df.trt <- df %>% filter(sitename == "Philadelphia, PA (PH)") %>%
    group_by(year) %>%
    dplyr::summarize(
      sitename = "Philadelphia, PA (PH)",
      age = weighted.mean(age, weight, na.rm = TRUE),
      bmi = weighted.mean(bmi, weight, na.rm = TRUE),
      overweight = weighted.mean(overweight, weight, na.rm = TRUE),
      obesity = weighted.mean(obesity, weight, na.rm = TRUE),
      female = weighted.mean(female, weight, na.rm = TRUE),
      black = weighted.mean(black, weight, na.rm = TRUE),
      hispanic = weighted.mean(hispanic, weight, na.rm = TRUE),
      white = weighted.mean(white, weight, na.rm = TRUE),
      other_race = weighted.mean(other_race, weight, na.rm = TRUE),
      soda = weighted.mean(soda, weight, na.rm = TRUE),
      .groups = 'drop'
    )
  
  return(bind_rows(df.ctrl, df.trt))
}

#-------------------------------------------------
# Run SCM
#-------------------------------------------------

# 1. Raw Data Approach
df_raw <- ctrl_avg(df, use_weights = FALSE)
df_raw_all <- trt_avg(df, df_raw)

raw_scm <- df_raw_all %>%
  synthetic_control(
    outcome = soda,
    unit = sitename,
    time = year,
    i_unit = "Philadelphia, PA (PH)",
    i_time = (year.end-2),
    generate_placebos = TRUE
  ) %>%
  generate_predictor(
    time_window = year.start:(year.end-2),
    age = mean(age, na.rm = TRUE),
    bmi = mean(bmi, na.rm = TRUE),
    overweight = mean(overweight, na.rm = TRUE),
    obesity = mean(obesity, na.rm = TRUE),
    female = mean(female, na.rm = TRUE),
    black = mean(black, na.rm = TRUE),
    hispanic = mean(hispanic, na.rm = TRUE),
    white = mean(white, na.rm = TRUE),
    other_race = mean(other_race, na.rm = TRUE)
  ) %>%
  generate_predictor(time_window = 2007, soda_2007 = soda) %>%
  generate_predictor(time_window = 2009, soda_2009 = soda) %>%
  generate_predictor(time_window = 2011, soda_2011 = soda) %>%
  generate_predictor(time_window = 2013, soda_2013 = soda) %>%
  generate_predictor(time_window = (year.end-2), soda_2015 = soda) %>%
  generate_weights(
    optimization_window = year.start:(year.end-2),
    optimization_method = "L-BFGS-B",
    margin_ipop = .01,
    sigf_ipop = 7,
    bound_ipop = 6
  ) %>%
  generate_control()

# 2. Survey-Weighted Approach
df_weighted <- ctrl_avg(df, use_weights = TRUE)
df_weighted_all <- trt_avg(df, df_weighted)

weighted_scm <- df_weighted_all %>%
  synthetic_control(
    outcome = soda,
    unit = sitename,
    time = year,
    i_unit = "Philadelphia, PA (PH)",
    i_time = (year.end-2),
    generate_placebos = TRUE
  ) %>%
  generate_predictor(
    time_window = year.start:(year.end-2),
    age = mean(age, na.rm = TRUE),
    bmi = mean(bmi, na.rm = TRUE),
    overweight = mean(overweight, na.rm = TRUE),
    obesity = mean(obesity, na.rm = TRUE),
    female = mean(female, na.rm = TRUE),
    black = mean(black, na.rm = TRUE),
    hispanic = mean(hispanic, na.rm = TRUE),
    white = mean(white, na.rm = TRUE),
    other_race = mean(other_race, na.rm = TRUE)
  ) %>%
  generate_predictor(time_window = 2007, soda_2007 = soda) %>%
  generate_predictor(time_window = 2009, soda_2009 = soda) %>%
  generate_predictor(time_window = 2011, soda_2011 = soda) %>%
  generate_predictor(time_window = 2013, soda_2013 = soda) %>%
  generate_predictor(time_window = (year.end-2), soda_2015 = soda) %>%
  generate_weights(
    optimization_window = year.start:(year.end-2),
    optimization_method = "L-BFGS-B",
    margin_ipop = .01,
    sigf_ipop = 7,
    bound_ipop = 6
  ) %>%
  generate_control()

#-------------------------------------------------
# Analyze Results
#-------------------------------------------------

# Plot results
raw_trend <- raw_scm %>% plot_trends() + ggtitle("Raw Data Approach: Observed vs. Synthetic Soda Consumption")
weighted_trend <- weighted_scm %>% plot_trends() + ggtitle("Survey-Weighted Approach: Observed vs. Synthetic Soda Consumption")

ggsave("Results/data_application/raw_scm_trend.png", raw_trend, width = 10, height = 6)
ggsave("Results/data_application/weighted_scm_trend.png", weighted_trend, width = 10, height = 6)

# Extract weights
raw_weights <- raw_scm %>% plot_weights()
weighted_weights <- weighted_scm %>% plot_weights()

ggsave("Results/data_application/raw_scm_weights.png", raw_weights, width = 10, height = 6)
ggsave("Results/data_application/weighted_scm_weights.png", weighted_weights, width = 10, height = 6)

# Extract treatment effects
raw_effect <- raw_scm %>% grab_synthetic_control(placebo = TRUE) %>%
  filter(.placebo == 0, time_unit == 2017) %>%
  dplyr::summarize(effect = real_y - synth_y) %>%
  pull(effect)

weighted_effect <- weighted_scm %>% grab_synthetic_control(placebo = TRUE) %>%
  filter(.placebo == 0, time_unit == 2017) %>%
  dplyr::summarize(effect = real_y - synth_y) %>%
  pull(effect)

# Extract p-values
raw_p <- raw_scm %>% grab_significance() %>%
  filter(type == 'Treated') %>%
  pull(fishers_exact_pvalue) 

weighted_p <- weighted_scm %>% grab_significance() %>%
  filter(type == 'Treated') %>%
  pull(fishers_exact_pvalue) 

# Create results table
results_df <- data.frame(
  Method = c("Raw Data Approach", "Survey-Weighted Approach"),
  `Treatment Effect` = c(raw_effect, weighted_effect),
  `P-value` = c(raw_p, weighted_p)
)

print(results_df)