library(dplyr)
library(data.table)
library(tidysynth)
library(stats)
library(ggplot2)
library(Hmisc)
set.seed(123456)

#setwd("~/Desktop/scm_indiv_survey")
SADCQ <- fread("Data/data_application/SADCQ.csv")

#-------------------------------------------------
# Preprocess Data
#-------------------------------------------------

# Define time parameters
year.start <- 2007
year.end <- 2017 # year.trt

# Select districts for analysis
select.site <- c("Broward County, FL (FT)",
                 "Los Angeles, CA (LO)", "Miami-Dade County, FL (MM)",
                 "Orange County, FL (OL)", "Palm Beach County, FL (PB)", "Philadelphia, PA (PH)",
                 "San Diego, CA (SA)")

# Filter data and recode variables
df <- SADCQ %>%
  filter(sitename %in% select.site, year >= year.start, year <= year.end) %>%
  mutate(
    age = age + 11,  
    soda = case_when(
      q75 == 1 ~ 0,
      q75 == 2 ~ 2,
      q75 == 3 ~ 5,
      q75 == 4 ~ 1*7,
      q75 == 5 ~ 2*7,
      q75 == 6 ~ 3*7,
      q75 == 7 ~ 4*7
    ),
    overweight = case_when(
      qnowt == 1 ~ 1,
      qnowt == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    obesity = case_when(
      qnobese == 1 ~ 1,
      qnobese == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    female = case_when(
      sex == 1 ~ 1,  # Female
      sex == 2 ~ 0,  # Male
      TRUE ~ NA_real_
    ),
    black = case_when(
      race4 == 2 ~ 1, 
      TRUE ~ 0
    ),
    hispanic = case_when(
      race4 == 3 ~ 1,  
      TRUE ~ 0
    ),
    white = case_when(
      race4 == 1 ~ 1,  
      TRUE ~ 0
    ),
    other_race = case_when(
      race4 == 4 ~ 1,  
      TRUE ~ 0
    )
  )

#-------------------------------------------------
# Create Summary Table
#-------------------------------------------------

# Group separation
df_phil <- df %>% filter(sitename == "Philadelphia, PA (PH)")
df_ctrl <- df %>% filter(sitename != "Philadelphia, PA (PH)")

# Sample sizes
n_vals <- list(
  phil_weighted = round(sum(df_phil$weight, na.rm = TRUE)),
  ctrl_weighted = round(sum(df_ctrl$weight, na.rm = TRUE)),
  phil_unweighted = nrow(df_phil),
  ctrl_unweighted = nrow(df_ctrl)
)

summary_stats <- function(df, var, weighted = FALSE) {
  m <- if (weighted) weighted.mean(df[[var]], df$weight, na.rm = TRUE) else mean(df[[var]], na.rm = TRUE)
  
  if (grepl("female|black|hispanic|white|other_race|overweight|obesity", var)) {
    # Format as percentage
    return(sprintf("%.2f", m * 100))
  } else {
    # Continuous: return mean (SD)
    s <- if (weighted) sqrt(wtd.var(df[[var]], df$weight, na.rm = TRUE)) else sd(df[[var]], na.rm = TRUE)
    return(sprintf("%.2f (%.2f)", m, s))
  }
}

# Variables to dplyr::summarize
vars <- c("age", "female", "black", "hispanic", "white", "other_race", "bmi", "overweight", "obesity")
labels <- c("Age (SD)", "Female (%)", "Black (%)", "Hispanic (%)", "White (%)", "Other Race (%)",
            "BMI (SD)", "Overweight (%)", "Obesity (%)")

biennial_years <- seq(year.start, year.end, by = 2)

soda_rows <- do.call(rbind, lapply(biennial_years, function(y) {
  phil_w <- df_phil %>% filter(year == y)
  ctrl_w <- df_ctrl %>% filter(year == y)
  
  data.frame(
    Characteristic = paste0("Soda (", y, ")"),
    `Philadelphia (n = ` = summary_stats(phil_w, "soda", TRUE),
    `Controls (n = ` = summary_stats(ctrl_w, "soda", TRUE),
    `Philadelphia (n = ` = summary_stats(phil_w, "soda", FALSE),
    `Controls (n = ` = summary_stats(ctrl_w, "soda", FALSE),
    stringsAsFactors = FALSE
  )
}))

other_rows <- data.frame(
  Characteristic = labels,
  `Philadelphia (n = ` = sapply(vars, function(v) summary_stats(df_phil, v, TRUE)),
  `Controls (n = ` = sapply(vars, function(v) summary_stats(df_ctrl, v, TRUE)),
  `Philadelphia (n = ` = sapply(vars, function(v) summary_stats(df_phil, v, FALSE)),
  `Controls (n = ` = sapply(vars, function(v) summary_stats(df_ctrl, v, FALSE)),
  stringsAsFactors = FALSE
)

table_data <- rbind(soda_rows, other_rows)

colnames(table_data) <- c(
  "Characteristic",
  paste0("Philadelphia (n = ", n_vals$phil_weighted, ")"),
  paste0("Controls (n = ", n_vals$ctrl_weighted, ")"),
  paste0("Philadelphia (n = ", n_vals$phil_unweighted, ")"),
  paste0("Controls (n = ", n_vals$ctrl_unweighted, ")")
)

print(table_data)

#-------------------------------------------------
# Plot soda consumption trend
#-------------------------------------------------

# Create unweighted (raw-data) summary
df_mean_soda_unweighted <- df %>%
  group_by(sitename, year) %>%
  dplyr::summarize(mean_soda = mean(soda, na.rm = TRUE), .groups = "drop") %>%
  mutate(type = "Unweighted")

# Create survey-weighted summary
df_mean_soda_weighted <- df %>%
  group_by(sitename, year) %>%
  dplyr::summarize(mean_soda = weighted.mean(soda, weight, na.rm = TRUE), .groups = "drop") %>%
  mutate(type = "Survey-weighted")

df_mean_soda_all <- bind_rows(df_mean_soda_unweighted, df_mean_soda_weighted) %>%
  rename(site = sitename)

# Plot: weighted vs. unweighted soda consumption
soda_trend_plot <- ggplot(df_mean_soda_all, 
       aes(x = year, y = mean_soda, 
           color = site, linetype = type, 
           group = interaction(site, type))) +
  geom_line(linewidth = 1.2) + 
  ggtitle("Mean Servings of Soda per Week (2007–2017)") +
  ylab("Servings of soda per week (mean)") +
  xlab("Year") +
  scale_x_continuous(limits = c(year.start, year.end), 
                     breaks = seq(year.start, year.end, 2)) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
  geom_vline(xintercept = year.end, linetype = "dotted", 
             color = "black", linewidth = 1.2) +
  annotate("text", x = 2013, y = 7.5, 
           label = "Philadelphia beverage tax,\nJan 1, 2017", 
           hjust = 0, size = 4) +
  labs(color = "District", linetype = "Population") +
  scale_linetype_discrete(labels = c("sampled population", "survey weighted population")) +
  theme_minimal(base_size = 14)
ggsave("Results/data_application/soda_trend_plot.png", soda_trend_plot, width = 10, height = 6)

save(df, year.start, year.end, file = "Data/data_application/processed_data.RData")