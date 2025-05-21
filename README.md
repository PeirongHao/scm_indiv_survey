# Synthetic Control Methods for Policy Evaluation with Individual-Level Survey Data

## Overview

Synthetic control methods (SCM) offer a data-driven alternative to difference-in-differences (DiD) for policy evaluation when the parallel trends assumption may not hold. However, implementations of the conventional SCM often focus on aggregate-level data and may not perform well with individual-level data. Moreover, valid inference using probability samples from the target population remains unclear.

We propose to construct synthetic controls for intervention districts that include multiple individual units. Our methods consider individual-level outcomes and covariates, as well as district-level covariates, both with and without survey weighting. We conduct simulation studies under a linear factor model to evaluate the performance of these approaches. We apply the proposed methods to evaluate the impact of Philadelphia's sweetened beverage tax on soda consumption among high school students.

This repository provides the code and data to reproduce the results in the paper.

## Data

The data application examines the impact of Philadelphia's sweetened beverage tax—implemented in January 2017 at 1.5 cents per ounce—on high school students' soda consumption using YRBS data.

### Variables used in the analysis

- `q75`: Soda consumption
- `sitename`: District name
- `year`: Survey year
- `weight`: Survey weight
- `age`: Age
- `sex`: Sex
- `race4`: Race/ethnicity
- `bmi`: Body mass index (BMI)
- `qnowt`: Overweight indicator
- `qnobese`: Obesity indicator

## Code

### Data Application

- `data_preprocess.R`: filters data from selected districts (2007–2017), recodes variables, creates summary statistics, and produces trend plots for soda consumption.
- `run_scm.R`: calculates district-level averages with/without survey weights, fits SCM, and outputs treatment effect estimates and corresponding p-values.

### Simulation Study

- `main.R`: defines parameters and runs all simulations.
- `data_gen.R`: generates individual-level data based on a linear factor model.
- `apply_scm.R`: computes region-level averages and applies SCM.
- `sim_func.R`: executes iterative simulations across various parameter configurations.
- `visual_func.R`: produces plots and tables summarizing simulation results.
