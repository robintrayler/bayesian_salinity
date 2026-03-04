# This script validates the re-code of the Kim et al. (2014) model using the 
# arctic sharks data set. 
# load required packages and data ---------------------------------------------
library(tidyverse)
source('./R/functions/bayesian_salinity.R')
theme_set(theme_classic())

# load the data ---------------------------------------------------------------
# load the salinity regression data from Kim et al. (2014)

salinity_regression_data <- read.csv(file = './data/salinity_regression_data.csv')

teeth        <- read.csv(file = './data/sharks/Carcharias taurus.csv', header = TRUE) 
d18Ofw_prior <- read.csv(file = './data/iCESM_x3CO2_d18Of.csv')
t_prior      <- read.csv(file ='./data/iCESM_x3CO2_TEMP_d18Osw.csv')

# set up data 
d18Op = teeth$d18Op

# set up vectors to use as prior distributions
salinity_prior_vector    <- runif(100000, 0, 40)
temperature_prior_vector = runif(100000, 12, 16)
d18Ofw_prior_vector = runif(100000, -8, -6.9)

# run the Bayesian model ------------------------------------------------------
parameters <- bayesian_salinity(
  d18Op = d18Op,
  water_eq = 'lecuyer',
  salinity_regression_data = salinity_regression_data,
  temperature_prior_vector = temperature_prior_vector,
  d18Ofw_prior_vector = d18Ofw_prior_vector,
  salinity_prior_vector = salinity_prior_vector,
  iterations = 100000)

# # plot histograms of the results ----------------------------------------------
# parameters$parameters |>
#   pivot_longer(cols = c('B0', 'B1', 'sigma_reg',
#                         'temperature',
#                         'salinity',
#                         'd18Ofw',
#                         'd18Oenv',
#                         'sigma_env')) |>
#   filter(iteration > parameters$burn) |>
#   mutate(name = factor(name, levels = c('B0', 'B1', 'sigma_reg',
#                                         'temperature',
#                                         'salinity',
#                                         'd18Ofw',
#                                         'd18Oenv',
#                                         'sigma_env'))) |>
#   ggplot(mapping = aes(x = value,
#                        fill = name)) +
#   geom_histogram(bins = 30,
#                  show.legend = FALSE) +
#   facet_wrap(name~.,
#              scales = 'free',
#              nrow = 2)
# 
# # take a look at trace plots --------------------------------------------------
# parameters$parameters |>
#   pivot_longer(cols = c('B0', 'B1', 'sigma_reg',
#                         'temperature',
#                         'salinity',
#                         'd18Ofw',
#                         'd18Oenv',
#                         'sigma_env')) |>
#   mutate(name = factor(name, levels = c('B0', 'B1', 
#                                         'sigma_reg',
#                                         'temperature',
#                                         'salinity',
#                                         'd18Ofw',
#                                         'd18Oenv',
#                                         'sigma_env'))) |>
#   ggplot(mapping = aes(x = iteration,
#                        y = value,
#                        color = name)) +
#   geom_line(show.legend = FALSE) +
#   facet_wrap(name~.,
#              scales = 'free',
#              nrow = 2)

# write the results to a file -------------------------------------------------
parameters$parameters |> 
  add_column(taxon = 'Carcharias taurus (modern)') |> 
  filter(iteration > parameters$burn) |> 
  select(-iteration) |> 
  write_csv(file = './results/Carcharias taurus_posterior.csv')

parameters$CI |> 
  add_column(taxon = 'Striatolamia') |> 
  write_csv(file = './results/Carcharias taurus_CI.csv')