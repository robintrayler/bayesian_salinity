bayesian_salinity <- function(
    d18Op,
    temperature_prior_vector = runif(100000, 0, 45),
    d18Ofw_prior_vector = runif(100000, -40, 0),
    salinity_prior_vector = runif(100000, 0, 50),
    salinity_regression_data,
    water_eq = 'longinelli',
    iterations = 100000,
    burn = round(iterations * 0.2, 0)) {
  
  # load libraries and helper functions
  library(tidyverse)
  library(progress)
  source('./R/functions/adaptive_update.R')
  source('./R/functions/truncated_random_normal.R')
  source('./R/functions/truncated_standard_normal.R')
  
  # INPUTS: 
  # d18Op = vector of phosphate oxygen isotope compositions from shark teeth
  # temperature_prior_vector = vector of empirical temperature data to use as a 
  # prior distribution. Defaults to a uniform distribution between 
  # 0 and 45C
  # d18Ofw_prior_vector = vector of empirical freshwater d18O values to use as a 
  # prior distribution. Defaults to a uniform distribution between -40 
  # and 0 permil.
  # salinity_prior_vector = vector of emperical salinity values to use as a prior 
  # distribution. Defaults to a uniform distribution between 0 and 50 PSU.
  # salinity_regression_data = empirical distribution of d18Ofw intercepts and
  # Δd18O / Δsalinity slopes. Defaults to the data from Kim et al. (2014). 
  # columns must be named exactly "slope", "intercept" 
  # water_eq = equation to use for converting d18Op to water values. Options are
  # "longinelli", "kolodny", "puceat". Defaults to "longinelli"
  
  # iterations = number of MCMC iterations to run for
  # burn = number of initial iterations to discard. defaults to 20% of iterations
  
  # OUTPUTS: 
  # outputs a list containing the following objects
  # parameters: data.frame containing raw posterior distributions (including burn-in)
  # B1 = posterior distribution of  slope for salinity_regression_data
  # B0 = posterior distribution of intercept for salinity_regression_data
  # sigma_reg = posterior distribution of dispersion parameter 
  # for salinity_regression_data.
  
  # temperature = posterior distribution of temperature
  # salinity = posterior distribution of salinity
  # d18Ofw = posterior distribution of the d18O of freshwater
  # d18Oenv = posterior distribution of the d18O of the environmental waters
  # sigma_env = posterior distribution of the dispersion parameter associated with
  # the environmental parameters
  # CI: data.frame containing 95% credible interval for `parameters` with the burn-in period removed
  # iteration and burn are also returned to help with data manipulation
  
  # set up parameter storgage ---------------------------------------------------
  parameters <- data.frame(temperature = vector(length = iterations),
                           salinity = vector(length = iterations),
                           d18Ofw = vector(length = iterations),
                           d18Oenv = vector(length = iterations),
                           sigma_env = vector(length = iterations),
                           B0 = vector(length = iterations),
                           B1 = vector(length = iterations),
                           sigma_reg = vector(length = iterations),
                           iteration = vector(length = iterations))
  
  
  # initialize parameters 
  parameters$temperature[1] <- runif(1, 
                                     min(temperature_prior_vector), 
                                     max(temperature_prior_vector))
  parameters$salinity[1]    <- runif(1, 
                                     min(salinity_prior_vector),
                                     max(salinity_prior_vector))
  parameters$d18Ofw[1]      <- runif(1, 
                                     min(d18Ofw_prior_vector), 
                                     max(d18Ofw_prior_vector))
  parameters$sigma_env[1]   <- runif(1, 0, 1)
  
  # start regression parameters with a linear regression
  coefs <- salinity_regression_data |> 
    lm(intercept ~ slope, data = _) |> 
    coef()
  parameters$B1[1]          <- coefs[2]
  parameters$B0[1]          <- coefs[1]
  parameters$sigma_reg[1]   <- runif(1, 0, 1)
  
  # start iteration counter at 1
  parameters$iteration[1]   <- 1
  
  # set up prior distribution functions for temperature and d18Ow
  t_dns <- density(temperature_prior_vector)
  t_prior_f <- approxfun(t_dns$x, t_dns$y)
  
  d18Ofw_dns <- density(d18Ofw_prior_vector)
  d18Ofw_prior_f <- approxfun(d18Ofw_dns$x, d18Ofw_dns$y)
  
  salinity_dns <- density(salinity_prior_vector)
  salinity_prior_f <- approxfun(salinity_dns$x, salinity_dns$y)
  
  # set up a progress bar
  pb <- progress_bar$new(total = iterations,
                         format = '[:bar] :percent eta: :eta')
  
  
  calculate_d18Op <- function(d18Oenv, temperature, water_eq) {
    switch(water_eq,
           'longinelli' = d18Oenv - (temperature - 111.4) / 4.3,
           'kolodny'    = d18Oenv - (temperature - 113.3) / 4.38,
           'puceat'     = d18Oenv - (temperature - 118.7) / 4.22,
           'lecuyer'    = d18Oenv - (temperature - 117.4) / 4.5)
  }
  
  # main model loop -----------------------------------------------------------
  for(i in 2:iterations) {
    pb$tick()
    # set old values to current iteration
    parameters[i, ] <- parameters[i - 1, ]
    parameters$iteration[i] <- i
    
    # update regression parameters first 
    B1        <- adaptive_update(chain = parameters$B1, i = i)
    B0        <- adaptive_update(chain = parameters$B0, i = i)
    sigma_reg <- adaptive_update(chain = parameters$sigma_reg, i = i, lower = 0)
    
    # calculate regression mean
    reg_mu_proposed <- salinity_regression_data$slope * B1 + B0
    reg_mu_current  <- salinity_regression_data$slope * parameters$B1[i - 1] + parameters$B0[i - 1]
    # calculate probabilities
    p_reg_proposed  <- 
      # probability of the data 
      sum( dnorm(salinity_regression_data$intercept, 
                 mean = reg_mu_proposed, 
                 sd = sigma_reg, 
                 log = TRUE) ) + 
      # prior probabilities using vague priors
      dnorm(B1, 
            mean = 0, 
            sd = 100,
            log = TRUE) + 
      dnorm(B0, 
            mean = 0,
            sd = 100,
            log = TRUE) + 
      dgamma(sigma_reg, 
             shape = 1,
             rate = 1,
             log = TRUE)
    
    
    p_reg_current  <- 
      # probability of the data 
      sum( dnorm(salinity_regression_data$intercept, 
                 mean = reg_mu_current, 
                 sd = parameters$sigma_reg[i - 1], 
                 log = TRUE) ) + 
      # prior probabilities using vague priors
      dnorm(parameters$B1[i - 1], 
            mean = 0, 
            sd = 100,
            log = TRUE) + 
      dnorm(parameters$B0[i - 1], 
            mean = 0,
            sd = 100,
            log = TRUE) + 
      dgamma(parameters$sigma_reg[i - 1], 
             shape = 1,
             rate = 1,
             log = TRUE)
    
    a <- p_reg_proposed - p_reg_current
    
    if(is.finite(a)) {
      if(!is.na(a)) {
        if(exp(a) > runif(1)) {
          parameters$B1[i]        <- B1
          parameters$B0[i]        <- B0
          parameters$sigma_reg[i] <- sigma_reg
          
        }
      }
    }
    
    # update the environmental parameters ---------------------------------------
    # update each parameter individually in a random order
    ord <- sample(1:4, 4, replace = FALSE)
    
    for(j in ord) {
     
       # propose new values
      if(j == 1) {
        temperature <- adaptive_update(chain = parameters$temperature, 
                                       i = i)
        salinity    <- parameters$salinity[i]
        d18Ofw      <- parameters$d18Ofw[i]
        sigma_env   <- parameters$sigma_env[i]
      } else if (j == 2) {
        salinity    <- adaptive_update(chain = parameters$salinity, 
                                       i = i)
        temperature <- parameters$temperature[i]
        d18Ofw      <- parameters$d18Ofw[i]
        sigma_env   <- parameters$sigma_env[i]
      } else if (j == 3) {
        d18Ofw      <- adaptive_update(chain = parameters$d18Ofw, 
                                       i = i)
        salinity    <- parameters$salinity[i]
        temperature <- parameters$temperature[i]
        sigma_env   <- parameters$sigma_env[i]
      } else if (j == 4) {
        sigma_env   <- adaptive_update(chain = parameters$sigma_env, 
                                       i = i, 
                                       lower = 0)
        salinity <- parameters$salinity[i-1]
        temperature <- parameters$temperature[i-1]
        d18Ofw <- parameters$d18Ofw[i-1]
        
      }
      
      # calculate mu proposed
      d18Oenv_proposed = d18Ofw + salinity * ((d18Ofw - B0) / B1)
      # convert to d18O phosphate for probability calculations
      
      mu_proposed <- calculate_d18Op(d18Oenv_proposed, temperature, water_eq)
      # mu_proposed <- d18Oenv_proposed - (temperature - 111.4) / 4.3
      
      # calculate mu current
      d18Oenv_current = parameters$d18Ofw[i - 1] + 
        parameters$salinity[i - 1] * ((parameters$d18Ofw[i - 1] - 
                                         parameters$B0[i]) / parameters$B1[i])
      
      # convert to d18O phosphate for probability calculations
      mu_current <- calculate_d18Op(d18Oenv_current, parameters$temperature[i - 1], water_eq)
      # mu_current <- d18Oenv_current - (parameters$temperature[i - 1] - 111.4) / 4.3
      
      # calculate probabilities
      p_env_proposed <- 
        # likelihood
        sum( dnorm(d18Op, 
                   mean = mu_proposed,
                   sd   = sigma_env,
                   log  = TRUE) ) + 
        # prior probability
        log(t_prior_f(temperature)) + 
        # dunif(temperature, min = 8,     max = 13,    log = TRUE) + 
        # dunif(salinity,    min = 0,     max = 50,    log = TRUE) +
        log(salinity_prior_f(salinity)) + 
        log(d18Ofw_prior_f(d18Ofw)) + 
        dgamma(sigma_env, shape = 1, rate = 1, log = TRUE)
      
      p_env_current <- 
        # likelihood
        sum( dnorm(d18Op, 
                   mean = mu_current,
                   sd   = parameters$sigma_env[i - 1],
                   log  = TRUE) ) + 
        # prior probability
        log(t_prior_f(parameters$temperature[i-1])) + 
        log(salinity_prior_f(parameters$salinity[i-1])) + 
        log(d18Ofw_prior_f(parameters$d18Ofw[i-1]) ) + 
        dgamma(parameters$sigma_env[i-1], shape = 1, rate = 1, log = TRUE)
      
      # calculate alpha
      a <- p_env_proposed - p_env_current
      
      if(is.finite(a)) {
        if(!is.na(a)) {
          if(exp(a) > runif(1)) {
            parameters$temperature[i] <- temperature
            parameters$salinity[i]    <- salinity
            parameters$d18Ofw[i]      <- d18Ofw
            parameters$sigma_env[i]   <- sigma_env
            parameters$d18Oenv[i]     <- d18Oenv_proposed
            
          }
        }
      }
    }
  }
  
  parameters_CI <- 
    parameters |> 
    filter(iteration > burn) |> 
    select(-iteration) |> 
    pivot_longer(cols = c(1:8)) |> 
    group_by(name) |> 
    reframe(quantiles = quantile(value, probs = c(0.025, 0.5, 0.975))) |> 
    ungroup() |> 
    add_column(quantile = rep(c('CI_2.5', 'CI_0.5', 'CI_97.5'), 8)) |> 
    pivot_wider(names_from = 'quantile',
                values_from = 'quantiles') |> 
    mutate_if(is.numeric, round, 1)
  
  return(list(parameters = parameters,
              CI = parameters_CI,
              iterations = iterations,
              burn = burn))
}




