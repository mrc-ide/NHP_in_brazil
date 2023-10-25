# a script for examining fitting a NHP in one region
# loading libraries
library(mcstate)
library(odin.dust)
library(dplyr)
library(tidyr)

#define the population size fo Sao Paulo
pop_dens <- 10 # based on density in Culot et al. Botacatu and others
area <- 248219 #based on Sao Paulo state

pop_size <- pop_dens * area #assume similar across the region *assumption* (we will assume this is equal in the 2 places)


# load example data
df_tidy <- read.csv(system.file("nested_sir_incidence.csv", package="mcstate"), 
                    stringsAsFactors = TRUE) # this is just example data from the package


# turn into data object
data <- mcstate::particle_filter_data(df_tidy, 
                                      time = "day",
                                      rate = 1,
                                      initial_time = 0,
                                      population="population")

plot(cases ~ day, df_tidy[df_tidy$population == "B", ],
     type = "o", xlab = "Day", ylab = "New cases", pch = 19)
lines(cases ~ day, df_tidy[df_tidy$population == "A", ],
      type = "o", xlab = "Day", ylab = "New cases", pch = 19, col = 2)
legend("topright", col = 2:1, legend = c("A", "B"), lwd = 1)

# sir model
sir <- odin.dust::odin_dust("sir_temp_multipatch_example/sir_temp_multipatch.R") # standard SIR model

# look at the information
sir_model <- sir$new(pars = list(I0=c(10,10), S0=c(1000,1000), Temp=c(20,24.5)), time = 0, n_particles = 1L)
sir_model$info()

# sourcing the compare function which defines what to fit to and how- here it is assuming a poisson likelihood
# sourcing the index function- this tells mcstate which model and data outputs to look for

# need a new index function and compare function
# Index
index <- function(info){
  list(run = c(cases = info$index$cases_inc),
       state = c(S = info$index$S,
                 I = info$index$I,
                 R = info$index$R))
}

# log-likelihood of Poisson count
ll_pois <- function(obs, model) {
  exp_noise <- 1e6
  if (is.na(obs)) {
    # Creates vector of zeros in ll with same length, if no data
    ll_obs <- numeric(length(model))
  } else {
    lambda <- model +
      rexp(n = length(model), rate = exp_noise)
    ll_obs <- dpois(x = obs, lambda = lambda, log = TRUE)
  }
  ll_obs
}

compare <- function(state, observed, pars=NULL){
  ll_pois(observed %>% tidyr::pivot_wider(names_from = population, values_from = cases, names_prefix = "cases_"), 
          state["cases", , drop = TRUE])
}

#parameteres
filter <- mcstate::particle_filter$new(data, 
                                       model = sir,
                                       n_particles=100,
                                       compare = compare,
                                       index = index) #just constructs the object
# 


# now for the fitting
priors <- list(
  mcstate::pmcmc_parameter("beta", 0.2, min = 0),
  mcstate::pmcmc_parameter("gamma", 0.1, min = 0, 
                           prior = function(p)  dgamma(p, shape = 1, scale = 0.2, log = TRUE)),
  mcstate::pmcmc_parameter("temp_scale", 1e-5, min=0))

transform <- function(theta) {
  as.list(theta)
} # this is where you may include a log transform or other

make_transform <- function(I0, Temp) {
  function(theta) {
    list(I0 =I0, #this runs everything with standardised I0 inital infecion
         Temp=Temp,
         beta=theta[["beta"]],
         gamma=theta[["gamma"]],
         temp_scale=theta[["temp_scale"]])
  }
}

pars <- list(beta=0.25, gamma=0.2, temp_scale = 1e-5)
no_param <- length(pars)

vcv <- (1e-2) ^ 2 * diag(no_param) / no_param # this is for the proposal distribution
transform <- make_transform(I0=c(10,10), Temp=c(20, 24.5)) #then this bounds I0 as 10

#getting proposal priors and starting parameters in the correct format
mcmc_pars <- mcstate::pmcmc_parameters$new(priors, vcv, transform) 

n_steps_in <- 1e3

control <- mcstate::pmcmc_control(
  n_steps = n_steps_in,
  progress = TRUE)

samples <- mcstate::pmcmc(mcmc_pars, filter, control = control)

#samples
plot(samples$probabilities[500:1000, "log_posterior"], type = "s",
     xlab = "Sample", ylab = "Log posterior") #this shows the chain

# plot result
pars <- list(beta = samples$pars[n_steps_in,1], gamma = samples$pars[n_steps_in,2], temp_scale=samples$pars[n_steps_in,3])
filter$run(pars, save_history = TRUE)
h <- filter$history()

matplot(h["t", 1, ], t(h["cases", , ]), type = "l", col = "#00000011", 
        xlab = "Day", ylab = "Cases", las = 1)
points(cases ~ day, df_tidy, pch = 19, col = "red")  
