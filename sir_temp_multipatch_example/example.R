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

df_tidy <- df_tidy %>%
  pivot_wider(names_from = population, values_from=cases, names_prefix = "cases_")


# turn into data object
data <- mcstate::particle_filter_data(df_tidy, 
                                      time = "day",
                                      rate = 1,
                                      initial_time = 0)

plot(cases_A ~ day, df_tidy,
     type = "o", xlab = "Day", ylab = "New cases", pch = 19)
lines(cases_B ~ day, df_tidy,
      type = "o", xlab = "Day", ylab = "New cases", pch = 19, col = 2)
legend("topright", col = 2:1, legend = c("A", "B"), lwd = 1)

# sir model
sir <- odin.dust::odin_dust("sir_temp_multipatch_example/sir_temp_multipatch.R", verbose=FALSE) # standard SIR model

# look at the information
sir_model <- sir$new(pars = list(I0=c(10,10), S0=c(1000,1000), Temp=c(20,24.5)), time = 0, n_particles = 1L)
sir_model$info()

# run and view
model <- sir$new(pars = list(N_patch=2,
                                 Temp=c(20,24.5),
                                 temp_scale=1e-3,
                                 I0=c(10,10),
                                 S0=c(pop_size,pop_size)),
                     time = 1,
                     n_particles = 10,
                     n_threads = 1L,
                     seed = 1L)

#Define how long the model runs for, number of time steps
n_times <- 100

x <- array(NA, dim = c(model$info()$len, 10, n_times))

# For loop to run the model iteratively
for (t in seq_len(n_times)) {
  x[ , , t] <- model$run(t)
}
time <- x[4, 1, ]
x <- x[-c(1:4), , ]

# Plotting the trajectories
par(mfrow = c(1,2), oma=c(2,3,0,0))
for (i in 1:2) {
  par(mar = c(3, 4, 2, 0.5))
  cols <- c(S = "#8c8cd9", I = "#cc0044", R = "#999966")
  matplot(time, t(x[i,, ]), type = "l", # Offset to access numbers in age compartment
          xlab = "", ylab = "", yaxt="none", 
          col = cols[["S"]], lty = 1)
  matlines(time, t(x[i + 2, , ]), col = cols[["I"]], lty = 1)
  matlines(time, t(x[i + 2*2, , ]), col = cols[["R"]], lty = 1)
  legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
  axis(2, las =2)
}
mtext("Number of individuals", side=2,line=1, outer=T)
mtext("Time", side = 1, line = 0, outer =T)

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

compare <- function(state, observed, pars=NULL){ ## NOTE: this is not nicely written- I will come back to this
  l1 <- ll_pois(observed$cases_A,state[1, , drop = TRUE])
  l2 <- ll_pois(observed$cases_B,state[2, , drop = TRUE])
  
  l1+l2
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
  mcstate::pmcmc_parameter("temp_scale", 1e-3, min=0))

transform <- function(theta) {
  as.list(theta)
} # this is where you may include a log transform or other

make_transform <- function(I0, S0, Temp) {
  function(theta) {
    list(I0 =I0, #this runs everything with standardised I0 inital infecion
         S0=S0,
         Temp=Temp,
         beta=theta[["beta"]],
         gamma=theta[["gamma"]],
         temp_scale=theta[["temp_scale"]])
  }
}

pars <- list(beta=0.25, gamma=0.2, temp_scale = 1e-5)
no_param <- length(pars)

vcv <- (1e-2) ^ 2 * diag(no_param) / no_param # this is for the proposal distribution
transform <- make_transform(I0=c(10,10), S0=c(pop_size, pop_size), Temp=c(20, 24.5)) #then this bounds I0 as 10

#getting proposal priors and starting parameters in the correct format
mcmc_pars <- mcstate::pmcmc_parameters$new(priors, vcv, transform) 

n_steps_in <- 1e4

control <- mcstate::pmcmc_control(
  n_steps = n_steps_in,
  progress = TRUE, save_trajectories = TRUE)

samples <- mcstate::pmcmc(mcmc_pars, filter, control = control)

par(mfrow = c(1,1))
#samples
plot(samples$probabilities[500:n_steps_in, "log_posterior"], type = "s",
     xlab = "Sample", ylab = "Log posterior") #this shows the chain


#### view the output
# run and view
model <- sir$new(pars = list(N_patch=2,
                             Temp=c(20,24.5),
                             temp_scale=samples$pars[1000,3],
                             beta = samples$pars[1000,1],
                             gamma = samples$pars[1000,2],
                             I0=c(10,10),
                             S0=c(pop_size,pop_size)),
                 time = 1,
                 n_particles = 10,
                 n_threads = 1L,
                 seed = 1L)

#Define how long the model runs for, number of time steps
n_times <- 100

x <- array(NA, dim = c(model$info()$len, 10, n_times))

# For loop to run the model iteratively
for (t in seq_len(n_times)) {
  x[ , , t] <- model$run(t)
}
time <- x[4, 1, ]
A <- x[13,,]
B <- x[14,,]

# Plotting the trajectories
par(mfrow = c(1,2))
matplot(time, t(A), type="l", col = "grey80", ylim = range(df_tidy$cases_A))
points(time, df_tidy$cases_A, col="red")

matplot(time, t(B), type="l", col = "grey80", ylim=range(df_tidy$cases_B))
points(time, df_tidy$cases_B, col="red")
