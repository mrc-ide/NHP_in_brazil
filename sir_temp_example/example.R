# a script for examining fitting a NHP in one region
# loading libraries
library(mcstate)
library(odin.dust)


#define the population size fo Sao Paulo
pop_dens <- 10 # based on density in Culot et al. Botacatu and others
area <- 248219 #based on Sao Paulo state

pop_size <- pop_dens * area #assume similar across the region *assumption*

# let's examine temperature
mean_Temp <- 24.5 # mean temperature in degrees for Sao Paulo

# load example data
df_tidy <- read.csv(system.file("sir_incidence.csv", package="mcstate")) # this is just example data from the package

# turn into data object
data <- mcstate::particle_filter_data(df_tidy, 
                                      time = "day",
                                      rate = 1,
                                      initial_time = 0)

# sir model
sir <- odin.dust::odin_dust("sir_temp_example/sir_temp.R") # standard SIR model

# sourcing the compare function which defines what to fit to and how- here it is assuming a poisson likelihood
# sourcing the index function- this tells mcstate which model and data outputs to look for
source("sir_temp_example/compare_sir_temp.R") # observation process could also appear here
source("sir_temp_example/index_sir_temp.R")

#parameteres
pars <- list( gamma = 0.2, S0 = pop_size, Temp=mean_Temp, temp_scale=3e-6) #leave others as defaults
mod <- sir$new(pars, 0, max(data$time_end)) 
y <- mod$simulate(c(0, data$time_end))

filter <- mcstate::particle_filter$new(data, 
                                       model = sir, 
                                       n_particles=100,
                                       compare = compare, 
                                       index = index) #just constructs the object

#now running just the particle filter- this will probably not align with the data
filter$run(pars, save_history = TRUE)
h <- filter$history(1) #pull just one particle out
dim(h) # dimension of the index and time dimension at the places where you told it to output

plot(h["I", , ]) # showing a plot of the particle filter only

h <- filter$history() #now pull all particles

matplot(h["t", 1, ], t(h["cases", , ]), type = "l", col = "#00000011", 
        xlab = "Day", ylab = "Cases", las = 1)
points(cases ~ day, df_tidy, pch = 19, col = "red")  #compare to the data

matplot(h["t", 1, ], t(h["I", , ]), type = "l", col = "#00000011", 
        xlab = "Day", ylab = "Number of infecteds (I)", las = 1) #now plotting underlying infections

# now for the fitting
priors <- list(
  mcstate::pmcmc_parameter("gamma", 0.1, min = 0, prior = function(p)
    dgamma(p, shape = 1, scale = 0.2, log = TRUE)),
  mcstate::pmcmc_parameter("temp_scale", 1e-5, min=0))

transform <- function(theta) {
  as.list(theta)
} # this is where you may include a log transform or other

make_transform <- function(I0) {
  function(theta) {
    list(I0 =I0, #this runs everything with standardised I0 inital infecion
         gamma=theta[["gamma"]],
         temp_scale=theta[["temp_scale"]])
  }
}

pars <- list( gamma=0.2, temp_scale = 1e-5)
no_param <- length(pars)

vcv <- (1e-2) ^ 2 * diag(no_param) / no_param # this is for the proposal distribution
vcv
transform <- make_transform(10) #then this bounds I0 as 10

mcmc_pars <- mcstate::pmcmc_parameters$new(priors, vcv, transform) #getting proposal priors and starting parameters in the correct format

n_steps_in <- 1e3

control <- mcstate::pmcmc_control(
  n_steps = n_steps_in,
  progress = TRUE)
samples <- mcstate::pmcmc(mcmc_pars, filter, control = control)

#samples
plot(samples$probabilities[500:1000, "log_posterior"], type = "s",
     xlab = "Sample", ylab = "Log posterior") #this shows the chain

# plot result
pars <- list( gamma = samples$pars[n_steps_in,1], temp_scale=samples$pars[n_steps_in,2])
filter$run(pars, save_history = TRUE)
h <- filter$history()

matplot(h["t", 1, ], t(h["cases", , ]), type = "l", col = "#00000011", 
        xlab = "Day", ylab = "Cases", las = 1)
points(cases ~ day, df_tidy, pch = 19, col = "red")  
