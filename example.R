# a script for examining fitting a NHP in one region
library(mcstate)
library(odin.dust)


#define the population size fo Sao Paulo
pop_dens <- 10 # based on density in Culot et al. Botacatu and others
area <- 248219 #based on Sao Paulo state

pop_size <- pop_dens * area

# load data
df_tidy <- read.csv(system.file("sir_incidence.csv", package="mcstate")) # this is just example data from the package

# move n to start so mcstate knows where to look
data <- mcstate::particle_filter_data(df_tidy, 
                                      time = "day",
                                      rate = 1,
                                      initial_time = 0)

# sir model
sir <- odin.dust::odin_dust("sir.R") # standard SIR model

source("compare_sir.R") # this is a simple poisson likelihood with random noise added - observation process could appear here
source("index_sir.R")

#parms
pars <- list(beta = 0.25, gamma = 0.2) #leave others as defaults
mod <- sir$new(pars, 0, 100)
y <- mod$simulate(c(0, data$time_end))

filter <- mcstate::particle_filter$new(data, model = sir, n_particles=100,
                                       compare = compare, index = index) #just constructs the object


filter$run(pars, save_history = TRUE)
h <- filter$history(1)
dim(h) # dimension of the index and time dimension at the places where you told it to output

plot(h["I", , ])

h <- filter$history()

matplot(h["t", 1, ], t(h["cases", , ]), type = "l", col = "#00000011", 
        xlab = "Day", ylab = "Cases", las = 1)
points(cases ~ day, df_tidy, pch = 19, col = "red")  

matplot(h["t", 1, ], t(h["I", , ]), type = "l", col = "#00000011", 
        xlab = "Day", ylab = "Number of infecteds (I)", las = 1)

# now for the fitting
priors <- list(
  mcstate::pmcmc_parameter("beta", 0.2, min = 0),
  mcstate::pmcmc_parameter("gamma", 0.1, min = 0, prior = function(p)
    dgamma(p, shape = 1, scale = 0.2, log = TRUE)))

transform <- function(theta) {
  as.list(theta)
} # this is about as simple as this could be

make_transform <- function(I0) {
  function(theta) {
    list(I0 =I0, #this runs everything with standardised I0 inital infecion
         beta=theta[["beta"]],
         gamma=theta[["gamma"]])
  }
}
no_param <- length(pars)

vcv <- (1e-2) ^ 2 * diag(no_param) / no_param 
vcv
transform <- make_transform(10) #then this bounds I0 as 10

mcmc_pars <- mcstate::pmcmc_parameters$new(priors, vcv, transform)

n_steps_in <- 1e3

control <- mcstate::pmcmc_control(
  n_steps = n_steps_in,
  progress = TRUE)
samples <- mcstate::pmcmc(mcmc_pars, filter, control = control)

#samples
plot(samples$probabilities[, "log_posterior"], type = "s",
     xlab = "Sample", ylab = "Log posterior")


# plot result
pars <- list(beta = samples$pars[n_steps_in,1], gamma = samples$pars[n_steps_in,2])
filter$run(pars, save_history = TRUE)
h <- filter$history()

matplot(h["t", 1, ], t(h["cases", , ]), type = "l", col = "#00000011", 
        xlab = "Day", ylab = "Cases", las = 1)
points(cases ~ day, df_tidy, pch = 19, col = "red")  
