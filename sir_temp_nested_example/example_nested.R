# a script for examining fitting a NHP in one region
# loading libraries
library(mcstate)
library(odin.dust)


#define the population size fo Sao Paulo
pop_dens <- 10 # based on density in Culot et al. Botacatu and others
area <- 248219 #based on Sao Paulo state

pop_size <- pop_dens * area #assume similar across the region *assumption* (we will assume this is equal in the 2 places)

# let's examine temperature
mean_Temp <- c(24.5, 20) # mean temperature in degrees for Sao Paulo and somewhere else

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
sir <- odin.dust::odin_dust("sir_temp_nested_example/sir_temp.R") # standard SIR model

# sourcing the compare function which defines what to fit to and how- here it is assuming a poisson likelihood
# sourcing the index function- this tells mcstate which model and data outputs to look for
source("sir_temp_nested_example/compare_sir_temp.R") # observation process could also appear here
source("sir_temp_nested_example/index_sir_temp.R")

n_particles <- 100
filter <- mcstate::particle_filter$new(data = data,
                                       model = sir,
                                       n_particles = n_particles,
                                       compare = compare,
                                       seed = 42L)

# now for the fitting
gamma <-   mcstate::pmcmc_parameter("gamma", 0.5, min = 0, prior = function(p)
    dgamma(p, shape = 1, scale = 0.2, log = TRUE))
temp_scale <-   mcstate::pmcmc_varied_parameter("temp_scale", c(1e-4, 1e-4), 
                                  populations = c("a", "b"),min=0) #this changes between populations

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

#proposal distributions
proposal_fixed <- (1e-2) ^ 2 * diag(1) / 2 # this is for the proposal distribution
colnames(proposal_fixed) <- c("gamma")
proposal_varied <- array(c(0.0001, 0.000095), c(1, 1, 2),
                         dimnames = list("temp_scale", "temp_scale", c("a", "b")))

transform <- make_transform(10) #then this bounds I0 as 10

#getting proposal priors and starting parameters in the correct format
mcmc_pars <- mcstate::pmcmc_parameters_nested$new(
  parameters = list(gamma = gamma, temp_scale=temp_scale),
  proposal_varied = proposal_varied,
  proposal_fixed = proposal_fixed,
  populations = c("a", "b")
) 

n_steps_in <- 1e3

control <- mcstate::pmcmc_control(
  n_steps = n_steps_in,
  progress = TRUE,
  save_trajectories = TRUE)

samples <- mcstate::pmcmc(mcmc_pars, filter, control = control)

#samples
par(mfrow = c(1, 2))
plot(samples$probabilities[1:1000, "log_posterior",1], type = "s",
     xlab = "Sample", ylab = "Log posterior") #this shows the chain
plot(samples$probabilities[1:1000, "log_posterior",2], type = "s",
     xlab = "Sample", ylab = "Log posterior") #this shows the chain

# plot result
t <- 0:100
observed_A=df_tidy[df_tidy$population == "A", ]
observed_B=df_tidy[df_tidy$population == "B", ]
output_A=t(samples$trajectories$state[6, 1, , ])
output_B=t(samples$trajectories$state[6, 2, , ])
par(mfrow = c(1, 2))
matplot(x=t, y=output_A, type = "l", lty = 1, col = "#00000005", xlab = "Day", ylab = "Daily incidence",
        main = "Daily Incidence - Pop. A",ylim=c(0,max(observed_A$cases,output_A)))
matplot(x=observed_A$day,observed_A$cases,type="p",pch=19,col="blue",add=TRUE)
matplot(x=t, y=output_B, type = "l", lty = 1, col = "#00000005", xlab = "Day", ylab = "Daily incidence",
        main = "Daily Incidence - Pop. B",ylim=c(0,max(observed_B$cases,output_B)))
matplot(x=observed_B$day,observed_B$cases,type="p",pch=19,col="blue",add=TRUE)
