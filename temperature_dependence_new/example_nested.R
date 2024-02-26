#New example using data created using temperature dependent model

# a script for examining fitting a NHP in one region
# loading libraries
library(mcstate)
library(odin.dust)


#define the population size for Sao Paolo
pop_dens <- 0.001 # adjusted based on this example dataset
area <- c(220000, 248219) #based on Sao Paolo state

pop_size_in <- pop_dens * area #assume similar across the region *assumption* (we will assume this is equal in the 2 places)

# let's examine temperature
mean_Temp <- c(24.5, 20) # mean temperature in degrees for Sao Paulo and somewhere else

# load example data
df_tidy <- read.csv("sir_temp_nested_example/nested_sir_incidence.csv",stringsAsFactors = TRUE) #New example data generated using sir_temp

# turn into data object
data <- mcstate::particle_filter_data(df_tidy, 
                                      time = "day",
                                      rate = 1,
                                      initial_time = 0,
                                      population="population")

par(mfrow = c(1, 1))
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
gamma <-   mcstate::pmcmc_parameter("gamma", 0.1, min = 0) #Does not vary between populations
log_temp_scale <-   mcstate::pmcmc_varied_parameter("log_temp_scale", c(log(1e-4), log(1.5e-4)), populations = c("a", "b")) #this changes between populations
Temp <-   mcstate::pmcmc_varied_parameter("Temp", mean_Temp, populations = c("a", "b"),min=0) #this changes between populations
S0 <- mcstate::pmcmc_varied_parameter("S0", pop_size_in, populations = c("a", "b")) #this changes between populations

#proposal distributions
proposal_fixed <- (1e-2) ^ 2 * diag(1) / 2 # this is for the proposal distribution
colnames(proposal_fixed) <- c("gamma")
proposal_varied <- array((1e-2) ^ 2 * diag(3) / 2, c(3, 3, 2),
                         dimnames = list(c("log_temp_scale", "S0", "Temp"), 
                                         c("log_temp_scale", "S0", "Temp"), 
                                         c("a", "b")))

#getting proposal priors and starting parameters in the correct format
mcmc_pars <- mcstate::pmcmc_parameters_nested$new(
  parameters = list(gamma = gamma, log_temp_scale=log_temp_scale, S0 = S0, Temp = Temp),
  proposal_varied = proposal_varied,
  proposal_fixed = proposal_fixed,
  populations = c("a", "b")
) 

# this is where you can tell mcstate not to estimate these parameters but
# keep them varied between patches (you can comment this out if you do 
# want to estimate population size)
mcmc_pars2 <- mcmc_pars$fix(fixed=cbind(a = c(S0=pop_size_in[1], Temp=mean_Temp[1]), 
                                        b = c(S0=pop_size_in[2], Temp=mean_Temp[2]))) 

n_steps_in <- 1e3

control <- mcstate::pmcmc_control(
  n_steps = n_steps_in,
  progress = TRUE,
  save_trajectories = TRUE,
  rerun_random = FALSE,
  rerun_every=10)

samples <- mcstate::pmcmc(mcmc_pars, filter, control = control)

burnin1=1
#samples
par(mfrow = c(1, 2))
plot(samples$probabilities[1:n_steps_in, "log_posterior",1], type = "s",
     xlab = "Sample", ylab = "Log posterior") #this shows the chain
plot(samples$probabilities[1:n_steps_in, "log_posterior",2], type = "s",
     xlab = "Sample", ylab = "Log posterior") #this shows the chain

# temp scaling
par(mfrow = c(1, 2))
plot(exp(samples$pars[c(burnin1:n_steps_in), "log_temp_scale",1]), type = "s",
     xlab = "Sample", ylab = "Temp scale") #this shows the chain
plot(exp(samples$pars[c(burnin1:n_steps_in), "log_temp_scale",2]), type = "s",
     xlab = "Sample", ylab = "Temp scale") #this shows the chain

# plot result
burnin2=500
t <- 0:100
observed_A=df_tidy[df_tidy$population == "A", ]
observed_B=df_tidy[df_tidy$population == "B", ]
output_A=t(samples$trajectories$state[6, 1, c(burnin2:n_steps_in), ])
output_B=t(samples$trajectories$state[6, 2, c(burnin2:n_steps_in), ])
par(mfrow = c(1, 2))
matplot(x=t, y=rowMeans(output_A), type = "l", lty = 1, col = 1, xlab = "Day", ylab = "Daily incidence",
        main = "Daily Incidence - Pop. A",ylim=c(0,max(df_tidy$cases)))
matplot(x=observed_A$day,observed_A$cases,type="p",pch=19,col="blue",add=TRUE)
matplot(x=t, y=rowMeans(output_B), type = "l", lty = 1, col = 1, xlab = "Day", ylab = "Daily incidence",
        main = "Daily Incidence - Pop. B",ylim=c(0,max(df_tidy$cases)))
matplot(x=observed_B$day,observed_B$cases,type="p",pch=19,col="blue",add=TRUE)
par(mfrow = c(1, 1))