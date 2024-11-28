# This is a demonstration of key features of the updated modelling packages from Imperial College London:
# monty - A new package for parameter estimation, an update of mcstate - https://mrc-ide.github.io/monty/
# odin2 - Update of odin; domain specific language (DSL) for differential equation implementation - https://mrc-ide.github.io/odin2/
# dust2 - Update of dust; engine for running dynamic systems - https://mrc-ide.github.io/dust2/
# Note that eventually odin2 and dust2 will be folded into the existing odin and dust packages
# 
# For tutorials, see https://mrc-ide.github.io/odin-monty/odin.html

# INSTALLATION------------------------------------------------------------------
# Requires devtools package (run install.packages("devtools") if not installed)
devtools::install_github("mrc-ide/monty")
devtools::install_github("mrc-ide/dust2")
devtools::install_github("mrc-ide/odin2")
# Note the order - odin2 requires dust2 and monty to be installed

# MODEL CODE--------------------------------------------------------------------
# This code is a simple SIR model with incidence as an output
# Differences from implementation of the same model in odin are noted in comments
# Note that not much has changed in the implementation of this model!
# This code is reproduced in the file sir_odin2.R in this folder
# To load from the separate file, run sir <- odin2::odin("new_packages_demo/sir_odin2.R")
sir <- odin2::odin({
  N <- parameter(1000) # In odin2, parameter() is used for inputs where user() was used in odin
  I0 <- parameter(10)
  beta <- parameter(0.2)
  gamma <- parameter(0.1)
  
  p_SI <- 1 - exp(-beta * I / N * dt) # Note parameter dt - this is a time interval set when the model is set up via dust_system_create
  p_IR <- 1 - exp(-gamma * dt)
  n_SI <- Binomial(S, p_SI) # In odin2, Binomial() is used instead of rbinom()
  # See https://mrc-ide.github.io/odin2/articles/functions.html#distribution-functions for available distributions
  n_IR <- Binomial(I, p_IR)
  
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  update(incidence) <- n_SI
  
  initial(S) <- N - I0
  initial(I) <- I0
  initial(R) <- 0
  initial(incidence) <- 0
  
  # Adding goodness of fit to model (in odin2, this is implemented within the model code)
  # Only used when comparing to observed data
  cases <- data() #Data to compare
  cases ~ Poisson(incidence) #Distribution
})

# SETTING UP AND RUNNING MODEL--------------------------------------------------
# This differs significantly from the old package versions but reading output data in particular is easier

# Creating a list of input parameter values 
pars <- list(N = 1000, I0 = 10, beta = 0.2, gamma = 0.1)

# Setting up the model
# Note the parameter dt - this defines the time interval of each calculation step
# dt can be in any units but must have value 1 or less and be the inverse of an integer; here dt = 1 day
# For more flexibility (e.g. calculation steps 5 days apart), dt can have alternative units with a separate time interval defined
# Model can be set to deterministic mode by setting deterministic = TRUE
# Parameter seed can be used to produce reproducible stochastic output; set to NULL to randomize
# Parameter preserve_particle_dimension is set to TRUE so that output has same dimensionality even if n_particles = 1
sys <- dust2::dust_system_create(sir, pars, time = 0, dt = 1, deterministic = FALSE, n_particles = 1, n_threads = 1, seed = 1,
                                 preserve_particle_dimension = TRUE)

# Set up time sequence for output data points
# The output time sequence has the same units as the calculation time interval dt
# Note that the separation of the output time points needs to be compatible with dt
# For example, if dt = 1, the output time sequence cannot have fractional points
t <- seq(1, 100, by = 1) 

# Initialize model to starting conditions
dust2::dust_system_set_state_initial(sys)

# Run model over time sequence
# Note that if you want to re-run with a different time sequence, you must go back and re-run from line 63
# Output is an array with dimensions (number of output values, number of particles, number of time points)
y <- dust2::dust_system_simulate(sys, t)
dim(y)

# Convert model output to labelled values
results <- dust2::dust_unpack_state(sys, y)

#Plot SIR on graph
matplot(x = t, y = t(results$S), type="p", pch = 16, col=2, xlab = "Day", ylab = "", ylim=c(0, pars$N))
matplot(x = t, y = t(results$I), type="p", pch = 16, col=3, add=TRUE)
matplot(x = t, y = t(results$R), type="p", pch = 16, col=4, add=TRUE)
legend("topright",c("S","I","R"),pch=c(16,16,16),col=c(2,3,4))

#Plot incidence on graph
matplot(x = t, y = t(results$incidence), type="p", pch = 16, col=1, xlab = "Day", ylab= "Incidence")

# Obtaining individual outputs
# Index of all outputs can be obtained; this may be useful if number of outputs is large and only certain values are needed
index = dust2::dust_unpack_index(sys)
index
S = y[index$S,,]
I = y[index$I,,]
R = y[index$R,,]
incidence = y[index$incidence,,]

# ESTIMATING PARAMETER VALUES---------------------------------------------------
# Load example incidence data to use for estimation
data = read.csv(file = "new_packages_demo/example_data.csv", header=TRUE)

# Create filter
filter <- dust2::dust_filter_create(sir, data = data, time_start = 0, n_particles = 20)

# Create packer - divide input parameters into estimated (beta + gamma) and fixed (N + I0)
# Significantly simplified from mcstate!
packer <- monty::monty_packer(c("beta","gamma"),fixed=list(N = 1000, I0 = 10))

# Set prior likelihood distributions for estimated parameters
# Here a simple uniform distribution is used for beta and gamma, with permitted minimum/maximum values
# See https://mrc-ide.github.io/odin2/articles/functions.html#distribution-functions for available distributions
prior <- monty::monty_dsl({
  beta ~ Uniform(min=0.05,max=0.5)
  gamma ~ Uniform(min=0.05,max=0.5)
})
prior$domain #Check limits (can adjust manually by adjusting values in prior$domain)

# Combine filter + packer
likelihood <- dust2::dust_likelihood_monty(filter, packer)

# Combine likelihood + prior to make posterior
posterior <- likelihood + prior

# variance-covariance matrix
vcv <- matrix(c(0.01, 0.005, 0.005, 0.01), 2, 2)
vcv

# Random walk sampler (other samplers are available)
sampler <- monty::monty_sampler_random_walk(vcv)

# Run samples to estimate parameters
# Initial values of estimated parameters set using parameter initial
n_chains=1
n_iterations=500
samples <- monty::monty_sample(posterior,sampler,n_iterations,initial=array(rep(c(0.05,0.25),n_chains),
                                                                            dim=c(2,n_chains)),n_chains=n_chains)

# Plot posterior over iterations
matplot(x = c(1:n_iterations), y = samples$density[,1], type = "l", xlab = "Iteration", ylab = "Posterior")

# Plot parameter estimates over iterations (solid lines) with true values used to generate example data (dotted lines)
matplot(x = c(1:n_iterations), y = samples$pars[1,,1], type = "l", col = 2, xlab = "Iteration", ylab = "", ylim = c(0, 0.5))
matplot(x = c(1:n_iterations), y = samples$pars[2,,1], type = "l", col = 3, add = TRUE)
matplot(x = c(1, n_iterations), y = rep(0.2, 2), type = "l", col = 2, lty = 2, add = TRUE)
matplot(x = c(1, n_iterations), y = rep(0.1, 2), type = "l", col = 3, lty = 2, add = TRUE)
legend("topright", c("beta", "gamma"), lty=c(1,1), col=c(2,3))

# Calculate incidence using final estimated values and compare with example data
pars2 <- list(N = 1000, I0 = 10, beta = as.numeric(samples$pars[1,n_iterations,1]), 
                                gamma = as.numeric(samples$pars[2,n_iterations,1]))
sys2 <- dust2::dust_system_create(sir, pars2, time = 0, dt = 1, deterministic = FALSE, n_particles = 10, n_threads = 1, seed = 1,
                                 preserve_particle_dimension = TRUE)
index2 = dust2::dust_unpack_index(sys2)
dust2::dust_system_set_state_initial(sys2)
y2 <- dust2::dust_system_simulate(sys2, t = data$time)

matplot(x = data$time, y = data$cases, type="l", col=1, xlab = "Day", ylab= "Incidence")
matplot(x = data$time, y = t(y2[index2$incidence,,]), type="p", pch = 1, col=2, add=TRUE)
legend("topright",c("Original data","Data from estimated parameters"),col=c(1,2),lty=c(1,0),pch=c(NA,1))