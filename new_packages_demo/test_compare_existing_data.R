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
  # See https://mrc-ide.github.io/odin2/articles/functions.html#distribution-functions for available distributions
})

# ESTIMATING PARAMETER VALUES---------------------------------------------------
# Load example incidence data to use for estimation
data_all = read.csv(file = "new_packages_demo/cases_per_state_pnh_sp_mg_pr.csv", header=TRUE)
data=data.frame(time=data_all$day,cases=data_all$cases_A)

# Creating a list of input parameter values 
pars <- list(N = 1000, I0 = 10, beta = 0.2, gamma = 0.1)

# Setting up the model
sys <- dust2::dust_system_create(generator = sir, pars = pars, time = 0, dt = 1, deterministic = FALSE, 
                                 n_particles = 1, n_threads = 1, seed = 1, preserve_particle_dimension = TRUE)

# Set up time sequence for output data points
t <- seq(1, nrow(data), by = 1) 

# Initialize model to starting conditions
dust2::dust_system_set_state_initial(sys = sys)

# Run model over time sequence
y <- dust2::dust_system_simulate(sys = sys, times = t)
dim(y)

# Obtaining individual outputs
# Alternative to dust_unpack_state() as means of extracting labelled values
# Index of all outputs can be obtained; this may be useful if number of outputs is large and only certain values are needed
index = dust2::dust_unpack_index(obj = sys)
index
S = y[index$S,,]
I = y[index$I,,]
R = y[index$R,,]
incidence = y[index$incidence,,]

# Create filter
filter <- dust2::dust_filter_create(generator = sir, data = data, time_start = 0, n_particles = 20)

# Create packer - divide input parameters into estimated (beta + gamma) and fixed (N + I0)
# Significantly simplified from mcstate!
packer <- monty::monty_packer(scalar = c("beta","gamma"), array = NULL, fixed = list(N = 1000, I0 = 10))

# Set prior likelihood distributions for estimated parameters
# Here a simple uniform distribution is used for beta and gamma, with permitted minimum/maximum values
# See https://mrc-ide.github.io/odin2/articles/functions.html#distribution-functions for available distributions
prior <- monty::monty_dsl({
  beta ~ Uniform(min=0.05,max=0.5)
  gamma ~ Uniform(min=0.05,max=0.5)
})
prior$domain #Check limits (can adjust manually by adjusting values in prior$domain)

# Combine filter + packer
likelihood <- dust2::dust_likelihood_monty(obj = filter, packer = packer)

# Combine likelihood + prior to make posterior
posterior <- likelihood + prior

# Variance-covariance matrix for parameters to be estimated (beta, gamma)
# This is a key parameter to adjust when trying to improve estimation
vcv <- matrix(c(0.01, 0.005, 0.005, 0.01), 2, 2)
vcv

# Random walk sampler (other samplers are available)
sampler <- monty::monty_sampler_random_walk(vcv = vcv)

# Run samples to estimate parameters
# Initial values of estimated parameters set using parameter initial
# Initial values are another thing to vary when trying to improve estimation
n_chains=2
n_iterations=500
samples <- monty::monty_sample(model = posterior,sample = sampler,n_steps = n_iterations,
                               initial = array(rep(c(0.05, 0.25), n_chains), dim=c(2, n_chains)), n_chains = n_chains)

pars2 <- list(N = 1000, I0 = 10, beta = as.numeric(samples$pars[1, n_iterations, 1]),
              gamma = as.numeric(samples$pars[2, n_iterations, 1]))
sys2 <- dust2::dust_system_create(generator = sir, pars = pars2, time = 0, dt = 1, deterministic = FALSE,
                                  n_particles = 10, n_threads = 1, seed = 1, preserve_particle_dimension = TRUE)
index2 = dust2::dust_unpack_index(sys2)
dust2::dust_system_set_state_initial(sys2)
y2 <- dust2::dust_system_simulate(sys = sys2, times = data$time)

matplot(x = data$time, y = data$cases, type="l", col=1, xlab = "Month", ylab= "Incidence")
matplot(x = data$time, y = t(y2[index2$incidence,,]), type="p", pch = 1, col=2, add=TRUE)
legend("topright", c("Original data", "Data from estimated parameters"), col=c(1,2), lty=c(1,0), pch=c(NA,1))