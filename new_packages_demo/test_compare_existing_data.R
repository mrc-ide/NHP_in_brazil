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
  update(incidence) <- n_SI + 0.01 #Small value added to facilitate likelihood calculation
  
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

data_all = read.csv(file = "new_packages_demo/cases_per_state_pnh_sp_mg_pr.csv", header=TRUE)
plot(data_all$day,data_all$cases_A)
pts_select=c(480:1500)
#pts_select=c(1:nrow(data_all))
data=data.frame(time=c(1:length(pts_select)),cases=data_all$cases_A[pts_select])
plot(data$time,data$cases)
t <- seq(1, nrow(data), by = 1) 

pars <- list(N = 1000, I0 = 10, beta = 0.2, gamma = 0.1)
sys <- dust2::dust_system_create(generator = sir, pars = pars, time = 0, dt = 1, deterministic = FALSE, 
                                 n_particles = 100, n_threads = 1, seed = 1, preserve_particle_dimension = TRUE)
dust2::dust_system_set_state_initial(sys = sys)
y <- dust2::dust_system_simulate(sys = sys, times = t)

index = dust2::dust_unpack_index(obj = sys)
S = y[index$S,,]
I = y[index$I,,]
R = y[index$R,,]
incidence = y[index$incidence,,]

filter <- dust2::dust_filter_create(generator = sir, data = data, time_start = 0, n_particles = 100)
packer <- monty::monty_packer(scalar = c("beta","gamma"), array = NULL, fixed = list(N = 1000, I0 = 10))
prior <- monty::monty_dsl({
  beta ~ Uniform(min=0.01,max=1)
  gamma ~ Uniform(min=0.01,max=1)
})
likelihood <- dust2::dust_likelihood_monty(obj = filter, packer = packer, save_trajectories = TRUE)
posterior <- likelihood + prior
vcv <- matrix(c(0.005, 0.0025, 0.0025, 0.005), 2, 2)
sampler <- monty::monty_sampler_random_walk(vcv = vcv)

n_chains=1
n_iterations=1000
samples <- monty::monty_sample(model = posterior,sample = sampler,n_steps = n_iterations,
                               initial = array(rep(c(0.05, 0.25), n_chains), dim=c(2, n_chains)), n_chains = n_chains)

trajectory_final_incidence=samples$observations$trajectories[4,,n_iterations,]
matplot(x = data$time, y = data$cases, type = "l",col=1,xlab="Day",ylab="Incidence")
matplot(x=data$time,y=trajectory_final_incidence, type="p",col=2,pch=1,add=TRUE)
legend("topright", c("Original data", "Final trajectory"), col=c(1,2), lty=c(1,0), pch=c(NA,1))

pars2 <- list(N = pars$N, I0 = pars$I0, beta = as.numeric(samples$pars[1, n_iterations, 1]),
              gamma = as.numeric(samples$pars[2, n_iterations, 1]))
sys2 <- dust2::dust_system_create(generator = sir, pars = pars2, time = 0, dt = 1, deterministic = FALSE,
                                  n_particles = 100, n_threads = 1, seed = 1, preserve_particle_dimension = TRUE)
index2 = dust2::dust_unpack_index(sys2)
dust2::dust_system_set_state_initial(sys2)
y2 <- dust2::dust_system_simulate(sys = sys2, times = data$time)

matplot(x = data$time, y = data$cases, type="l", col=1, xlab = "Day", ylab= "Incidence")
matplot(x = data$time, y = t(y2[index2$incidence,,]), type="l", pch = 1, col=2, add=TRUE)
legend("topright", c("Original data", "Data from estimated parameters"), col=c(1,2), lty=c(1,0), pch=c(NA,1))