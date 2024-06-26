#Extrapolating gamma and temperature coefficient values from nested sets of data generated by temperature dependent model
#NOT CURRENTLY WORKING DUE TO ISSUES WITH VARIABLE TEMPERATURE

library(mcstate)
library(odin.dust)

# SIR model - temperature dependent beta with gamma and temperature coefficient input as log values to aid fitting
sir <- odin.dust::odin_dust("time_varying_input/sir_var_temp.R")

# load example data generated by same model
data_example <- readRDS(file = "time_varying_input/example_data.Rds")

# load example case data
df_tidy <- data_example$example_data_nested

# Get parameters used to generate data
gen_params <- data_example$parameters_nested

# turn into data object
data <- particle_filter_data(df_tidy, time = "day", rate = 1, initial_time = 0, population="population")

# sourcing the compare function which defines what to fit to and how- here it is assuming a poisson likelihood
# sourcing the index function- this tells mcstate which model and data outputs to look for
source("sir_temp_example/compare_sir_temp.R") # observation process could also appear here
source("sir_temp_example/index_sir_temp.R")

#Create particle filter
filter <- particle_filter$new(data,  model = sir, n_particles=100, compare = compare, index = index)

#Set up prior probability distributions and initial, minimum and maximum values for parameters
priors <- list(pmcmc_parameter("log_gamma", initial = log(0.25), min = log(0.05), max = log(0.5), prior = function(p)
                               pnorm(p,mean=log(0.1),sd=1,log=TRUE)),
               pmcmc_parameter("log_temp_scale", initial = log(1.5e-5), min = log(1.0e-5), max = log(1.0e-3), prior = function(p)
                               pnorm(p,mean=log(1.0e-4),sd=1,log=TRUE)),
               pmcmc_varied_parameter("S0", initial = gen_params$S0, populations = c("a", "b")),
               pmcmc_varied_parameter("Temp", initial = gen_params$Temp[,1], populations = c("a", "b")))

#Transform function for other parameters
make_transform <- function(I0, Temp_0, Temp_m, dt, t_pts) {
  function(theta) {
    c(list(I0 = I0, Temp_0 = Temp_0, Temp_m = Temp_m, dt = dt, t_pts = t_pts), 
      as.list(theta))
  }
}

#Proposal distribution
proposal_fixed <- (1e-1) ^ 2 * diag(2)
colnames(proposal_fixed) <- c("log_gamma", "log_temp_scale")
proposal_varied <- array((1e-1) ^ 2 * diag(2), c(2, 2, 2),
                         dimnames = list(c("S0", "Temp"), 
                                         c("S0", "Temp"), 
                                         c("a", "b")))

#MCMC parameters
mcmc_pars <- pmcmc_parameters_nested$new(
  parameters = priors,
  proposal_varied = proposal_varied,
  proposal_fixed = proposal_fixed,
  populations = c("a", "b"),
  transform=make_transform(I0 = gen_params$I0,Temp_0 = gen_params$Temp_0, Temp_m = gen_params$Temp_m, dt = gen_params$dt, t_pts = gen_params$t_pts)
) 

#Set population and temperature for each region to be constant
mcmc_pars2 <- mcmc_pars$fix(fixed=cbind(a = c(S0=gen_params$S0[1],Temp=gen_params$Temp[1,]), 
                                        b = c(S0=gen_params$S0[2],Temp=gen_params$Temp[2,]))) 

#MCMC estimation
n_steps_in <- 1e3
control <- pmcmc_control(n_steps = n_steps_in, progress = TRUE, save_trajectories=TRUE)
samples <- pmcmc(mcmc_pars2, filter, control = control)

#Plot posterior likelihood progression over steps
matplot(x=c(1:nrow(samples$probabilities)),y=rowSums(samples$probabilities[,3,]),type="l",xlab="Iteration",ylab="Posterior likelihood (A+B)")
title("Posterior likelihood progression")

burnin=0.5*n_steps_in #Value from which to begin taking trajectory case values to plot

#Plot post-burn-in modelled values compared with observed data
par(mfrow=c(1,2))
matplot(x=c(1:101),y=t(samples$trajectories$state["cases",1,c(burnin:n_steps_in),]),type="p",pch=16,cex=0.25,col=1,
        ylim=c(0,max(c(samples$trajectories$state["cases",1,c(burnin:n_steps_in),],df_tidy$cases))),xlab="Day",ylab="Cases")
matplot(x=df_tidy$day[c(1:100)],y=df_tidy$cases[df_tidy$population=="A"],type="l",col=2,lwd=2.0,add=TRUE)
legend("topleft",c("Modelled", "Observed"),col=c(1,2),pch=c(16,NA),lty=c(0,1))
title("A")
matplot(x=c(1:101),y=t(samples$trajectories$state["cases",2,c(burnin:n_steps_in),]),type="p",pch=16,cex=0.25,col=1,
        ylim=c(0,max(c(samples$trajectories$state["cases",2,c(burnin:n_steps_in),],df_tidy$cases))),xlab="Day",ylab="Cases")
matplot(x=df_tidy$day[c(1:100)],y=df_tidy$cases[df_tidy$population=="B"],type="l",col=2,lwd=2.0,add=TRUE)
legend("topleft",c("Modelled", "Observed"),col=c(1,2),pch=c(16,NA),lty=c(0,1))
title("B")
par(mfrow=c(1,1))

parameters_out=list(log_gamma=samples$pars[,"log_gamma",1],log_temp_scale=samples$pars[,"log_temp_scale",])

#Plot convergence of log gamma values towards original input value
matplot(x=c(1:n_steps_in),y=parameters_out$log_gamma,type="l",col=1,xlab="Iteration",ylab="Value",
        ylim=range(c(parameters_out$log_gamma,gen_params$log_gamma)))
matplot(x=c(1,n_steps_in),y=rep(gen_params$log_gamma,2),type="l",lwd=2.0,lty=2,col=1,add=TRUE)
legend("topright",c("Log gamma (estimated)", "Log gamma (input)"),lty=c(1,2),col=c(1,1))

#Plot convergence of log temp_scale values towards original input values
matplot(x=c(1:n_steps_in),y=parameters_out$log_temp_scale[,1],type="l",col=1,xlab="Iteration",ylab="Value",
        ylim=range(c(parameters_out$log_temp_scale,gen_params$log_temp_scale)))
matplot(x=c(1:n_steps_in),y=parameters_out$log_temp_scale[,2],type="l",col=2,add=TRUE)
matplot(x=c(1,n_steps_in),y=rep(gen_params$log_temp_scale[1],2),type="l",lwd=2.0,lty=2,col=1,add=TRUE)
matplot(x=c(1,n_steps_in),y=rep(gen_params$log_temp_scale[2],2),type="l",lwd=2.0,lty=2,col=2,add=TRUE)
legend("right",c("Log temperature coefficient (estimated, A)", "Log temperature coefficient (input, A",
       "Log temperature coefficient (estimated, B)", "Log temperature coefficient (input, B)"),
       lty=c(1,2,1,2),col=c(1,1,2,2))
