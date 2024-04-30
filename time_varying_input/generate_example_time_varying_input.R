#Generate example nested data for fitting, with settings designed to give sparse data with many zero values
#Population and temperature vary between 2 regions
#Temperature varies with time
#gamma and temp_scale (logarithmic transform used to aid fitting) are constant across regions
#temp_scale is a value used to calculate beta (multiplier of exponential function of temperature)

# SIR model
sir <- odin.dust::odin_dust("time_varying_input/sir_var_temp.R") # standard SIR model

# Set parameters for all regions
t_pts=200 #Number of time points to run
pop_dens <- 0.01  #Population density
area <- c(220000, 248219) 
pop_size <- pop_dens * area 
I0 <- 10 #Initial number of infectious cases in each area
log_gamma = log(0.1) #Gamma parameter (defined logarithmically to aid fitting)
log_temp_scale=log(1.2e-4) #Temperature scale parameters (defined logarithmically to aid fitting) used to calculate beta parameter
dt <- 1.0 #Time scale (days)

#Set temperature values
Temp <- array(NA,dim=c(2,t_pts))
Temp[1,] <- 24.5*(0.25*sin(((c(1:t_pts)-(t_pts/2))*(2*pi))/t_pts) + 0.75)
Temp[2,] <- 20.0*(0.25*sin((c(1:t_pts)*(2*pi))/t_pts) + 0.75)
matplot(t(Temp),type="l",lty=1,xlab="Date",ylab="Temperature")
legend("topright",legend=c("A","B"),lty=c(1,1),col=c(1,2))

#Check R0 values
Temp_0 <- 2.9285 
Temp_m <- 40.1368 
beta <- exp(log_temp_scale)*Temp*(Temp-Temp_0)*(Temp_m-Temp)^0.5
R0=beta/exp(log_gamma)
matplot(t(R0),type="l",lty=1,xlab="Date",ylab="R0")
legend("topright",legend=c("A","B"),lty=c(1,1),col=c(1,2))

#Set up nested parameters
n_param_sets=2
assertthat::assert_that(n_param_sets==length(pop_size)) #Check input number of parameter sets matches number of population values
assertthat::assert_that(n_param_sets==dim(Temp)[1]) #Check input number of parameter sets matches number of sets of temperature values
pars=list()
for(i in 1:n_param_sets){
  pars[[i]]=list(log_temp_scale = log_temp_scale, log_gamma = log_gamma,I0 = I0, S0 = pop_size[i],
                 Temp = Temp[i,], Temp_0 = Temp_0, Temp_m = Temp_m, t_pts = t_pts, dt = dt)
}

x <- sir$new(pars=pars, time = 0, n_particles = 1, n_threads = 1, deterministic = FALSE, pars_multi = TRUE)
x_res <- array(NA, dim = c(6, n_param_sets, t_pts)) #Outputs from SIR model - time, S, I, R, cumulative cases, case incidence
for(step in 1:t_pts){
  x_res[,,step] <- x$run(step)
}

output_frame=data.frame(cases=rep(NA,2*t_pts),day=rep(c(1:t_pts),2),population=as.factor(c(rep("A",t_pts),rep("B",t_pts))))
output_frame$cases[output_frame$population=="A"]=x_res[6,1,]
output_frame$cases[output_frame$population=="B"]=x_res[6,2,]

matplot(x=output_frame$day[c(1:t_pts)],y=output_frame$cases[c(1:t_pts)],type="b",col=1,pch=19,
        xlab="Day",ylab="Cases",ylim=c(0,max(output_frame$cases)))
matplot(x=output_frame$day[c(1:t_pts)+t_pts],y=output_frame$cases[c(1:t_pts)+t_pts],type="b",col=2,pch=19,add=TRUE)

#Save inputs and outputs - one set for region A to use with example_single.R, one nested set to use with example_nested.R
#Population and temperature needed for fitting, log_gamma and log_temp_scale used to check success of fit
#NB - PAST FITTING MAY HAVE BEEN AFFECTED BY DISCREPANCIES IN I0
outputs <- list(example_data_nested = output_frame, 
                parameters_nested = list(log_temp_scale = log_temp_scale, log_gamma = log_gamma, S0 = pop_size, I0 = I0,
                                         Temp = Temp, Temp_0 = Temp_0, Temp_m = Temp_m, dt = dt, t_pts = t_pts),
                example_data_single=output_frame[output_frame$population == "A", c(1:2)],
                parameters_single = list(log_temp_scale = log_temp_scale, log_gamma = log_gamma, S0 = pop_size[1], I0 = I0,
                                         Temp = Temp[1,], Temp_0 = Temp_0, Temp_m = Temp_m, dt = dt, t_pts = t_pts))
#saveRDS(outputs,file = "time_varying_input/example_data.Rds")


