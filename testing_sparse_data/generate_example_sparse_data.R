#Generate example nested data for fitting, with settings designed to give sparse data with many zero values
#Population and temperature vary between 2 regions
#gamma (calculated from log(gamma) to aid fitting) is constant across regions
#temp_scale is a value used to calculate beta (multiplier of exponential function of temperature)
#temp_scale (calculated from log(temp_scale) to aid fitting) varies between regions

#df_tidy <- read.csv(system.file("nested_sir_incidence.csv", package="mcstate"), stringsAsFactors = TRUE)

# SIR model
sir <- odin.dust::odin_dust("temperature_dependence_new/sir_temp_new.R") # standard SIR model

# Set parameters for all regions
t_pts=200 #Number of time points to run
pop_dens <- 0.001  #Population density
area <- c(220000, 248219) 
pop_size <- pop_dens * area 
log_gamma = log(0.1) #Gamma parameter (defined logarithmically to aid fitting)
log_temp_scale=log(c(1.2e-4, 9.5e-5)) #Temperature scale parameters (defined logarithmically to aid fitting) used to calculate beta parameter

#Check R0 values
Temp <- c(24.5, 20)
Temp_0 <- 2.9285 
Temp_m <- 40.1368 
beta <- exp(log_temp_scale)*Temp*(Temp-Temp_0)*(Temp_m-Temp)^0.5
R0=beta/exp(log_gamma)

#Set up nested parameters
n_param_sets=2
pars=list()
for(i in 1:n_param_sets){
  pars[[i]]=list(log_temp_scale = log_temp_scale[i], log_gamma = log_gamma,I0 = 2, S0 = pop_size[i],
                 Temp = Temp[i], Temp_0 = Temp_0, Temp_m = Temp_m)
}

x <- sir$new(pars=pars, time = 0, n_particles = 1, n_threads = 1, deterministic = FALSE, pars_multi = TRUE)

x_res <- array(NA, dim = c(6, n_param_sets, t_pts)) #6 outputs from SIR model - time, S, I, R, cumulative cases, case incidence
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
outputs <- list(example_data_nested = output_frame, 
                parameters_nested = list(log_temp_scale = log_temp_scale, log_gamma = log_gamma, I0 = 10, S0 = pop_size,
                                                                        Temp = Temp, Temp_0 = Temp_0, Temp_m = Temp_m),
                example_data_single=output_frame[output_frame$population == "A", c(1:2)],
                parameters_single = list(log_temp_scale = log_temp_scale[1], log_gamma = log_gamma, I0 = 10, S0 = pop_size[1],
                                         Temp = Temp[1], Temp_0 = Temp_0, Temp_m = Temp_m))
#saveRDS(outputs,file = "testing_sparse_data/example_data.Rds")


