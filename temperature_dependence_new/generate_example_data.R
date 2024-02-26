#Generate example nested data for fitting
#Population and temperature vary between regions, gamma and temp_scale (multiplier used to calculate beta) constant

# Old example data as template
df_tidy <- read.csv(system.file("nested_sir_incidence.csv", package="mcstate"), stringsAsFactors = TRUE)
t_pts=max(df_tidy$day)

# SIR model
sir <- odin.dust::odin_dust("temperature_dependence_new/sir_temp_new.R") # standard SIR model

# Set parameters for all regions
pop_dens <- 0.1 
area <- c(220000, 248219) 
pop_size <- pop_dens * area 
log_gamma = log(0.1)
log_temp_scale=log(c(1.0e-4, 1.5e-4))

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
  pars[[i]]=list(log_temp_scale = log_temp_scale[i], log_gamma = log_gamma,I0 = 10, S0 = pop_size[i],
                 Temp = Temp[i], Temp_0 = Temp_0, Temp_m = Temp_m)
}

x <- sir$new(pars=pars, time = 0, n_particles = 1, n_threads = 1, deterministic = FALSE, pars_multi = TRUE)

x_res <- array(NA, dim = c(6, n_param_sets, t_pts)) #6 outputs from SIR model - time, S, I, R, cumulative cases, case incidence
for(step in 1:t_pts){
  x_res[,,step] <- x$run(step)
}

df_tidy_new=df_tidy
df_tidy_new$cases[df_tidy_new$population=="A"]=x_res[6,1,]
df_tidy_new$cases[df_tidy_new$population=="B"]=x_res[6,2,]

matplot(x=df_tidy_new$day[c(1:t_pts)],y=df_tidy_new$cases[c(1:t_pts)],type="b",col=1,pch=19,
        xlab="Day",ylab="Cases",ylim=c(0,max(df_tidy_new$cases)))
matplot(x=df_tidy_new$day[c(1:t_pts)+t_pts],y=df_tidy_new$cases[c(1:t_pts)+t_pts],type="b",col=2,pch=19,add=TRUE)

outputs <- list(example_data_nested = df_tidy, 
                parameters_nested = list(log_temp_scale = log_temp_scale, log_gamma = log_gamma, I0 = 10, S0 = pop_size,
                                                                        Temp = Temp, Temp_0 = Temp_0, Temp_m = Temp_m),
                example_data_single=df_tidy_new[df_tidy_new$population == "A", c(1:2)],
                parameters_single = list(log_temp_scale = log_temp_scale[1], log_gamma = log_gamma, I0 = 10, S0 = pop_size[1],
                                         Temp = Temp[1], Temp_0 = Temp_0, Temp_m = Temp_m))
saveRDS(outputs,file = "temperature_dependence_new/example_data.Rds")


