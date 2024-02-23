#Generate example nested data for fitting
#Population and temperature vary between regions, gamma and temp_scale (multiplier used to calculate beta) constant

# Old example data as template
df_tidy <- read.csv(system.file("nested_sir_incidence.csv", package="mcstate"), stringsAsFactors = TRUE)
t_pts=max(df_tidy$day)

# SIR model
sir <- odin.dust::odin_dust("sir_temp_nested_example/sir_temp.R") # standard SIR model

# Set parameters for all regions
pop_dens <- 0.01 
area <- c(220000, 248219) 
pop_size <- pop_dens * area 
gamma = 0.2
temp_scale=c(4.0e-4, 5.0e-4)

#Check R0 values
Temp <- c(24.5, 20)
Temp_0 <- 2.9285 
Temp_m <- 40.1368 
beta <- temp_scale*Temp*(Temp-Temp_0)*(Temp_m-Temp)^0.5
R0=beta/gamma

#Set up nested parameters
n_param_sets=2
pars=list()
for(i in 1:n_param_sets){
  pars[[i]]=list(temp_scale = temp_scale[i], gamma = gamma,I0 = 10, S0 = pop_size[i], Temp = Temp[i], Temp_0 = Temp_0, Temp_m = Temp_m)
}

x <- sir$new(pars=pars, time = 0, n_particles = 1, n_threads = 1, deterministic = FALSE, pars_multi = TRUE)

x_res <- array(NA, dim = c(6, n_param_sets, t_pts)) #6 outputs from SIR model - time, S, I, R, cumulative cases, case incidence
for(step in 1:t_pts){
  x_res[,,step] <- x$run(step)
}

df_tidy_new=df_tidy
df_tidy_new$cases[df_tidy_new$population=="A"]=x_res[6,1,]
df_tidy_new$cases[df_tidy_new$population=="B"]=x_res[6,2,]

plot(cases ~ day, df_tidy_new[df_tidy_new$population == "B", ],
     type = "o", xlab = "Day", ylab = "New cases", pch = 19)
lines(cases ~ day, df_tidy_new[df_tidy_new$population == "A", ],
      type = "o", xlab = "Day", ylab = "New cases", pch = 19, col = 2)
legend("topright", col = 2:1, legend = c("A", "B"), lwd = 1)

#write.csv(df_tidy_new, file = "sir_temp_nested_example/nested_sir_incidence.csv", row.names=FALSE)

