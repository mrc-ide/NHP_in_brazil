N <- S + I + R #total population size
p_SI <- 1 - exp(-(beta) * I / N) #probability of new infection
p_IR <- 1 - exp(-(gamma))        #probability of recovering 
n_IR <- rbinom(I, p_IR * dt)     #number of new recovered
n_SI <- rbinom(S, p_SI * dt)     #number of new infections

update(time) <- (step + 1) * dt 
update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
update(cases_cumul) <- cases_cumul + n_SI #cumulatively accounting for cases
update(cases_inc) <- n_SI
#update(cases_inc) <- if (step %% dt == 0) n_SI else cases_inc + n_SI #this pulls out the cases at the right points

initial(time) <- 0
initial(S) <- S0
initial(R) <- 0
initial(I) <- I0
initial(cases_cumul) <- 0
initial(cases_inc) <- 0

beta <- exp(log_temp_scale)*Temp[step+1]*(Temp[step+1]-Temp_0)*(Temp_m-Temp[step+1])^0.5 #Briere formulation of Mordecai et al 2017
gamma <- exp(log_gamma)

Temp[] <- user() #Input temperature at each time point
Temp_0 <- user(2.9285)  #1st parameter of Briere formulation
Temp_m <- user(40.1368) #2nd parameter of Briere formulation
log_temp_scale <- user() #logarithm of temperature scaling parameter, to be estimated 
log_gamma <- user() #logarithm of gamma, to be estimated
I0 <- user() #Starting number of infectious cases
S0 <- user(1000) #Starting susceptible population
t_pts <- user() #Total number of time points to be evaluated (length of temperature vector)
dt <- user() #Days per time point

dim(Temp) <- t_pts