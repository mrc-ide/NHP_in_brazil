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
update(cases_inc) <- if (step %% freq == 0) n_SI else cases_inc + n_SI #this pulls out the cases at the right points

initial(time) <- 0
initial(S) <- S0
initial(R) <- 0
initial(I) <- I0
initial(cases_cumul) <- 0
initial(cases_inc) <- 0

beta <- -temp_scale*Temp*(Temp–Temp_0)*(Temp_m–Temp)^0.5 #Briere formulation of Mordecai et al 2017

Temp <- user(25) #as input
Temp_0 <- user(2.9285) #to be estimated (but would need more data - this value is from Gaythorpe et al. 2020)
Temp_m <- user(40.1368) #to be estimated (but would need more data - this value is from Gaythorpe et al. 2020)
temp_scale <- user(3e-4) #to be estimated (but would need more data - this value is from Gaythorpe et al. 2020)

gamma <- user(0.1) #to be estimated
I0 <- user(10) #as input
S0 <- user(1000) #as input 

freq <- user(1)
dt <- 1.0 / freq