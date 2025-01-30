# Non-human primate SIR model, written in format for compilation by odin.dust package
# This code can be modified to change the NHP model (e.g. to have it accept environmental variables or
# have different default values of user-defined quantities)

# User-defined quantities
R0[] <- user()              #Basic reproduction number as a function of time
gamma <- user(0.1)          #Rate of recovery (inverse of average infectious period)
I0 <- user(10)              #Initial number of infectious individuals
S0 <- user(1000)            #Initial number of susceptible individuals
dt <- user(1)               #Time increment in days
mortality <- user(0.1/365)  #Mortality rate per day used for demographic turnover (inverse of average lifespan in days) 
n_pts <- user()             #Number of time points (dimension of R0 vector)

N <- S + I + R                   #Total population size
beta <- R0[step]*gamma           #Infection rate
p_SI <- 1 - exp(-(beta) * I / N) #Probability of new infection
p_IR <- 1 - exp(-(gamma))        #Probability of recovering 
n_IR <- rbinom(I, p_IR * dt)     #Number of new recovered
n_SI <- rbinom(S, p_SI * dt)     #Number of new infections

update(time) <- (step + 1) * dt 
update(S) <- S - n_SI + ((N-S)*mortality*dt)
update(I) <- I + n_SI - n_IR - (I*mortality*dt)
update(R) <- R + n_IR - (R*mortality*dt)

initial(time) <- 0
initial(S) <- S0
initial(R) <- 0
initial(I) <- I0

dim(R0) <- n_pts