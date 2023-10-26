

R0[] <- user()
gamma <- user(0.1)
I0 <- user(10)
S0 <- user(1000)
dt <- user(1)
mortality <- user(0.1/365)
#freq <- 1.0/dt
n_pts <- user()

N <- S + I + R #total population size
beta <- R0[step]*gamma
p_SI <- 1 - exp(-(beta) * I / N) #probability of new infection
p_IR <- 1 - exp(-(gamma))        #probability of recovering 
n_IR <- rbinom(I, p_IR * dt)     #number of new recovered
n_SI <- rbinom(S, p_SI * dt)     #number of new infections

update(time) <- (step + 1) * dt 
update(S) <- S - n_SI + ((N-S)*mortality*dt)
update(I) <- I + n_SI - n_IR - (I*mortality*dt)
update(R) <- R + n_IR - (R*mortality*dt)
#update(cases_cumul) <- cases_cumul + n_SI #cumulatively accounting for cases
#update(cases_inc) <- if (step %% freq == 0) n_SI else cases_inc + n_SI #this pulls out the cases at the right points

initial(time) <- 0
initial(S) <- S0
initial(R) <- 0
initial(I) <- I0
#initial(cases_cumul) <- 0
#initial(cases_inc) <- 0

dim(R0) <- n_pts