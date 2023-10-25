update(S_tot) <- S_tot - sum(n_SI) # these are total populations
update(I_tot) <- I_tot + sum(n_SI) - sum(n_IR)
update(R_tot) <- R_tot + sum(n_IR)


N <- S_tot + I_tot + R_tot #total population size
p_SI[] <- 1 - exp(-(beta[i]) * I[i] / N) #probability of new infection
p_IR <- 1 - exp(-(gamma))        #probability of recovering 

n_SI[] <- rbinom(S[i], p_SI[i] * dt)     #number of new infections
n_IR[] <- rbinom(I[i], p_IR * dt)     #number of new recovered

update(time) <- (step + 1) * dt 

update(S[]) <- S[i] - n_SI[i]
update(I[]) <- I[i] + n_SI[i] - n_IR[i]
update(R[]) <- R[i] + n_IR[i]

update(cases_cumul[]) <- cases_cumul[i] + n_SI[i] #cumulatively accounting for cases
update(cases_inc[]) <- if (step %% freq == 0) n_SI[i] else cases_inc[i] + n_SI[i] #this pulls out the cases at the right points


## Initial states:
initial(S_tot) <- sum(S0)
initial(I_tot) <- sum(I0)
initial(R_tot) <- 0

initial(time) <- 0
initial(S[]) <- S0[i]
initial(R[]) <- 0
initial(I[]) <- I0[i]
initial(cases_cumul[]) <- 0
initial(cases_inc[]) <- 0

beta[] <- temp_scale*Temp[i]*(Temp[i]-Temp_0)*(Temp_m-Temp[i])^0.5 #Briere formulation of Mordecai et al 2017

Temp[] <- user() #as input
Temp_0 <- user(2.9285) #to be estimated (but would need more data - this value is from Gaythorpe et al. 2020)
Temp_m <- user(40.1368) #to be estimated (but would need more data - this value is from Gaythorpe et al. 2020)
temp_scale <- user(3e-4) #to be estimated 
gamma <- user(0.1) #to be estimated
I0[] <- user() #as input
S0[] <- user() #as input 

freq <- user(1)
dt <- 1.0 / freq

# dimensions
N_patch <- user(2)
dim(S0) <- N_patch
dim(I0) <- N_patch
dim(S) <-  N_patch
dim(I) <-  N_patch
dim(R) <-  N_patch
dim(cases_inc) <- N_patch
dim(n_SI) <- N_patch
dim(n_IR) <- N_patch
dim(p_SI) <- N_patch
dim(Temp) <- N_patch
dim(cases_cumul) <- N_patch
dim(beta) <- N_patch