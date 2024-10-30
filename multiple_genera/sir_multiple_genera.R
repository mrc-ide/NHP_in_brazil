# SEIR (Susceptible, Exposed, Infectious, Recovered) yellow fever model for non-human primates (NHPs),
# with multiple genera

# TO DO (to add or to include in future version):
# -Incorporate deaths due to YF (infection fatality rate parameter)
# -Update population modelling to allow for population extinction or recovery (look at existing ecological models)
# -Set number of births/deaths at each time point using binomial expression like new infections, with
# -Incorporate introductions of individuals from outside
# -Make force of infection genus dependent
# -Add vector model including functionality for multiple vector species
# -Make force of infection time-varying (by making beta and/or FOI_in time-varying)
# -Make force of infection (via beta and/or FOI_in) driven by temperature/precipitation input

dt <- user() #Time increment in days
initial(time) <- 0 #Initial value of time in days
update(time) <- time + dt

#Parameters---------------------------------------------------------------------
beta <- user() #beta parameter for transmission [TODO: MAKE TIME-VARYING/CLIMATE-DRIVEN]
t_incubation <- user() #Length in days of yellow fever incubation period in mosquito vectors [TODO: MAKE VAR BY SPECIES]
t_latent[] <- user() #Length in days of latent period in NHPs exposed to yellow fever
t_infectious[] <- user() #Length of infectious period in NHPs with yellow fever
FOI_in <- user() #Force of infection (per day) for outside introduction [TODO: MAKE TIME-VARYING/CLIMATE-DRIVEN]
R_dP[] <- user() #Rate of population turnover by day (birth rate/death rate, assumed equal for stable population)
n_gen <- user() #Number of genera to look at

#Initial conditions-------------------------------------------------------------
S_0[] <- user() #Susceptible NHP population by genus at start
E_0[] <- user() #Exposed NHP population by genus at start
I_0[] <- user() #Infectious NHP population by genus at start
R_0[] <- user() #Recovered NHP population by genus at start

Pmin <- 1.0e-99 #Minimum population setting to avoid negative numbers
FOI_max <- 1.0 #Upper threshold for total force of infection to avoid more infections than individuals
rate1[1:n_gen] <- dt/(t_incubation+t_latent[i]) #Rate of transference E->I
rate2[1:n_gen] <- dt/t_infectious[i] #Rate of transference I->R

FOI_sum <-  min(FOI_max,dt*((beta*sum(I)*sum(inv_P)) + FOI_in)) #Total force of infection
E_new[1:n_gen] <- rbinom(as.integer(S[i]), FOI_sum) #New exposed NHPs by genus
I_new[1:n_gen] <- E[i]*rate1[i]     #New infectious NHPs by genus
R_new[1:n_gen] <- I[i]*rate2[i]     #New recovered NHPs by genus
P[1:n_gen] <- S[i] + R[i] #Total NHP population by genus (excluding E+I)
inv_P[1:n_gen] <- 1.0/P[i]

#Updates to output values at each time increment--------------------------------
update(FOI_total) <- FOI_sum
update(S[1:n_gen]) <- max(Pmin,S[i] - E_new[i] + R_dP[i] - (R_dP[i]*S[i]*inv_P[i]))
update(E[1:n_gen]) <- max(Pmin,E[i] + E_new[i] - I_new[i])
update(I[1:n_gen]) <- max(Pmin,I[i] + I_new[i] - R_new[i])
update(R[1:n_gen]) <- max(Pmin,R[i] + R_new[i] - (R_dP[i]*R[i]*inv_P[i]))
update(C[1:n_gen]) <- I_new[i]

#Initial values-----------------------------------------------------------------
initial(FOI_total) <- FOI_in
initial(S[1:n_gen]) <- S_0[i]
initial(E[1:n_gen]) <- E_0[i]
initial(I[1:n_gen]) <- I_0[i]
initial(R[1:n_gen]) <- R_0[i]
initial(C[1:n_gen]) <- 0

#Dimensions---------------------------------------------------------------------
dim(S) <- n_gen
dim(E) <- n_gen
dim(I) <- n_gen
dim(R) <- n_gen
dim(C) <- n_gen
dim(E_new) <- n_gen
dim(I_new) <- n_gen
dim(R_new) <- n_gen

dim(P) <- n_gen
dim(inv_P) <- n_gen

dim(S_0) <- n_gen
dim(E_0) <- n_gen
dim(I_0) <- n_gen
dim(R_0) <- n_gen

dim(t_latent) <- n_gen
dim(t_infectious) <- n_gen
dim(R_dP) <- n_gen
dim(rate1) <- n_gen
dim(rate2) <- n_gen

