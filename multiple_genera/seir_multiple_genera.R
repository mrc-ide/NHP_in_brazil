# SEIR (Susceptible, Exposed, Infectious, Recovered) yellow fever model for non-human primates (NHPs),
# with multiple genera

# TO DO (to add or to include in future version):
# -Update population modelling to allow for population extinction or recovery (look at existing ecological models, ongoing)
# -Incorporate introductions of individuals from outside (and outward migration)
# -Make force of infection genus dependent
# -Add vector model including functionality for multiple vector species
# -Make force of infection (via beta and/or FOI_in) driven by temperature/precipitation input
# -Change model of force of infection, latency, etc. for better short-term timescale modelling

dt <- user() #Time increment in days
initial(time) <- 0 #Initial value of time in days
update(time) <- time + dt

#Parameters---------------------------------------------------------------------
n_t_pts <- user() #Number of time points to be run
temp[] <- user() #Time-varying temperature used to calculate beta and FOI_in
precip[] <- user() #Time-varying precipitation used to calculate beta and FOI_in
c_beta_temp <- user() #Coefficient used to calculate temperature component of beta parameter
c_beta_precip <- user() #Coefficient used to calculate precipitation component of beta parameter
c_FOI_in_temp <- user() #Coefficient used to calculate temperature component of FOI_in parameter
c_FOI_in_precip <- user() #Coefficient used to calculate precipitation component of FOI_in parameter
t_incubation <- user() #Length in days of yellow fever incubation period in mosquito vectors [TODO: MAKE VAR BY SPECIES]
t_latent[] <- user() #Length in days of latent period in NHPs exposed to yellow fever
t_infectious[] <- user() #Length of infectious period in NHPs with yellow fever
ifr[] <- user() #Infection fatality rate
FOI_in[] <- user() 
mu[] <- user() #Rate of population turnover by day (birth rate/death rate, assumed equal for stable population)
n_gen <- user() #Number of genera to look at

#Initial conditions-------------------------------------------------------------
S_0[] <- user() #Susceptible NHP population by genus at start
E_0[] <- user() #Exposed NHP population by genus at start
I_0[] <- user() #Infectious NHP population by genus at start
R_0[] <- user() #Recovered NHP population by genus at start

Pmin <- 0.0 #Minimum population setting to avoid negative numbers
FOI_max <- 1.0 #Upper threshold for total force of infection to avoid more infections than individuals
rate1[1:n_gen] <- dt/(t_incubation+t_latent[i]) #Rate of transference E->I
rate2[1:n_gen] <- dt/t_infectious[i] #Rate of transference I->R

beta <- 0 #TBA #Time-varying beta parameter for transmission
FOI_in <- 0 #TBA #Time-varying force of infection (per day) for outside introduction
FOI_sum <-  min(FOI_max,dt*(((beta[step+1]*sum(I))/sum(P)) + FOI_in[step+1])) #Total force of infection
births[1:n_gen] <- rbinom(as.integer(P[i]), mu[i]) #New births by genus
deaths_nat[1:n_gen] <- rbinom(as.integer(P[i]), mu[i]) #New deaths by genus (excluding YF deaths)
E_new[1:n_gen] <- rbinom(as.integer(S[i]), FOI_sum) #New exposed NHPs by genus
I_new[1:n_gen] <- E[i]*rate1[i]     #New infectious NHPs by genus
R_new0[1:n_gen] <- I[i]*rate2[i]     #Number of NHPS removed from I by genus (recovery from or death due to YF)
R_new[1:n_gen] <- R_new0[i]*(1.0-ifr[i]) #New recovered NHPs by genus
deaths_YF[1:n_gen] <- R_new0[i]*ifr[i]  #New deaths of NHPs by genus
P[1:n_gen] <- S[i] + E[i] + I[i] + R[i] #Total NHP population by genus
inv_P[1:n_gen] <- 1.0/P[i]

#Updates to output values at each time increment--------------------------------
update(FOI_total) <- FOI_sum
update(S[1:n_gen]) <- max(Pmin,S[i] - E_new[i] + births[i] - (deaths_nat[i]*S[i]*inv_P[i]))
update(E[1:n_gen]) <- max(Pmin,E[i] + E_new[i] - I_new[i] - (deaths_nat[i]*E[i]*inv_P[i]))
update(I[1:n_gen]) <- max(Pmin,I[i] + I_new[i] - R_new0[i] - (deaths_nat[i]*I[i]*inv_P[i]))
update(R[1:n_gen]) <- max(Pmin,R[i] + R_new[i] - (deaths_nat[i]*R[i]*inv_P[i]))
update(C[1:n_gen]) <- I_new[i]
update(D[1:n_gen]) <- deaths_YF[i]

#Initial values-----------------------------------------------------------------
initial(FOI_total) <- FOI_in[1]
initial(S[1:n_gen]) <- S_0[i]
initial(E[1:n_gen]) <- E_0[i]
initial(I[1:n_gen]) <- I_0[i]
initial(R[1:n_gen]) <- R_0[i]
initial(C[1:n_gen]) <- 0
initial(D[1:n_gen]) <- 0

#Dimensions---------------------------------------------------------------------
dim(S) <- n_gen
dim(E) <- n_gen
dim(I) <- n_gen
dim(R) <- n_gen
dim(C) <- n_gen
dim(D) <- n_gen

dim(births) <- n_gen
dim(deaths_nat) <- n_gen
dim(E_new) <- n_gen
dim(I_new) <- n_gen
dim(R_new0) <- n_gen
dim(R_new) <- n_gen
dim(deaths_YF) <- n_gen

dim(P) <- n_gen
dim(inv_P) <- n_gen

dim(S_0) <- n_gen
dim(E_0) <- n_gen
dim(I_0) <- n_gen
dim(R_0) <- n_gen

dim(beta) <- n_t_pts
dim(FOI_in) <- n_t_pts
dim(t_latent) <- n_gen
dim(t_infectious) <- n_gen
dim(ifr) <- n_gen
dim(mu) <- n_gen
dim(rate1) <- n_gen
dim(rate2) <- n_gen
