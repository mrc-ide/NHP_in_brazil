# SEIR (Susceptible, Exposed, Infectious, Recovered) yellow fever model for non-human primates (NHPs),
# with multiple genera

# TO DO (to add or to include in future version):
# -Update population modelling to allow for population extinction or recovery (look at existing ecological models, ongoing)
# -Incorporate introductions of individuals from outside (and outward migration)
# -Add vector model including functionality for multiple vector species (framework in place, need parameters)
# -Change model of force of infection, latency, etc. for better short-term timescale modelling and better 
#  modelling of multiple interacting populations (Diekmann et al 2010)
# -Make coefficients governing calculation of epi parameters from temp/precip genus dependent
# -Make calculation of epi parameters from temp/precip non-linear

#Parameters---------------------------------------------------------------------
dt <- user() #Time increment in days
c_beta_temp <- user() #Coefficient used to calculate temperature component of beta parameter
c_beta_precip <- user() #Coefficient used to calculate precipitation component of beta parameter
c_FOI_in_temp <- user() #Coefficient used to calculate temperature component of FOI_in parameter
c_FOI_in_precip <- user() #Coefficient used to calculate precipitation component of FOI_in parameter
t_incubation <- user() #Length in days of yellow fever incubation period in mosquito vectors [TODO: MAKE VAR BY SPECIES]
#
n_t_pts <- user() #Number of time points to be run
temp[] <- user() #Time-varying temperature used to calculate beta and FOI_in
precip[] <- user() #Time-varying precipitation used to calculate beta and FOI_in
#
n_gen <- user() #Number of genera
t_latent[] <- user() #Length in days of latent period in NHPs exposed to yellow fever, by genus
t_infectious[] <- user() #Length of infectious period in NHPs with yellow fever, by genus
ifr[] <- user() #Infection fatality rate by genus
mu[] <- user() #Rate of population turnover by day (birth rate/death rate, assumed equal for stable population) by genus
S_0[] <- user() #Susceptible NHP population by genus at start
E_0[] <- user() #Exposed NHP population by genus at start
I_0[] <- user() #Infectious NHP population by genus at start
R_0[] <- user() #Recovered NHP population by genus at start

#Additional/derived values (constant)
Pmin <- 0.0 #Minimum population setting to avoid negative numbers
FOI_max <- 1.0 #Upper threshold for total force of infection to avoid more infections than individuals
rate1[1:n_gen] <- dt/(t_incubation+t_latent[i]) #Rate of transference E->I
rate2[1:n_gen] <- dt/t_infectious[i] #Rate of transference I->R

#Derived values (variable)
beta <- (c_beta_temp*temp[step+1])+(c_beta_precip*precip[step+1]) #Time-varying beta parameter for transmission
FOI_in <- (c_FOI_in_temp*temp[step+1])+(c_FOI_in_precip*precip[step+1]) #Time-varying force of infection (per day) for outside introduction
FOI_sum <-  min(FOI_max,dt*(((beta*sum(I))/sum(P)) + FOI_in)) #Total force of infection
E_new[1:n_gen] <- rbinom(as.integer(S[i]), FOI_sum) #New exposed NHPs by genus
I_new[1:n_gen] <- E[i]*rate1[i]     #New infectious NHPs by genus
R_new0[1:n_gen] <- I[i]*rate2[i]     #Number of NHPS removed from I by genus (recovery from or death due to YF)
R_new[1:n_gen] <- R_new0[i]*(1.0-ifr[i]) #New recovered NHPs by genus
deaths_YF[1:n_gen] <- R_new0[i]*ifr[i]  #New deaths of NHPs by genus due to YF
P[1:n_gen] <- S[i] + E[i] + I[i] + R[i] #Total NHP population by genus

#Updates to output values at each time increment--------------------------------
update(time) <- time + dt
update(FOI_total) <- FOI_sum
update(S[1:n_gen]) <- max(Pmin,S[i] - E_new[i] + rbinom(as.integer(P[i]), mu[i]) - rbinom(as.integer(S[i]), mu[i]))
update(E[1:n_gen]) <- max(Pmin,E[i] + E_new[i] - I_new[i] - rbinom(as.integer(E[i]), mu[i]))
update(I[1:n_gen]) <- max(Pmin,I[i] + I_new[i] - R_new0[i] - rbinom(as.integer(I[i]), mu[i]))
update(R[1:n_gen]) <- max(Pmin,R[i] + R_new[i] - rbinom(as.integer(R[i]), mu[i]))
update(C[1:n_gen]) <- I_new[i]
update(D[1:n_gen]) <- deaths_YF[i]

#Initial values of outputs------------------------------------------------------
initial(time) <- 0
initial(FOI_total) <- 0
initial(S[1:n_gen]) <- S_0[i]
initial(E[1:n_gen]) <- E_0[i]
initial(I[1:n_gen]) <- I_0[i]
initial(R[1:n_gen]) <- R_0[i]
initial(C[1:n_gen]) <- 0
initial(D[1:n_gen]) <- 0

#Dimensions - inputs------------------------------------------------------------
dim(temp) <- n_t_pts
dim(precip) <- n_t_pts
dim(S_0) <- n_gen
dim(E_0) <- n_gen
dim(I_0) <- n_gen
dim(R_0) <- n_gen
dim(t_latent) <- n_gen
dim(t_infectious) <- n_gen
dim(ifr) <- n_gen
dim(mu) <- n_gen

#Dimensions - derived values----------------------------------------------------
dim(rate1) <- n_gen
dim(rate2) <- n_gen
dim(E_new) <- n_gen
dim(I_new) <- n_gen
dim(R_new0) <- n_gen
dim(R_new) <- n_gen
dim(deaths_YF) <- n_gen
dim(P) <- n_gen

#Dimensions - outputs-----------------------------------------------------------
dim(S) <- n_gen
dim(E) <- n_gen
dim(I) <- n_gen
dim(R) <- n_gen
dim(C) <- n_gen
dim(D) <- n_gen
