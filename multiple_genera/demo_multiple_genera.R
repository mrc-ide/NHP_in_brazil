# Demonstration of SEIR non-human primate model with multiple genera
# Epidemiological (epi) parameters calculated linearly from temperature and precipitation:
# - beta - standard beta parameter governing new infections due to infectious individuals
# - FOI_in - daily force of infection for new infections imported from outside
# Epi parameters calculated from linear combination of temperature and precipitation using coefficients
# Future version will incorporate non-linear calculation of epi parameters from temperature and precipitation

#Load odin.dust and compile model-----------------------------------------------
library(odin.dust)
model <- odin.dust::odin_dust("multiple_genera/seir_multiple_genera.R")

#Set up parameters--------------------------------------------------------------
dt=1 #Time increment in days
n_t_pts=7300 #Number of time points to run (beta and FOI_in must have lengths equal to n_t_pts)
temp=rep(25.0,n_t_pts) #Temperature at each time point, to be used to calculate epi parameters beta and FOI_in
precip=rep(100,n_t_pts) #Precipitation at each time point, to be used to calculate epi parameters beta and FOI_in
c_beta_temp=1e-4 #Coefficient used to calculate temperature component of beta parameter
c_beta_precip=1e-6 #Coefficient used to calculate precipitation component of beta parameter
c_FOI_in_temp=1e-4 #Coefficient used to calculate temperature component of FOI_in parameter
c_FOI_in_precip=1e-6 #Coefficient used to calculate precipitation component of FOI_in parameter
t_incubation=5 #Incubation time in vectors (currently no vector species separation)
t_latent=c(5,5) #Latent period in NHPs by genus
t_infectious=c(5,5)  #Infectious period in NHPs by genus
ifr=c(0.1,0.5) #Infection fatality rate by genus
mu=c(0.1,0.1)/365 #Birth/death rate per unit population per day (1/lifespan in years divided by 365) by genus
n_gen=2 #Number of genera (vectors of input values must have lengths equal to n_gen)
S_0=c(1000,1000) #Initial number of susceptible individuals by genus
E_0=c(0,0) #Initial number of exposed individuals by genus
I_0=c(0,0) #Initial number of infectious individuals by genus
R_0=c(0,0) #Initial number of recovered individuals by genus
pars <- list(dt = dt, n_t_pts = n_t_pts, temp = temp, precip = precip, 
             c_beta_temp = c_beta_temp, c_beta_precip = c_beta_precip, c_FOI_in_temp = c_FOI_in_temp, c_FOI_in_precip = c_FOI_in_precip, 
             t_incubation = t_incubation, t_latent = t_latent, t_infectious = t_infectious, 
             ifr = ifr, mu = mu, n_gen = n_gen, S_0 = S_0, E_0 = E_0, I_0 = I_0, R_0 = R_0)

n_particles=1
run <- model$new(pars = pars, time=1, n_particles=n_particles, n_threads=1, deterministic=TRUE)
output=array(NA,dim=c(pars$n_t_pts,n_particles,2+(6*pars$n_gen)))
for(i in 1:pars$n_t_pts){
  output[i,,]=run$run(i)
}

time=output[,1,1]
FOI=output[,1,2]
S=output[,1,2+c(1:pars$n_gen)]
E=output[,1,2+pars$n_gen+c(1:pars$n_gen)]
I=output[,1,2+(2*pars$n_gen)+c(1:pars$n_gen)]
R=output[,1,2+(3*pars$n_gen)+c(1:pars$n_gen)]
C=output[,1,2+(4*pars$n_gen)+c(1:pars$n_gen)]
D=output[,1,2+(5*pars$n_gen)+c(1:pars$n_gen)]
C_c=D_c=C
for(i in 1:pars$n_gen){
  C_c[,i]=cumsum(C[,i])
  D_c[,i]=cumsum(D[,i])
}

par(mar=c(4,4,1,1))
FOI_ticks=10^c(-8:-1)
matplot(x=time,y=log(FOI),type="l",col=1,xlab="Time",ylab="Total force of infection",yaxt="n")
axis(side=2,at=log(FOI_ticks),labels=FOI_ticks)

par(mfrow=c(1,pars$n_gen),mar=c(4,2,1,1))
for(i in 1:pars$n_gen){
  matplot(x=c(1,max(time)),y=c(0,1.2*(pars$S_0[i]+pars$E_0[i]+pars$I_0[i]+pars$R_0[i])),type="p",col=0,xlab="Time",ylab="")
  matplot(x=time,y=S[,i],type="l",col=1,add=TRUE)
  matplot(x=time,y=R[,i],type="l",col=2,add=TRUE)  
  matplot(x=time,y=S[,i]+E[,i]+I[,i]+R[,i],type="l",col=3,add=TRUE)  
  legend("topright",legend=c("Susceptible","Recovered","Total population"),lty=c(1,1,1),col=c(1,2,3))
}
par(mfrow=c(1,1))

par(mfrow=c(1,pars$n_gen),mar=c(4,2,1,1))
for(i in 1:pars$n_gen){
  matplot(x=c(1,max(time)),y=c(0,max(C_c[,i])),type="p",col=0,xlab="Time",ylab="")
  matplot(x=time,y=C_c[,i],type="l",col=1,add=TRUE)
  matplot(x=time,y=D_c[,i],type="l",col=2,add=TRUE)  
  legend("topleft",legend=c("Infections","Deaths"),lty=c(1,1),col=c(1,2))
}
par(mfrow=c(1,1))