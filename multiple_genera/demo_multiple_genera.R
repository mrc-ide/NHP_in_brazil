library(odin.dust)

model <- odin.dust::odin_dust("multiple_genera/sir_multiple_genera.R")

dt=1 #Time increment in days
beta=0.075 #Beta parameter governing transmission from infectious individuals
t_incubation=5 #Incubation time in vectors (currently no vector species separation)
t_latent=c(5,5) #Latent period in NHPs by genus
t_infectious=c(5,5)  #Infectious period in NHPs by genus
FOI_in=1.0e-9 #Force of infection (per day) for importation of infections from outside modelled populations
R_dP=c(100,200)/365 #Population turnover (birth rate/death rate per day) by genus (assumed stable population)
n_gen=2 #Number of genera (vectors of input values must have lengths equal to n_gen)
S_0=c(1000,1000) #Initial number of susceptible individuals by genus
E_0=c(0,0) #Initial number of exposed individuals by genus
I_0=c(0,0) #Initial number of infectious individuals by genus
R_0=c(0,0) #Initial number of recovered individuals by genus
pars <- list(dt=dt,beta=beta,t_incubation=t_incubation,t_latent=t_latent,t_infectious=t_infectious,
             FOI_in=FOI_in,R_dP=R_dP,n_gen=n_gen,S_0=S_0,E_0=E_0,I_0=I_0,R_0=R_0)

tmax=36500
n_particles=1
run <- model$new(pars = pars, time=1, n_particles=n_particles, n_threads=1, deterministic=TRUE)
output=array(NA,dim=c(tmax,n_particles,12))
for(i in 1:tmax){
  output[i,,]=run$run(i)
}

time=output[,1,1]
FOI=output[,1,2]
S=output[,1,2+c(1:pars$n_gen)]
E=output[,1,2+pars$n_gen+c(1:pars$n_gen)]
I=output[,1,2+(2*pars$n_gen)+c(1:pars$n_gen)]
R=output[,1,2+(3*pars$n_gen)+c(1:pars$n_gen)]
C=output[,1,2+(3*pars$n_gen)+c(1:pars$n_gen)]

par(mar=c(4,4,1,1))
FOI_ticks=10^c(-8:-1)
matplot(x=time,y=log(FOI),type="l",col=1,xlab="Time",ylab="Total force of infection",yaxt="n")
axis(side=2,at=log(FOI_ticks),labels=FOI_ticks)

par(mfrow=c(1,pars$n_gen),mar=c(4,2,1,1))
for(i in 1:pars$n_gen){
  matplot(x=c(1,max(time)),y=c(0,pars$S_0[i]+pars$E_0[i]+pars$I_0[i]+pars$R_0[i]),type="p",col=0,xlab="Time",ylab="")
  matplot(x=time,y=S[,i],type="l",col=1,add=TRUE)
  matplot(x=time,y=R[,i],type="l",col=2,add=TRUE)  
  legend("topright",legend=c("S","R"),lty=c(1,1),col=c(1,2))
}
par(mfrow=c(1,1))