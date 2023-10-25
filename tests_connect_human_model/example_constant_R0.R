# a script for examining fitting a NHP in one region
# loading libraries
library(mcstate)
library(odin.dust)
library(YEPaux)

# sir model
OD_NHPs <- odin.dust::odin_dust("tests_connect_human_model/OD_NHPs.R") # Non-human primate SIR model
OD_human <- odin.dust::odin_dust("tests_connect_human_model/OD_human.R") # Human SEIRV model
source("tests_connect_human_model/R_human.R") #R functions for running human model

#Common parameters
years_data=c(1940:2000) #Years of interest
dt=1.0 #Length of interval between time points in days
n_pts=length(years_data)*(365.0/dt) #Number of time points to run

#Parameters - NHPs
R0_NHPs <- rep(1.3,n_pts)
pop_dens <- 10 # based on density in Culot et al. Botacatu and others
area <- 248219 #based on Sao Paulo state
pop_size <- pop_dens * area #assume similar across the region *assumption*
pars_NHPs <- list(R0 = R0_NHPs, gamma = 0.1, S0 = pop_size, dt = dt, n_pts = n_pts) #leave others as defaults

#Run NHP model and save and plot output
{
  mod_NHPs <- OD_NHPs$new(pars=pars_NHPs,time=0,n_particles=1,n_threads=1,deterministic=TRUE) 
  mod_NHPs_out=data.frame(time=rep(NA,n_pts),S=rep(NA,n_pts),I=rep(NA,n_pts),R=rep(NA,n_pts))
  for(i in 1:n_pts){
    output_step=mod_NHPs$run(i)
    mod_NHPs_out[i,]=output_step[c(1,2,4,3),1]
  }
  
  matplot(x=mod_NHPs_out$time[c(1,n_pts)],y=c(0,max(mod_NHPs_out[c(2:4),])),col=0,xlab="Time",ylab="")
  matplot(x=mod_NHPs_out$time,y=mod_NHPs_out$S,type="l",col=1,add=TRUE)
  matplot(x=mod_NHPs_out$time,y=mod_NHPs_out$I,type="l",col=2,add=TRUE)
  matplot(x=mod_NHPs_out$time,y=mod_NHPs_out$R,type="l",col=3,add=TRUE)
  legend("bottomleft",legend=c("S","I","R"),col=c(1:3),lty=c(1,1,1))
  
  NHP_I_fraction=mod_NHPs_out$I/pop_size
  matplot(x=mod_NHPs_out$time,y=NHP_I_fraction,type="l",xlab="Time",ylab="NHP infectious fraction")  
}

#Parameters - human 
FOI_spillover_min=1.0e-8 #Baseline FOI_spillover to which value calculated from NHPs added
FOI_coeff=1e-5 #Coefficient by which infectious fraction of NHPs multiplied
FOI_spillover=(NHP_I_fraction*FOI_coeff)+FOI_spillover_min #Total spillover FOI
R0_human=rep(1.0,n_pts)
year0=years_data[1]
human_input_data <- readRDS(file = "tests_connect_human_model/input_data_example.Rds")
vacc_data <- human_input_data$vacc_data[1, human_input_data$years_labels %in% c(year0:(max(years_data)+1)), ]
pop_data <- human_input_data$pop_data[1, human_input_data$years_labels %in% c(year0:(max(years_data)+1)), ]
mode_start=0
vaccine_efficacy=1.0
start_SEIRV=NULL
dt=1.0
output_type="full"

#Run human model
human_model_data <- Human_Model_Run(FOI_spillover,R0_human,vacc_data,pop_data,years_data,start_SEIRV,output_type,
                                    year0,mode_start,vaccine_efficacy,dt,n_particles=1,n_threads=1,deterministic=TRUE)
YEPaux::plot_model_output(human_model_data)

#Plot spillover and total FOI
dates=year0+(c(1:n_pts)*(dt/365))
matplot(x=dates,y=log(FOI_spillover),type="l",col=1,xlab="Date",ylab="Force of infection",
        ylim=c(log(min(FOI_spillover)),log(max(human_model_data$FOI_total[1,]))),yaxt="n")
matplot(x=dates,y=log(human_model_data$FOI_total[1,]),type="l",col=2,add=TRUE)
axis(side=2,at=log(10^(c(-8:0))),labels=10^c(-8:0))
legend("topright",legend=c("Spillover","Total"),col=c(1,2),lty=c(1,1))
