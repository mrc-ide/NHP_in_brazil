library(assertthat)
t_incubation=5.0
t_latent=5.0
t_infectious=5.0

Human_Model_Run <- function(FOI_spillover = c(),R0 = c(),vacc_data = list(),pop_data = list(),years_data = c(1940:1941),
                      start_SEIRV = list(), output_type = "full", year0 = 1940, mode_start = 0,
                      vaccine_efficacy = 1.0, dt = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE) {
  
  #TODO Add assert_that functions (NB - Some checks carried out in human_pars_setup)
  assert_that(n_particles<=20,msg="Number of particles must be 20 or less")
  
  n_nv=3 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  n_data_pts=(6*N_age)+n_nv #Number of data values per time point in output
  step_begin=((years_data[1]-year0)*(365/dt)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/dt))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data
  
  x <- OD_human$new(pars=human_pars_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,
                                            vaccine_efficacy,start_SEIRV,dt),
                       time = 0, n_particles = n_particles, n_threads = n_threads, deterministic = deterministic)
  
  x_res <- array(NA, dim = c(n_data_pts, n_particles, t_pts_out))
  for(step in step_begin:step_end){
    x_res[,,step-step_begin+1] <- x$run(step)
  }
  if(step_begin==0){x_res[2,,1]=rep(year0,n_particles)}
  
  if(output_type=="full"){
    dimensions=c(N_age,n_particles,t_pts_out)
    output_data=list(day=x_res[1,1,],year=x_res[2,1,])
    output_data=list(day=x_res[1,1,],year=x_res[2,1,])
    output_data$FOI_total=array(x_res[3,,]/dt,dim=c(n_particles,t_pts_out))
    output_data$S=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=dimensions)
    output_data$E=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=dimensions)
    output_data$I=array(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),,],dim=dimensions)
    output_data$R=array(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),,],dim=dimensions)
    output_data$V=array(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),,],dim=dimensions)
    output_data$C=array(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),,],dim=dimensions)
  } else {
    if(output_type=="case_alt2"){
      output_data=list(day=x_res[1,1,],year=x_res[2,1,])
      output_data$C=array(0,dim=c(n_particles,t_pts_out))
      for(pt in 1:t_pts_out){
        for(n_p in 1:n_particles){
          output_data$C[n_p,pt]=sum(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),n_p,pt])
        }
      }
    }  else {
      n_years=length(years_data)
      output_data=list(year=years_data)
      if(output_type=="case+sero" || output_type=="sero"){
        output_data$V=output_data$R=output_data$I=output_data$E=output_data$S=array(0,dim=c(N_age,n_particles,n_years))
        for(n_year in 1:n_years){
          pts=c(1:t_pts_out)[x_res[2,1,]==years_data[n_year]]
          for(n_p in 1:n_particles){
            output_data$S[,n_p,n_year]=rowMeans(x_res[c((1+n_nv):(N_age+n_nv)),n_p,pts])
            output_data$E[,n_p,n_year]=rowMeans(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),n_p,pts])
            output_data$I[,n_p,n_year]=rowMeans(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),n_p,pts])
            output_data$R[,n_p,n_year]=rowMeans(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),n_p,pts])
            output_data$V[,n_p,n_year]=rowMeans(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),n_p,pts])
          }
        }
      }
      if(output_type=="case+sero" || output_type=="case"){
        output_data$C=array(0,dim=c(n_particles,n_years))
        for(n_year in 1:n_years){
          pts=c(1:t_pts_out)[x_res[2,1,]==years_data[n_year]]
          for(n_p in 1:n_particles){
            output_data$C[n_p,n_year]=sum(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),n_p,pts])
          }
        }
      }
      if(output_type=="case_alt"){
        output_data$C=array(0,dim=c(N_age,n_particles,n_years))
        for(n_year in 1:n_years){
          pts=c(1:t_pts_out)[x_res[2,1,]==years_data[n_year]]
          for(n_p in 1:n_particles){
            output_data$C[,n_p,n_year]=rowSums(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),n_p,pts])
          }
        }
      }
    }
  }
  
  return(output_data)
}
#-------------------------------------------------------------------------------
human_pars_setup <- function(FOI_spillover=c(),R0=c(),vacc_data=list(),pop_data=list(),year0=1940,
                            years_data=c(1941:1942),mode_start=0,vaccine_efficacy=1.0,start_SEIRV=list(),dt=1.0){
  
  assert_that(length(pop_data[,1])>1,msg="Need population data for multiple years")
  assert_that(length(pop_data[1,])>1,msg="Need population data for multiple age groups")
  n_years=length(pop_data[,1])-1
  N_age=length(pop_data[1,])
  assert_that(length(vacc_data[,1])==n_years+1,msg="Population and vaccination data must be for same time periods")
  assert_that(length(vacc_data[1,])==N_age,msg="No. age groups in population and vaccination data must match")
  assert_that(mode_start %in% c(0,1,2,3),msg="mode_start must have value 0, 1, 2 or 3")
  assert_that(vaccine_efficacy<=1.0 && vaccine_efficacy>=0.0,msg="Vaccine efficacy must be between 0 and 1")
  if(mode_start==2){
    assert_that(is.null(start_SEIRV$S)==FALSE,msg="When mode_start=2, start_SEIRV data is required")
  }
  assert_that(years_data[1]>=year0,msg="First data year must be greater than or equal to year0")
  assert_that(max(years_data)+1-year0<=n_years,msg="Period of years_data must lie within population data")
  vacc_initial=vacc_data[1,]
  assert_that(dt %in% c(1,2.5,5),msg="dt must have value 1, 2.5 or 5 days (must have integer no. points/year)")
  inv_365=1.0/365.0
  assert_that(length(FOI_spillover)==n_years*(365/dt),msg="Spillover FOI vector must have length equal to no. time points")
  assert_that(length(R0)==n_years*(365/dt),msg="R0 vector must have length equal to no. time points")
  
  P0=S_0=E_0=I_0=R_0=V_0=rep(0,N_age)
  dP1_all=dP2_all=vacc_rates=array(NA,dim=c(N_age,n_years))
  for(i in 1:N_age){
    P0[i]=max(1.0,pop_data[1,i]) #Set all population values to nonzero minimum to avoid NaN values
  }
  for(n_year in 1:n_years){
    for(i in 1:N_age){
      dP1_all[i,n_year]=max(1.0,pop_data[n_year+1,i])*inv_365
      dP2_all[i,n_year]=max(1.0,pop_data[n_year,i])*inv_365
      if(i==1){
        vacc_rates[i,n_year]=vacc_data[n_year+1,i]*inv_365
      } else {
        vacc_rates[i,n_year]=max(0.0,vacc_data[n_year+1,i]-vacc_data[n_year,i-1])*inv_365
      }
    }
  }
  
  #-----------------------------------------------------------------------------
  if(mode_start==2){
    S_0=start_SEIRV$S
    E_0=start_SEIRV$E
    I_0=start_SEIRV$I
    R_0=start_SEIRV$R
    V_0=start_SEIRV$V
  } else {
    V_0=P0*vacc_initial
    #-----------------------------------------------------------------------------
    if(mode_start==0){
      S_0=P0*(1.0-vacc_initial)
    }
    #-----------------------------------------------------------------------------
    if(mode_start==1){ #Herd immunity, uniform by age
      if(R0>1.0){
        herd_immunity=1.0-(1.0/R0)
      } else {
        herd_immunity=0.0
      }
      for(i in 1:N_age){
        if(vacc_initial[i]<herd_immunity){
          R_0[i]=P0[i]*(herd_immunity-vacc_initial[i])
          S_0[i]=P0[i]*(1.0-herd_immunity)
        } else {
          S_0[i]=P0[i]*(1.0-vacc_initial[i])
        }
      }
    }
    #-----------------------------------------------------------------------------
    if(mode_start==3){ #New herd immunity calculation to give age-stratified immunity profile based on notional FOI
      ages=c(1:N_age)-1
      if(R0[1]<=1.0){
        FOI_estimate=FOI_spillover*365.0
      } else {
        estimation_results=nlm(imm_fraction_function,p=-4,R0[1],ages,P0/sum(P0))
        FOI_estimate=min(0.1,(FOI_spillover*365.0)+exp(estimation_results$estimate))
      }
      herd_immunity=1.0-(exp(-FOI_estimate*(ages+0.5)))
      
      for(i in 1:N_age){
        if(vacc_initial[i]<herd_immunity[i]){
          R_0[i]=P0[i]*(herd_immunity[i]-vacc_initial[i])
          S_0[i]=P0[i]*(1.0-herd_immunity[i])
        } else {
          S_0[i]=P0[i]*(1.0-vacc_initial[i])
        }
      }
    }
  }
  
  return(list(FOI_spillover=FOI_spillover,R0=R0,vacc_rate_daily=vacc_rates,N_age=N_age,
              S_0=S_0,E_0=E_0,I_0=I_0,R_0=R_0,V_0=V_0,dP1_all=dP1_all,dP2_all=dP2_all,n_years=n_years,
              year0=year0,vaccine_efficacy=vaccine_efficacy,dt=dt,
              t_incubation=t_incubation,t_latent=t_latent,t_infectious=t_infectious))
}