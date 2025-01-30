N <- parameter(1000) 
I0 <- parameter(10)
c_beta_temp <- parameter() #Linear coefficient used to calculate beta from temperature
c_gamma_temp <- parameter() #Linear coefficient used to calculate gamma from temperature
temp <- parameter() #Time-varying temperature

t_pt <- time/dt
beta <- c_beta_temp*temp[t_pt]
gamma <- c_gamma_temp*temp[t_pt]
p_SI <- 1 - exp(-beta * I / N * dt)
p_IR <- 1 - exp(-gamma * dt)
n_SI <- Binomial(S, p_SI) 
n_IR <- Binomial(I, p_IR)

update(time) <- time + dt
update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
update(incidence) <- n_SI

initial(time) <- 0
initial(S) <- N - I0
initial(I) <- I0
initial(R) <- 0
initial(incidence) <- 0

# Adding goodness of fit to model 
# Only used when comparing to observed data
cases <- data() #Data to compare
cases ~ Poisson(incidence) #Distribution
