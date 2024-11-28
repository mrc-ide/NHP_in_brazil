N <- parameter(1000) # In odin2, parameter() is used for inputs where user() was used in odin
I0 <- parameter(10)
beta <- parameter(0.2)
gamma <- parameter(0.1)

p_SI <- 1 - exp(-beta * I / N * dt) # Note parameter dt - this is a time interval set when the model is set up via dust_system_create
p_IR <- 1 - exp(-gamma * dt)
n_SI <- Binomial(S, p_SI) # In odin2, Binomial() is used instead of rbinom()
# See https://mrc-ide.github.io/odin2/articles/functions.html#distribution-functions for available distributions
n_IR <- Binomial(I, p_IR)

update(S) <- S - n_SI
update(I) <- I + n_SI - n_IR
update(R) <- R + n_IR
update(incidence) <- n_SI

initial(S) <- N - I0
initial(I) <- I0
initial(R) <- 0
initial(incidence) <- 0

# Adding goodness of fit to model (in odin2, this is implemented within the model code)
# Only used when comparing to observed data
cases <- data() #Data to compare
cases ~ Poisson(incidence) #Distribution