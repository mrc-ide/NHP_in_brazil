compare <- function(state, observed, pars = NULL) {
  if (is.na(observed$cases)) {
    return(NULL)
  }
  exp_noise <- 1e6
  incidence_modelled <- state[5, , drop = TRUE]
  incidence_observed <- observed$cases
  lambda <- incidence_modelled +
    rexp(n = length(incidence_modelled), rate = exp_noise)
  dpois(x = incidence_observed, lambda = lambda, log = TRUE)
}