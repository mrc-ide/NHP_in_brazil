compare <- function(state, observed, pars = NULL) {
  if (any(is.na(observed$cases))) {
    return(NULL)
  } else {
    exp_noise <- 1e6
    incidence_modelled <- state[5, , drop = TRUE]
    incidence_observed <- observed$cases
    lambda <- incidence_modelled +
      rexp(n = length(incidence_modelled), rate = exp_noise)
    return(dpois(x = incidence_observed, lambda = lambda, log = TRUE))
  }
}