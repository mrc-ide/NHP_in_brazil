# Some ideas for expanding the model examples 

## Introduction

We have included a few different examples in this repository on how to start estimating an SIR model. However, it is likely that SIR is not the right structure, and that more environmental covariates are needed. This doc talks to some possible extensions.

## Other model structures

By this, we mean including other compartments. We have only included SIR models here for simplicity but it is likely that it will be important to include a latent or exposed class as well as disease related mortality and maybe even modelling vector populations explicitly. 

### Including vectors

There are stepping stones for including vectors. In an SIR model, they may be included implicitly in the definition of the transmission parameter. They may also be included in some of the calculations of the incubation periods. The most complicated option is to include another set of compartments to track infection in vectors; however, this would be quickly complicated by including different species and require a variety of new data. It would also make estimation of the model significantly more complex.

## Other environmental variables

In the examples here we include only temperature dependence through a simple Briere function (see [Mordecai et al. 2017](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0005568) for functional form and [Gaythorpe et al. 2020](https://elifesciences.org/articles/55619) for application/estimation for yellow fever).

However, it may be important to include other environmental variables such as precipitation, vegetation etc. These can be included in a similar way ie. through the transmission parameter and we would advise, that where there is no information on a fucntional form for the relationship, using a simple linear relationship as a starting point.

