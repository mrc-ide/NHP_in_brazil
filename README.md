# NHP_in_brazil

This repository will contain relevant examples of using odin.dust and mcstate for defining and estimating a model of transmission in non-human primates. These examples are based on the SIR model shared with odin/odin.dust and mcstate packages.

# SIR example
This is almost exactly the same as the vignette for mcstate. Here we have 4 files in a folder- the index, compare, sir and example scripts. The index and compare are the same as the mcstate example and are functions to pick out the relevant points from the data/model, or, perform the likelihood/ observation process calculation respectively. The sir model is a susceptible-infectious-recovered model with a 1 day timestep as default. In the example.R script is where you find an example of running the particle filter on this model, and an example of estimating beta and gamma (transmission parameter and recovery parameter respectively) from dummy data (also taken from the MCstate package).

# SIR_temp example
This is almost identical to the above but the transmission parameter beta now depends nonlinearly on temperature. It means that rather than estimating beta directly, you estimate of the parameters of the function for how beta varies with temperature. This will be more useful as other regions/ municipalities are added.
