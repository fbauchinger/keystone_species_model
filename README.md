# keystone_species_model
calculate training data and train a linear model to identify keystone species

adapted from "Deciphering microbial interactions and detecting keystone species with co-occurrence networks" (Berry and Widder; 2014)
https://doi.org/10.3389/fmicb.2014.00219

make_community.R
runs Lotka-Volterra multispecies model

keystone_runs.R
re-runs multispecies Lotka-Volterra model as specified in make_community.R
in each run, abundance of one species is set to 0
the distance between the initial model run, with the species present, and the model run with the species absent, is calculated
this procedure is repeated multiple times and an average distance for each species is calculated
this average distance is taken as a proxy for a species keystone potential

make_network.R
compute correlation networks based on output from Lotka-Volterra multispecies model (make_community.R)
correlation coefficients are compared to a null distribution to assess significance

linear_model.R
train a linear model to predict keystoneness of species based on correlation network statistics

trainingData.csv
training data used for linear model

trainingData_model_stats.csv
specification of input parameters used for model runs
