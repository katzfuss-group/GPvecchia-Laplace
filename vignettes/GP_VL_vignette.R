rm(list = ls())
#setwd("/Users/danielzilber/Desktop/Katzfuss/code_CPP/")

###  load GPvecchia package
# currently, we don't have an actual package yet
# so the "package" is just a bunch of functions in the folder GPvecchia
# library(GPvecchia)
for (nm in list.files("R",pattern = "\\.[RrSsQq]$")) {
  cat(nm,":"); source(file.path("R",nm)); cat("\n")
}
Rcpp::sourceCpp('src/NZentries_arma_omp1.cpp')


#####################  data generation  #######################

dimen=1 # number of spatial dimensions
samp_size=10^2  # number of observed locs
domn = 1# domain over which to generate points

sig2=1; range=.1; smooth=1.5
covparms=c(sig2,range,smooth)

covfun <- function(locs) sig2*Matern(rdist(locs),range=range,smoothness=smooth)

# current data generation is over a uniform lattice
ls=logistic_sample(samp_size,covfun, seed = 12, dom = domn, dimen=dimen)



#####################  posterior generation  #######################
neighbor_count= 4
likelihood_model = ls
posterior = calculate_posterior_VL(likelihood_model, covparms, m=neighbor_count)


#####################  Compare posterior to latent  #######################
plot(likelihood_model$locs, likelihood_model$y)
points(likelihood_model$locs, posterior$mean, type= "l", col=3)
