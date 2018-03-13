library(fields)

# helper method to create a uniform lattice over arbitrary dimension
# Probably implemented elsewhere, should check
create_locs = function(n_tot=100, dimen=1, dom=1){
  #n=5; dimen=3; dom=1
  # for higher dimension, generate the right number of locations
  n=round(n_tot^(1/dimen))
  if (n_tot^(1/dimen) != n){
    # Number of points doesnt yield an equal sided lattice for dimen
    stop
  }
  locs = seq(0,dom,length.out=n)
  loc_hd = matrix(0, nrow = n^dimen, ncol=dimen)
  for (d in 1:dimen){
    col_for_dim = c()
    for (ele in locs){
      col_for_dim = c(col_for_dim, rep(ele, n^(d-1)))
    }
    # change index just to match intuition
    loc_hd[,dimen-d+1] = rep(col_for_dim, n^(dimen-d))
  }
  return(loc_hd)
}



#################  LOGISTIC functions #########################
logistic_sample = function(n, covfun, seed = 125, dom = 1, dimen = 1){
  #n=400; dimen=2; dom=1
  set.seed(seed)
  # create locations for a domain [0, dom]
  locs = create_locs(n, dimen, dom)
  C = covfun(locs)
  y = t(chol(C))%*%rnorm(n)

  # generate observations with link
  z = rbinom(n,1,prob = exp(y)/(1+exp(y)))
  # store score and hessian functions
  logistic_llh = function(y_o, z) sum(z*y_o-log(1+exp(y_o)))
  logistic_hess = function(y_o, z) diag(array(exp(y_o)/(1+exp(y_o))^2))
  logistic_score = function(y_o, z) z - exp(y_o)/(1+exp(y_o))


  # return object with all components of model
  return(list("type"="logistic", "locs" =  matrix(locs,ncol=dimen), "z"=z, "y"=y,  "C"=C, "hess" = logistic_hess, "score"=logistic_score,"llh" = logistic_llh))
}

#################  Poisson function #########################
pois_sample = function(n, covfun, seed_val = 124, dom=1, dimen=1){

  set.seed(seed_val)
  locs = create_locs(n, dimen, dom)

  C = covfun(locs)
  y = t(chol(C))%*%rnorm(n)
  #link
  z = rpois(n, exp(y))
  pois_llh = function(y_o, z) sum(z*y_o -exp(y_o)-log(factorial(z)))
  pois_hess =function(y_o, z) diag(array(exp(y_o)))
  pois_score = function(y_o, z) z-exp(y_o)
  return(list("type"="Poisson", "locs" = matrix(locs,ncol=dimen), "z"=z, "y"=y, "C"=C, "hess" = pois_hess, "score"=pois_score,"llh" = pois_llh))
}

#################  Gaussian functions #########################
gauss_sample = function(n, covfun, seed_val=125, dom=1, dimen=1 , sigma = .3){
  set.seed(seed_val)
  locs = create_locs(n, dimen, dom)
  C = covfun(locs)
  y = t(chol(C))%*%rnorm(n)
  #link
  z = y+rnorm(n, sd = sigma)
  gauss_llh = function(y_o, z) sum(-.5*(z-y_o)^2/sigma^2) -n*(log(sigma)+log(2*pi)/2)
  gauss_hess = function(y_o, z) diag(array(rep(1/sigma^2, length(y_o))))
  gauss_score = function(y_o, z) (z-y_o)/sigma^2
  return(list("type"="Gaussian", "locs" =  matrix(locs,ncol=dimen), "z"=z, "y"=y, "C" = C, "hess" = gauss_hess, "score"=gauss_score, "llh"=gauss_llh))
}





################# Gamma  #########################
## one parameter: assume alpha known (alpha = 1: exponential)
## two parameter: what does it mean?
#beta^alpha / Gamma(alpha) * x^{alpha - 1} * exp(-bx)
#exp(-bx + (alpha-1)log(x) +alpha*log(beta)-log(Gamma(alpha)))
gamma_sample = function(n, covfun, seed_val=125, dom=1, dimen=1, alpha = 2 ){
  set.seed(seed_val)
  locs =create_locs(n, dimen, dom)
  C = covfun(locs)
  y_beta = t(chol(C))%*%rnorm(n)

  #link
  z = rgamma(n, shape = alpha, rate = exp(y_beta))

#  d_beta =function(y_o, z) = -z+y_o[,1]/y_o[,2]
#  d_alpha = function(y_o, z) = log(z)+log(y_o[,2])-digamma(y_o[,1])
#  d_bb = function(y_o, z) -y_o[,1]/y_o[,2]^2 # matrix?  cross matrix?
#  d_aa =function(y_o, z) derivative of digamma(y_o[,1])
#  d_ab = function(y_o, z) 1/y_o[,2]

  gamma_hess = function(y_o, z) diag(array(z*exp(y_o)))
  gamma_score = function(y_o, z) -z*exp(y_o)+ alpha
  gamma_llh = function(y_o, z) sum(-y_o*z + (alpha-1)*log(z) +alpha*log(y_o)-n*log(Gamma(alpha)))
  return(list("type"="Gamma", "locs" =  matrix(locs,ncol=dimen), "z"=z, "y"=y_beta, "C" = C, "hess" = gamma_hess, "score"=gamma_score, "alpha"=alpha, "llh" = gamma_llh))
}




################# Negative Binomial (NOT FINISHED) #########################
#ncr(k+r-1,k)*q^r*p^k
#(k successes,r fail)
# not implemented yet
negbin_sample = function(n, covfun, seed = 125, dom = 1){

  set.seed(seed)
  # create locations for a domain [0, dom]
  locs = seq(0,dom,length.out=n)
  C = covfun(locs)
  y = t(chol(C))%*%rnorm(n)
  #link

  eta = exp(y)/(1+exp(y))
  # P(y = 1) = cdf(y)
  # generate observations
  z = rnbinom(n,1,prob = eta)

  # store score and hessian functions, update
  negbin_hess = 0#function(y_o, z) diag(array(exp(y_o)/(1+exp(y_o))^2))
  negbin_score = 0 # function(y_o, z) z - exp(y_o)/(1+exp(y_o))

  # return object with all components of model
  return(list("locs" =  matrix(locs,ncol=1), "z"=z, "y"=y, "C"=C, "hess" = logistic_hess, "score"=logistic_score))

}

##########################################################################################
###################### Testing:  Replace with unit tests? ################################
##########################################################################################
#covfun <- function(locs) Matern(rdist(locs), range = .1, smoothness = 1)
#gms = gamma_sample(400, covfun, dimen=2)
# par(mfrow=c(1,1))
# plot(gms$locs, gms$z)
# plot(gms$locs, gms$y)
# pred_y_lv = estimate_laplace_vec(gms, covfun, m=3, get_sd = TRUE)
# pred_y_l = estimate_laplace(gms)
# points(gms$locs, pred_y_lv$mean, type = "l", col=3)
# points(gms$locs, pred_y_l$mean, type = "l", col=2)
#

# image(gms$locs[1:20,2],gms$locs[1:20,2], matrix(gms$z, nrow = n_locs_per_dim, byrow = TRUE))

#
# ls = logistic_sample(400, covfun, dimen=2, seed = 32 )
# image(ls$locs[1:20,2],ls$locs[1:20,2], matrix(ls$z, nrow = n_locs_per_dim, byrow = TRUE))
# image(ls$locs[1:20,2],ls$locs[1:20,2], matrix(ls$y, nrow = n_locs_per_dim, byrow = TRUE))
#
# par(mfrow=c(1,2))
# test_hd  = estimate_laplace_vec(ls, covfun, m=3, get_sd=TRUE)
# image(ls$locs[1:20,2],ls$locs[1:20,2], matrix(test_hd$mean, nrow = n_locs_per_dim, byrow = TRUE))
# image(ls$locs[1:20,2],ls$locs[1:20,2], matrix(ls$y, nrow = n_locs_per_dim, byrow = TRUE))
#
