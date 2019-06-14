# This code evaluates the likelihood for various models and methods on
# a grid of parameter values.

source("VL_scripts/importer.R") #library(GPVecchia)
require(ggplot2)


#####################################################################
######################    Parameter Estimation  #####################
#####################################################################

# estimate Vecchia log likelihood
get_vllh = function(data_gen, pred_lv, ms = 1, ps = 0){
  pseudo_marginal_loglik_vecchia=pred_lv$vec_lh
  # for testing approximation, allow perturbation of mean
  #n=length(pred_lv$mean)
  #perturbed = rnorm(n, sd = ps)
  post_mean = pred_lv$mean#*ms + perturbed
  pseudo_cond_loglik = dmvnorm(x=t(pred_lv$t), mean = post_mean, sigma = diag(pred_lv$D), log = TRUE)
  true_loglike = data_gen$llh(post_mean, data_gen$z)
  loglik_vecchia = pseudo_marginal_loglik_vecchia + true_loglike - pseudo_cond_loglik
  return(loglik_vecchia)
}

# Estimate laplace-only llh
get_lap_llh = function(data_gen, pred_l, ms = 1){
  pseudo_marginal_loglik = dmvnorm(x=t(array(pred_l$t)),  sigma = as.matrix(pred_l$C+pred_l$D), log = TRUE)# mean = t(array(pred_l$mean)),
  pseudo_cond_loglik = dmvnorm(x=t(array(pred_l$t)), mean = t(array(pred_l$mean))*ms, sigma = as.matrix(pred_l$D), log = TRUE)
  true_loglike = data_gen$llh(pred_l$mean*ms, data_gen$z)
  loglik_laplace = pseudo_marginal_loglik + true_loglike - pseudo_cond_loglik
  return(loglik_laplace)
}

# calculate the exact llh, for Gaussian case
gauss_true_ll=function(z, K){
  require(mvtnorm)
  #assumes likelihood model had sigma = 1
  dmvnorm(t(array(z)), sigma = K, log = TRUE)
}

# estimat the exact llh, for non-gaussian case?
numerical_true_ll = function(data_gen, covparms){
  #fix seed
  #set.seed(110101)
  # alternately:  fix iid N(0,1) input and multiply by Cholesky

  # for small sample, numerically calculate likelihood int z|y y dy
  print("Sampling MC true llh")
  llh = data_gen$llh
  z = data_gen$z
  locs = data_gen$locs
  covfun <- function(locs) Matern(rdist(locs), range = covparms[2], smoothness = covparms[3])
  #  instead:  use MC
  sig = covfun(locs)
  L = t(chol(sig))
  mc = function(y) exp(llh(y,z))
  L%*%rnorm(length(z))
  marg_x = mean(apply(rmvnorm(n=100, sigma = sig), 1, mc))
  return(log(marg_x))
}

nbrs = 1
# Loop through grid, estimating llh for each parameter point and model
loop_parameters = function(smooth_param_seq, range_param_seq, nbrs, data_gen){
  require(mvtnorm)
  llh_vals = c()
  #Lap_llh_vals = c()
  #exact_llh_vals = c()
  #options(show.error.messages = FALSE)
  for (smth in smooth_param_seq){
    y_init = NA
    y_init_LR=  NA
    for (rng in range_param_seq){
      covparms=c(1,rng,smth)
      covfun <- function(locs) Matern(rdist(locs), range = rng, smoothness = smth)
      pred_l = calculate_posterior_laplace(data_gen$z, data_gen$type, covfun(data_gen$locs), return_all = TRUE )

      # Laplace likelihood, no vecchia
      laplace_llh=NA
      try(( laplace_llh = get_lap_llh(data_gen, pred_l)))

      # MC based marginal
      num_llh = NA#numerical_true_ll(data_gen, covparms)


      # for vecchia, evaluate all required neighbor counts
      vecc_vals=c()
      lr_vals = c()
      for (m in nbrs){
        # approximate the latent y
        use_zy = TRUE
        use_lr = TRUE
        vecchia.approx=vecchia_specify(data_gen$locs, m, cond.yz = "zy")
        pred_lv = calculate_posterior_VL(data_gen$z, vecchia.approx, data_gen$type, covparms, return_all = TRUE, y_init = y_init)
        y_init = pred_lv$mean
        vecchia_llh=NA
        if(use_zy){
          # switch to SGV for prediction
          vecchia.approx.ZY=vecchia_specify(data_gen$locs, m)
          pred_lv$vec_lh = vecchia_likelihood(pred_lv$t, vecchia.approx.ZY,covparms,pred_lv$D)
        }
        try((vecchia_llh = get_vllh(data_gen, pred_lv, ms=1)))
        vecc_vals = c(vecc_vals,vecchia_llh)

        # Evlauate low rank, if needed
        if(use_lr){
          vecchia.approx.lr=vecchia_specify(data_gen$locs, m, conditioning = "firstm")
          pred_lr = calculate_posterior_VL(data_gen$z,vecchia.approx.lr, data_gen$type, covparms, return_all = TRUE, y_init = y_init_LR)
          y_init_LR = pred_lr$mean
          pred_lr$vec_lh = vecchia_likelihood(pred_lr$t, vecchia.approx = vecchia.approx.lr , covparms, pred_lr$D)
          vecchia_lr_llh = get_vllh(data_gen, pred_lr, ms=1)
          lr_vals = c(lr_vals, vecchia_lr_llh)
        }

      }


      # for gaussian, compare to truth:
      K=as.matrix(covfun(data_gen$locs)+pred_lv$D)
      exact_gauss_llh = NA#gauss_true_ll(pred_lv$t, K)

      # store results
      param_results = c(smth, rng,laplace_llh, num_llh, exact_gauss_llh, vecc_vals, lr_vals)
      param_results
      llh_vals = rbind(llh_vals, param_results)
    }
  }
  #options(show.error.messages = TRUE)
  return(llh_vals)
}
####################    param_estimation  ####################

covfun <- function(locs) Matern(rdist(locs), range = .05, smoothness = .5)
data_gen = pois_sample(25^2, covfun,seed = 12421, dimen = 2)#sample(1e5,1)
smooth_param_seq =seq(.05, 2, .05)
range_param_seq = seq(.01, .75, .04)
m = c(5,10,20)

t_start = Sys.time()
x = loop_parameters(smooth_param_seq, range_param_seq , m , data_gen)
t_end = Sys.time()

colnames(x) <- c("Smooth", "Range", "Laplace", "MC", "Exact", "Vecchia-5",  "Vecchia-10",  "Vecchia-20", "LR-5",  "LR-10",  "LR-20")
param_est_data = list("data" = x, "likelihood" = "poisson", "seed"=12421, "dimen"=2, "n" = 625)
save(param_est_data, file = "VL_scripts/saved_data/param_est_pois_VLZY_LR_yinit_2.Rdata")
write.csv(x = x,file = "VL_scripts/saved_data/param_est_pois_2D_VLZY_LR_yinit_2.csv", row.names = FALSE)

