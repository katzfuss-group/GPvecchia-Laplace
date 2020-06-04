
source("VL_scripts/importer.R")
library(mvtnorm)
#library(RandomFields)

# factory method to estimate Vecchia or LowRank log likelihood
get_vllh_optim = function(data_gen, vecchia.approx, vecchia.approx.ZY){
  fit_VL <- function(params){
    log_rng=params[1]; log_smth= params[2] 
    rng = exp(log_rng)
    smth = exp(log_smth)
    covparms=c(1, rng, smth)
    pred_lv = calculate_posterior_VL(data_gen$z, vecchia.approx, data_gen$type, covparms, return_all = TRUE)
    pred_lv$vec_lh = vecchia_likelihood(pred_lv$t, vecchia.approx.ZY,covparms,pred_lv$D)
    pseudo_marginal_loglik_vecchia=pred_lv$vec_lh
    post_mean = pred_lv$mean
    pseudo_cond_loglik = dmvnorm(x=t(pred_lv$t), mean = post_mean, sigma = diag(pred_lv$D), log = TRUE)
    true_loglike = data_gen$llh(post_mean, data_gen$z)
    loglik_vecchia = pseudo_marginal_loglik_vecchia + true_loglike - pseudo_cond_loglik
    if (abs(loglik_vecchia)>1e7) loglik_vecchia=sign(loglik_vecchia)*1e7
    -loglik_vecchia
  }
  
  return(fit_VL)
  
}



# Factory method to estimate laplace-only log likelihood
get_lap_llh_optim = function(data_gen){
  fit_laplace <- function(params){
    log_rng=params[1]; log_smth= params[2] 
    rng = exp(log_rng)
    smth = exp(log_smth)
    covparms=c(1, rng, smth)
    covfun <- function(locs) Matern(rdist(locs), range = rng, smoothness = smth)
    pred_l = calculate_posterior_laplace(data_gen$z, data_gen$type, covfun( data_gen$locs), return_all = TRUE )
    pseudo_marginal_loglik = dmvnorm(x=t(array(pred_l$t)),  sigma = as.matrix(pred_l$C+pred_l$D), log = TRUE)# mean = t(array(pred_l$mean)),
    pseudo_cond_loglik = dmvnorm(x=t(array(pred_l$t)), mean = t(array(pred_l$mean)), sigma = as.matrix(pred_l$D), log = TRUE)
    true_loglike = data_gen$llh(pred_l$mean, data_gen$z)
    loglik_laplace = pseudo_marginal_loglik + true_loglike - pseudo_cond_loglik
    -loglik_laplace
  }
  return(fit_laplace)
}

# data generating covariance
covfun <- function(locs) Matern(rdist(locs), range = .05, smoothness = .5)
covmodel <- RMmatern(.5,var = 1,scale = .05)

## The following methods generate data according to a seed
## then estimate the optimal parameters for the respective method

# parameter estimation for Laplace 
run_Laplace_optim = function(seed_value){
  data_gen = pois_sample(25^2, covmodel,seed = seed_value, dimen = 2)
  laplace_llh_function = get_lap_llh_optim(data_gen)
  params = c(-3,-.69)
  Laplace_res = optim(params, laplace_llh_function)#, control = list(trace = 1))
  sprintf("Laplace range: %.3f  smooth: %.3f ", exp(Laplace_res$par)[1],  exp(Laplace_res$par)[2])
  exp(Laplace_res$par)
}

# parameter estimation for LowRank
run_LR_optim = function(seed_value, m, bdd=T){
  data_gen = pois_sample(25^2, covmodel,seed = seed_value, dimen = 2)
  vecchia.approx.lr=vecchia_specify(data_gen$locs, m, conditioning = "firstm")
  lr_llh_function = get_vllh_optim(data_gen, vecchia.approx.lr, vecchia.approx.lr)
  params = c(-3,-.69)
  LR_res = NA
  if(bdd){
    LR_res = optim(params, lr_llh_function, method = "L-BFGS-B", lower = c(-6.9, -6.9), upper=c(log(10), log(10)))
  }
  if(!bdd){
    LR_res = optim(params, lr_llh_function)#, control = list(trace =1))
  }
  message(sprintf("%i LR range: %.3f  smooth: %.3f ", seed_value, exp(LR_res$par)[1],  exp(LR_res$par)[2]))
  exp(LR_res$par)
}

# parameter estimation for VL
run_VL_optim = function(seed_value, m){
  data_gen = pois_sample(25^2, covmodel,seed = seed_value, dimen = 2)
  vecchia.approx=vecchia_specify(data_gen$locs, m, cond.yz = "zy")
  vecchia.approx.ZY=vecchia_specify(data_gen$locs, m)
  vl_llh_function = get_vllh_optim(data_gen, vecchia.approx, vecchia.approx.ZY)
  params = c(-3,-.69)
  VL_res = optim(params, vl_llh_function)#,  control = list(trace = 1))
  message(sprintf("%i VL range: %.3f  smooth: %.3f ", seed_value, exp(VL_res$par)[1],  exp(VL_res$par)[2]))
  exp(VL_res$par)
}

nseeds = 101
seed_set =  sample( 1e6, nseeds, replace= F)

# laplace
param_laplace  = lapply(seed_set, FUN = function(x) try(run_Laplace_optim(x), TRUE))
laplace_results = cbind(seed_set, -1,"Laplace", matrix(unlist(param_laplace), ncol=2, byrow = T) )
save(laplace_results, file = "results_laplace.RData")


## m=5
param_VL_m5  = lapply(seed_set, FUN = function(x) try(run_VL_optim(x, 5), TRUE))
VL_results_m5 = cbind(seed_set, 5,"VL", matrix(unlist(param_VL_m5), ncol=2, byrow = T) )
param_LR_m5 = lapply(seed_set, FUN = function(x) try(run_LR_optim(x, 5), TRUE))
success_seeds_5 = unlist(lapply(1:nseeds, FUN = function(i) length(param_LR_m5[[i]])==2))
LR_results_m5 =  cbind(seed_set[success_seeds_5], 
                       5,  "LR", 
                       matrix(unlist(param_LR_m5[success_seeds_5]), ncol=2, byrow=T))

#save(LR_results_m5,VL_results_m5, file = "results_m5_bdd.RData")


#m=10 
param_VL_m10  = lapply(seed_set, FUN = function(x) try(run_VL_optim(x, 10), TRUE))
VL_results_m10 = cbind(seed_set, 10, "VL", matrix(unlist(param_VL_m10), ncol=2, byrow = T) )
param_LowRank_m10 = lapply(seed_set, FUN = function(x) try(run_LR_optim(x, 10), TRUE))
success_seeds_10 = unlist(lapply(1:nseeds, FUN = function(i) length(param_LowRank_m10[[i]])==2))
LR_results_m10 =  cbind(seed_set[success_seeds_10],
                        10,"LR", 
                        matrix(unlist(param_LowRank_m10[success_seeds_10]), ncol=2, byrow=T))

#save(LR_results_m10,VL_results_m10, file = "results_m10_bdd.RData")

# m=20
param_VL_m20  = lapply(seed_set, FUN = function(x) try(run_VL_optim(x, 20), TRUE))
VL_results_m20 = cbind(seed_set, 20, "VL", matrix(unlist(param_VL_m20), ncol=2, byrow = T) )
param_LowRank_m20 = lapply(seed_set, FUN = function(x) try(run_LR_optim(x, 20), TRUE))
success_seeds_20 = unlist(lapply(1:nseeds, FUN = function(i) length(param_LowRank_m20[[i]])==2))
LR_results_m20 =  cbind(seed_set[success_seeds_20],20, "LR", 
                        matrix(unlist(param_LowRank_m20[success_seeds_20]), ncol=2, byrow=T))

#save(LR_results_m20,VL_results_m20, file = "results_m20_bdd.RData")


# m=40
param_VL_m40  = lapply(seed_set, FUN = function(x) try(run_VL_optim(x, 40), TRUE))
VL_results_m40 = cbind(seed_set, 40, "VL", matrix(unlist(param_VL_m40), ncol=2, byrow = T) )
param_LowRank_m40 = lapply(seed_set, FUN = function(x) try(run_LR_optim(x, 40), TRUE))
success_seeds_40 = unlist(lapply(1:nseeds, FUN = function(i) length(param_LowRank_m40[[i]])==2))
LR_results_m40 =  cbind(seed_set[success_seeds_40],40, "LR", 
                        matrix(unlist(param_LowRank_m40[success_seeds_40]), ncol=2, byrow=T))
#save(LR_results_m40,VL_results_m40, file = "results_m40.RData")

# combined
my_df = data.frame(rbind(
  VL_results_m5, LR_results_m5, 
  VL_results_m10, LR_results_m10, 
  VL_results_m20, LR_results_m20, 
  VL_results_m40, LR_results_m40))
lap_df = data.frame(laplace_results)

names(my_df) = c("seed", "m", "method", "Range", "Smoothness")
names(lap_df) = c("seed", "m", "method", "Range", "Smoothness")

# data files are aggregated and used to make plots elsewhere.  

# in some cases, the LowRank method would fail so we would run
# a bounded version of the optimization for the missing values. 

LR_part = my_df[my_df$method == "LR",]
missing_seeds = data.frame()
for (m in c(5,10,20,40)) {
  for( idx in 1:101) {
    seedum = seed_set[idx]
    if(!seedum %in% LR_part[LR_part$m == m,]$seed){
        missing_seeds = rbind(missing_seeds, c(m, seedum))
    }
  }
}

# then for the missing seeds 
my_df_bdd = apply(missing_seeds, MARGIN=1, FUN = function(x) try(run_LR_optim(x[2], x[1]), TRUE))
success_seeds_miss = unlist(lapply(1:nrow(missing_seeds), FUN = function(i) length(my_df_bdd[[i]])==2))
LR_results_miss =  cbind(missing_seeds[success_seeds_miss, 2], 
                       missing_seeds[success_seeds_miss, 1],  "LR", 
                       matrix(unlist(my_df_bdd[success_seeds_miss]), ncol=2, byrow=T))

