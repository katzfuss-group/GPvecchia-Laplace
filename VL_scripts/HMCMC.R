
#**********************************************************************************
###################   Compare HMC to VL           ################################
#*********************************************************************************


# run HMC for relatively small sample
# run VL or laplace for comparable time


#####################  data generation  #######################
source("VL_scripts/importer.R")
dimen=2 # number of spatial dimensions
samp_size=625  # number of observed locs
domn = 1# domain over which to generate points

sig2=1; range=.05; smooth=.5
covparms=c(sig2,range,smooth)

covfun <- function(locs) sig2*Matern(rdist(locs),range=range,smoothness=smooth)

# current data generation is over a uniform lattice
#source(file.path("VL_scripts/data_generating_functions.R"))
ls=logistic_sample(samp_size,covfun, seed = 12, dom = domn, dimen=dimen)

m=5
vecchia.approx=vecchia_specify(ls$locs, m, cond.yz='zy')#, cond.yz = "z"  )
vecchia.approx_LL=vecchia_specify(ls$locs, m, conditioning = "firstm")#, cond.yz = "z"  )
#####################  comparison  #######################
mse_fun = function(x) (mean((x - ls$y)^2))
crps_fun = function(x) scoringRules::crps_sample(ls$y[,1], x)
time_st = Sys.time()
# for quick evaluation, compare posteriors
ll_pred = calculate_posterior_VL(ls$z, vecchia.approx_LL, ls$type, covparms)
vl_pred = calculate_posterior_VL(ls$z, vecchia.approx, ls$type, covparms)
lap_pred = calculate_posterior_laplace(ls$z, ls$type, covfun(ls$locs), return_all = TRUE)

#### Timings and accuracy for VL and Laplace methods ####
# To give smooth curve for VL, average over different orderings from vecchia specify
mse_fun(vl_pred$mean);mse_fun(lap_pred$mean)
mini_dt = c(); mini_dt_crps = c()
set.seed(15)
for(m in c(1,5,10,20, 40)){#
  run_avg = c()
  for (L in 1:30){
    vecchia.approx=vecchia_specify(ls$locs, m, cond.yz = "zy")
    time_st = Sys.time()
    vl_pred = calculate_posterior_VL(ls$z,vecchia.approx, ls$type, covparms, return_all = T)
    time_end = Sys.time()
    run_avg = rbind(run_avg,c(as.double(difftime(time_end, time_st, units = "secs")), mse_fun(vl_pred$mean)))
    
    sd_VL = sqrt(diag(solve(vl_pred$W)))
    crps_avg= rbind(crps_avg,c(as.double(difftime(time_end, time_st, units = "secs")), 
                              mean(scoringRules::crps_norm(y = ls$y, mean =vl_pred$mean, sd = sd_VL))))
  }
  res = apply(X = run_avg,MARGIN = 2,FUN = median)
  mini_dt = rbind(mini_dt, c(m, res))
  res_crps = apply(X = crps_avg,MARGIN = 2,FUN = median)
  mini_dt_crps = rbind(mini_dt_crps, c(m, res_crps))
}
time_st = Sys.time()
lap_pred = calculate_posterior_laplace(ls$z, ls$type, covfun(ls$locs), return_all = TRUE)
time_end = Sys.time()
lap_tim = as.double(difftime(time_end, time_st, units = "secs"))

# Rdata files generated in HMC section, time consuming
load("VL_scripts/saved_data/hmc_big_summary.RData")
load("VL_scripts/saved_data/hmc_big_path.RData")
time_inc = hmc_big_res$runtime_seconds/length(hmc_large_means_time$HMC_mse_by_sample)


#calculate improvements in HMC
chunks = (1000000-10000)/100  # thin every 10, chunk = 10 thinned samples
mse_by_sample = c()
for(i in 1:chunks){
  mode_approx = colMeans(hmc_big_path$thinned[1:(10*i),])
  mse_by_sample = c(mse_by_sample, mse_fun(mode_approx))
}
time_inc = hmc_big_res$runtime_seconds/chunks

#normalize:
laplace_mse = mse_fun(lap_pred$mean)
mse_by_sample = sqrt(mse_by_sample/laplace_mse)
mini_dt[,3]=sqrt(mini_dt[,3]/laplace_mse)
laplace_mse= 1


#  Show time and accuracy comparison
pdf("VL_scripts/plots/HMC_over_time_big.pdf", width = 9*3/4, height = 5*3/4)
plot(1:chunks*time_inc, mse_by_sample, log = "x", xlim = c(0.01,30000),
     type = "l", col=4, xlab = "Time (seconds)", ylab = "RRMSE")
abline(h=1, col="grey70")
points(mini_dt[, 2], mini_dt[, 3], type  ="l", col=3)
#points(c(mini_dt[4, 2], 1e4), rep(mini_dt[4, 3],2), type= "l", lty=2, lwd=1, col=3)
points(c(mini_dt[, 2]), mini_dt[, 3], pch=18, col=3)
points(lap_tim,laplace_mse)
points(1:chunks*time_inc, mse_by_sample, type="l", col=4)
#points(c(.05, 1e4), c(laplace_mse,laplace_mse), type= "l", lty=3, col=11)
legend("topleft", legend=c("HMC", "VL", "Laplace"), col = c(4,3,1), lty = c(1,2,NA), pch=c(NA, 18, 1))
dev.off()



# Conduct same analysis as MSE with CRPS
crps_by_sample = c()
for(i in 1:chunks){
  crps_chunk = mean(crps_fun(t(hmc_big_path$thinned[1:(10*i),])))
  crps_by_sample = c(crps_by_sample, crps_chunk)
}
time_inc = hmc_big_res$runtime_seconds/chunks
#normalize to laplace
laplace_crps = mean(scoringRules::crps_norm(ls$y, lap_pred$mean[,1], lap_pred$sd))
crps_by_sample_normed = sqrt(crps_by_sample/laplace_crps)
mini_dt_crps[,3]=sqrt(mini_dt_crps[,3]/laplace_crps)
laplace_crps= 1


#  Show time and accuracy comparison
pdf("HMC_over_time_CRPS_E6.pdf", width = 9*3/4, height = 5*3/4)
plot(1:chunks*time_inc, t(crps_by_sample_normed), log = "x", xlim = c(0.01,30000),
     type = "l", col=4, xlab = "Time (seconds)", ylab = "CRPS")
abline(h=1, col="grey70")
points(mini_dt_crps[, 2], mini_dt_crps[, 3], type  ="l", col=3)
points(c(mini_dt_crps[, 2]), mini_dt_crps[, 3], pch=18, col=3)
points(lap_tim,laplace_crps)
points(1:chunks*time_inc, crps_by_sample_normed, type="l", col=4)
legend("topleft", legend=c("HMC", "VL", "Laplace"), col = c(4,3,1), lty = c(1,2,NA), pch=c(NA, 18, 1))
dev.off()



### HMC GG trace plots
par(mfrow=c(4,1))
thinned_path = data.frame(cbind(1:nrow(hmc_big_path$thinned),hmc_big_path$thinned[,c(10,45,245,445)]))
colnames(thinned_path)<- c("Sample", "10","45","245","445")
melty_thin = melt(thinned_path, id = "Sample")
colnames(melty_thin)[3] <- "Value"
ggplot(data = melty_thin, aes(y = Value, x=Sample))+
  geom_point(pch=".")+theme_bw()+
  facet_wrap(.~ variable, scales = "free_y", ncol=1)
ggsave("VL_scripts/plots/HMC_trace_GG.pdf", device= "pdf",width = 8, height = 4)

acf(hmc_big_path$thinned[,10])



##### Further evaluaton for HMC ####
load("VL_scripts/saved_data/hmc_big.RData")

#### Trace and Auto correlation data ####
pdf("VL_scripts/plots/HMC_trace.pdf", height = 3, width = 7)
par(mfrow=c(3,1))
for (i in c(88, 175, 225)) {
  plot(hmc_pred_big$burned[,i][seq(1,200000-5000,500)],
       ylab = NA,
       main = paste("HMC Trace, index", i), pch=".")
}

acf(hmc_pred_big$burned[,145][seq(1,200000-5000,500)])
count = 0
for(i in 2:195001){
  diff_val = any(hmc_pred_big$burned[i,] != hmc_pred_big$burned[i-1,])
  if(!diff_val) count =count + 1
}
dev.off()


pdf("VL_scripts/plots/HMC_over_time.pdf", width = 6, height = 4)
par(mfrow=c(1,1), mar =c(3,3,1,1))
plot(1:39*2, mse_by_sample, type = "l", ylim = c(.64, .9),
     ylab = "MSE", xlab = "Minutes")
abline(h = mse_fun(lap_pred$mean), col=1, lwd =2)
abline(h = mse_fun(vl_pred$mean), col=3, lwd = 2)
points(.1, mse_fun(lap_pred$mean), col =1)
points(.1, mse_fun(vl_pred$mean), col =3)
points(1:39*2, mse_by_sample, pch = 1, ylim = c(.64, .9))
legend("topright", legend = c("HMC", "Laplace", "VL"), lwd = c(1,2,2), col=c(1,1,3))
dev.off()

#*****************************************************************************************
###################    Hamiltonian MCMC:   Implementation      ##########################
#*****************************************************************************************
library(fields)



## Leapfrog steps for generating the next point
gen_HMC_sample = function(U, grad_U, epsilon, lf_steps, q_current){
  # q is the current state / location
  q = q_current
  # generate random p (velocity) to initialize
  m=1
  p = rnorm(length(q), 0, m)
  p_current  = p

  # Leapfrog 1: take initial half step of velocity given random start velocity vector p
  p = p-epsilon*grad_U(q)/2

  # Leapfrog 2:  alternate full steps between momentum (random steps) and potential (likelihood)
  for(i in 1:lf_steps){
    q = q + epsilon * p
    # do not take a full step at the end
    if (i != lf_steps) p = p-epsilon*grad_U(q)
  }
  # Leapfrog 3: take final half step of momentum given last momentum p
  p = p-epsilon*grad_U(q)/2

  # negate momentum; given KE is even, doesnt affect momentum but makes proposal symmetric (ie reverse path)
  p=-p
  if (length(q_current)==1) print("ERRROR")
  #prepare to accept or reject
  U_prop = U(q)
  K_prop = sum(p^2) / (2*m^2)  # gaussian (0,1)
  U_old = U(q_current)
  K_old = sum(p_current^2) / (2*m^2)

  if(runif(1)< exp(U_old + K_old - U_prop - K_prop))
    return(q)
  return(q_current)
}


# Metropolis-Hastings loop
run_HMC<- function(lik_mod, niter =100000, burnin=5000 ){

  nlocs = length(lik_mod$z)
  C_inv = solve(lik_mod$C)
  # full posterior up to normalizing constant
  U_potential = function(y) -(lik_mod$llh(y, lik_mod$z)-1/2*t(as.vector(y))%*%C_inv%*%as.vector(y))
  grad_U_potential = function(y) -(lik_mod$score(y, lik_mod$z)-C_inv%*%as.vector(y))

  samples_run = niter
  path = matrix(nrow = samples_run, ncol = nlocs)
  y_prev = rep(0,nlocs)
  path[1,]= y_prev
  epsilon = .001
  lf_steps = 50

  count_accepts = 0
  for (i in 2:samples_run){
    y_prev =  gen_HMC_sample(U_potential, grad_U_potential, epsilon, lf_steps, y_prev)
    if (any(y_prev!= path[i-1,])) count_accepts= count_accepts+1
    path[i,]<- y_prev
  }

  count_accepts/samples_run

  burn = path[burnin:samples_run,]
  hthin = burn[seq(1,samples_run-burnin,10),]
  y_hmcmc = colMeans(hthin)
  y_hvar = sqrt(apply(hthin, MARGIN = 2, var))
  return(list("y_HMC"=y_hmcmc, "sd_HMC"=y_hvar, "burned"=burn, "thinned" = hthin))
}

lkhd_seed = 12
#source("VL_scripts/data_generating_functions.R")
ls=logistic_sample(samp_size,covfun, seed = lkhd_seed, dom = domn, dimen=dimen)
mse_fun = function(x) mean((x - ls$y)^2)

set.seed(1234099)
t_start_big = Sys.time()
hmc_pred_big = run_HMC(ls, niter = 300000, burnin = 10000)
t_end_big = Sys.time()
HMC_time_big = as.double(difftime(t_end_big, t_start_big, units = "secs"))


hmc_big_res = list("runtime_seconds"=HMC_time_big, "hmc_result" = hmc_pred_big$y_HMC,
                   "sample_dim_dom" = c(samp_size, dimen, domn), "covparms"=covparms, "hmc_seed"=1234099,
                   "lklhd" = "logistic", "lkhd_seed" =lkhd_seed, "mse" = mse_fun(hmc_pred_big$y_HMC))
save(hmc_big_res, file = "saved_data/hmc_big_summary.RData")

hmc_big_path = list("thinned" = hmc_pred_big$thinned, "lklhd" = "logistic", "lkhd_seed" =lkhd_seed,
                    "hmc_seed"=1234099, "accept_rate"=hmc_pred_big$accept)
save(hmc_big_path, file = "saved_data/hmc_big_path.RData")


##### Simulation for 20 samples of 8K HMC ####
HMC_CRPS_samps = rep(0,20)
for(i in 2:20){
  hmc_pred_8k = run_HMC(ls, niter = 8000, burnin = 5000)
  HMC_CRPS_samps[i] = mean(crps_fun(t(hmc_pred_8k$thinned)))
}


#***************************************************************************************************
###################    MCMC (No HMC)           ##################################################
#***************************************************************************************************

# For comparison to HMC
library(mvtnorm)

samples_run = 100000
path = matrix(nrow = samples_run, ncol = nlocs)
y_prev = rep(0,nlocs)
cond_z = function(z,y) prod(exp(-exp(y))*exp(y)^z)
path[1,] = y_prev
reject_alpha_hist = rep(0,samples_run)
for (i in 2:samples_run){
  y_prop = rnorm(nlocs, mean = y_prev, sd = .07)
  #alpha = cond_z(z1,y_prop)*cond_z(z2,y_prop)*dmvnorm(y_prop,mean = rep(0,10),C)/
  # (cond_z(z1,y_prev)*cond_z(z2,y_prev)*dmvnorm(y_prev,mean = rep(0,10),C))
  alpha = cond_z(z,y_prop)*dmvnorm(y_prop,mean = rep(0,nlocs),C) /
    (cond_z(z,y_prev)*dmvnorm(y_prev,mean = rep(0,nlocs),C))
  reject_alpha = runif(1)>alpha
  y_prev <- if(reject_alpha) y_prev else y_prop
  reject_alpha_hist[i]=reject_alpha
  path[i,]<- y_prev
}
mean(reject_alpha_hist)
## After running, remove burn in and thin samples
burn = path[5000:samples_run,]
thin = burn[seq(1,samples_run-5000,100),]

# estimate pointwise mean and variance
y_mcmc = colMeans(thin)
y_var = sqrt(apply(thin, MARGIN = 2, var))

#plot
pdf("MCMC_mix.pdf")
par(mfrow=c(2,2))
plot(thin[,1])
plot(thin[,12])
plot(thin[,30])
plot(thin[,45])
dev.off()


pdf("MCMC_f.pdf")
plot(locs,y, type = "p", ylim=c(-.5, 2), main = "Posterior Estimate, f(s)|y", lwd = 3)
points(locs,y_mcmc, type = "l", col=4, lwd = 3)
points(locs,y_mcmc+1.65*y_var, type = "l", pch = 2,col=4, lty = 2, lwd = 3)
points(locs,y_mcmc-1.65*y_var, type = "l", pch = 2, col=4, lty = 2, lwd = 3)
dev.off()

# explore
cov(burn)
mean(reject_alpha_hist)
plot(thin[,1])
f = hist(thin[,1], breaks = 50)
abline(v = y[1], col=2)
f$mids[which(f$counts==max(f$counts))]



#***************************************************************************************************
####################   Example plots, 1D data     ##################################################
#***************************************************************************************************
source("VL_scripts/data_generating_functions.R")

dimen=1 # number of spatial dimensions
domn=1
samp_size = 250
sig2=1; range=.05; smooth=1.5
covparms=c(sig2,range,smooth)

covfun <- function(locs) sig2*Matern(rdist(locs),range=range,smoothness=smooth)

# current data generation is over a uniform lattice
#source(file.path("VL_scripts/data_generating_functions.R"))
ls=logistic_sample(samp_size,covfun, seed = 126, dom = domn, dimen=dimen)

m=5
vecchia.approx=vecchia_specify(ls$locs, m)#, cond.yz = "z"  )
vecchia.approx_LL=vecchia_specify(ls$locs, m, conditioning = "firstm")#, cond.yz = "z"  )
ll_pred = calculate_posterior_VL(ls$z, vecchia.approx_LL, ls$type, covparms)
vl_pred = calculate_posterior_VL(ls$z, vecchia.approx, ls$type, covparms)
lap_pred = calculate_posterior_laplace(ls$z, ls$type, covfun(ls$locs), return_all = TRUE)


hmc_pred_big = run_HMC(ls)
hmc_fit_small1 = run_HMC(ls, niter = 6000)
hmc_fit_small2 = run_HMC(ls, niter = 6000)
hmc_fit_small3 = run_HMC(ls, niter = 6000)
realizations_list = cbind(hmc_fit_small1$y_HMC, hmc_fit_small2$y_HMC, hmc_fit_small3$y_HMC)

pdf("VL_scripts/plots/Simple_comparison_HMC.pdf", height = 7, width = 7)
par(mfrow=c(1,1))#,mar=c(2,2.1,2,1
plot(ls$locs, ls$y, main = "Posterior Predictions, Logistic Obs", type = "p", xlab = "Location", ylab = "Latent Value", pch = ".")
points(ls$locs, lap_pred$mean, type = "l",col=4, lwd = 1)
for(i in 1:3) points(ls$locs, realizations_list[,i], col=6, type = "l")
points(ls$locs, hmc_pred_big$y_HMC, col=5, type = "l")
legend("bottomright", legend = c("Truth", "Laplace", "HMC 6k", "HMC 100k"),
       col = c(1,4,6,5), lty = c(NA,1,1,1), pch = c(".",NA, NA, NA),bty = "y")
dev.off()


pdf("VL_scripts/plots/Simple_comparison.pdf", height = 7, width = 7)
par(mfrow=c(1,1))#,mar=c(2,2.1,2,1)
plot(ls$locs, ls$y, main = "Posterior Predictions, Logistic Obs", type = "p", xlab = "Location", ylab = "Latent Value", pch = ".")
points(ls$locs, vl_pred$mean, type = "l",col=3, lwd = 3)
points(ls$locs, ll_pred$mean, type = "l",col=2, lwd = 1)
points(ls$locs, lap_pred$mean, type = "l",col=4, lwd = 1)
points(ls$locs, hmc_pred_big$y_HMC, col=5, type = "l", pch = 1)
legend("bottom", legend = c("Truth", "Laplace", "LowRank (m=5)", "VL (m=5)", "HMC 100k"),
       col = c(1, 4, 2,3,5), pch = c(".",NA, NA, NA, NA), lty = c(NA, 1,1,1,1), bty = "y")
dev.off()
