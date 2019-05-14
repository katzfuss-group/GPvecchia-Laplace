source("server/importer.R")

load("server/MODIS_analysis/water_vapor_20190328.RData")
locs = as.matrix(locs)
n = length(locs[,1])

global_m = 20 #5 10 20 40
n_sample = 2.5e3

# test on subset
set.seed(123)
sub_idx = sample(length(z), n_sample, replace = FALSE)
pred_idx = sample(setdiff(1:length(z), sub_idx), n_sample, replace = FALSE)

sub_locs = locs[sub_idx,]
z_sub = z[sub_idx]
quilt.plot(sub_locs[,1], sub_locs[,2], z_sub, nx=128, ny =128)



#### trend estimation ####
X = sub_locs
X[,1]= 1



# glm: IRLS method from taking deriv of llh
#IRLS
# t_start = Sys.time()
# beta = c(1,0.001)
# for(i in 1:10){
#   XB = X%*% beta
#   W = -diag.spam(array(exp(-XB)*z_sub))
#   A  = exp(-XB)*z_sub-1
#   U  = W%*% XB - A
#   beta = solve( t(X) %*% W %*% X , t(X) %*% U)
# }
# t_end = Sys.time()
# time_dur = as.double(difftime(t_end, t_start, units = "mins")); print(time_dur)
# print(beta)

# use beta from n=1e6 run
beta = matrix(c(-1.515079549, 0.000765648))
XB = X%*% beta


### Covparm and shape estimation ####
# shape estimation: maximize conditional likelihood.  Integrated likelihood diverges
update_a = function(a_init, covparms,  vecchia.approx, vecchia.approx.IW, XB ){
  a_prev =a_init
  for(i in 1:10){
    t_start = Sys.time()
    posterior = calculate_posterior_VL(z_sub,
                                       vecchia.approx,
                                       likelihood_model="gamma",
                                       covparms=covparms,
                                       likparms = list("alpha"=a_prev),
                                       prior_mean = XB)
    mu = posterior$mean + XB
    llh = function(a) -sum(-a*exp(-mu)*z_sub + (a-1)*log(z_sub)+a*log(a)-a*mu-log(gamma(a))) # concave in a and XB
    param_est = optim(a_prev, llh, method = "Brent", lower = .01, upper = 1e2)
    a = param_est$par; print(a)
    if(abs(a-a_prev) < 1e-5) {print("convergence criteria met (fitting shape parameter)"); break}
    a_prev = a
  }
  return(a)
}

## Covar param estimation:  integrated likelihood
fit_covparms = function(a, covparms_init, vecchia.approx, vecchia.approx.IW, XB){
  vl_likelihood = function(x0){
    theta = exp(x0)
    covparms=c(theta[1], theta[2], theta[3]) # sigma range smoothness
    #sprintf("Evaluating covparms = (%.4f %.4f %.4f)", covparms[1], covparms[2],covparms[3])
    default_lh_params = list("alpha"=a, "sigma"=sqrt(.1))
    # Perform inference on latent mean with Vecchia Laplace approximation
    vll = vecchia_laplace_likelihood(z_sub,
                                     vecchia.approx,
                                     likelihood_model="gamma",
                                     covparms = covparms,
                                     return_all = FALSE,
                                     likparms = default_lh_params,
                                     prior_mean = XB,
                                     vecchia.approx.IW=vecchia.approx.IW,
                                     y_init = NA )
    sprintf("Likelihood for covparms = (%.4f %.4f %.4f): %.4f",
            covparms[1], covparms[2],covparms[3], vll)
    return(-vll)
  }
  x0 = log(covparms_init)
  vl_likelihood(x0)
  res = optim(x0, vl_likelihood, method = "Nelder-Mead", control = list("trace" = 4, "maxit" = 500, "reltol" = 1e-5))
  print(res$convergence); print(exp(res$par))
  return(exp(res$par))
}




## Do parameter estimation for multiple m values
if(FALSE){
  #param_table = matrix(0, nrow =length(m_vals), ncol =5)
  m_rep = global_m
  print(paste("Estimating parameters for m=", m_rep))

  print("Step 1, generating vecchia approximations")
  vecchia.approx = vecchia_specify(sub_locs, m=m_rep, cond.yz = "zy")
  vecchia.approx.IW = vecchia_specify(sub_locs, m=m_rep)

  ## Iterative method:  estimate a, then covparms, then a again
  print("Step 2, optimizing parameters")
  a_prev = 0.9
  covparms_prev = c(.8, 120,   1.3)
  print(a_prev);print(covparms_prev)
  t_start = Sys.time()
  iter_count = 1
  for(i in 1:5){
    a = update_a(a_prev, covparms_prev, vecchia.approx, vecchia.approx.IW, XB)
    covparms = fit_covparms(a, covparms_prev, vecchia.approx, vecchia.approx.IW, XB)
    print(paste("Found shape parameter a=", a))
    print(paste("Found covariance parameters (sig, range, smooth) = ",covparms))

    if(abs(a-a_prev)<1e-3 &
       abs(covparms_prev[1] - covparms[1]) < 1e-2 &
       abs(covparms_prev[2] - covparms[2])/covparms_prev[2] < 1e-2 &
       abs(covparms_prev[3] - covparms[3]) < 1e-2){
      print("Convergence criteria met (fitting all parameters)")
      iter_count = i
      break
    }
    a_prev = a
    covparms_prev = covparms
  }
  t_end = Sys.time()
  time_dur = as.double(difftime(t_end, t_start, units = "mins"))
  print(time_dur)
  cat(a, covparms, m_rep, time_dur, iter_count)
  #param_table[m_idx,] = c(a,covparms, iter_count)



  #param_df = data.frame(cbind(m_vals,param_table))
  #colnames(param_df) <- c("m", "a", "sig","rho","nu", "iters" )
  found_param = c(a, covparms, m_rep, time_dur, iter_count)
  save(found_param, file = "server/MODIS_analysis/saved_data/covparms_m40.RData")
}





#### Posterior estimation and predictions  ####
# set parameters found in exploration
a = 0.88
covparms = c(0.9, 110,   1.5)
global_m_lr = round(global_m^(3/2))
default_lh_params = list("alpha"=a, "sigma"=sqrt(.1))

## Calculate posterior using VL
t_start = Sys.time()  #TIME VL
vecchia.approx = vecchia_specify(sub_locs, m=global_m, cond.yz = "zy")
post_zy = calculate_posterior_VL(z_sub, vecchia.approx, "gamma" , covparms, likparms = default_lh_params, prior_mean = XB)
t_end = Sys.time();  time_dur = as.double(difftime(t_end, t_start, units = "mins"))
paste("Finished VL posterior in", time_dur, "minutes")

## Calculate posterior using LR
t_start = Sys.time() # Time LR
vecchia.approx.lr = vecchia_specify(sub_locs, m=global_m_lr, conditioning = "firstm")
post_lr = calculate_posterior_VL(z_sub, vecchia.approx.lr, "gamma" , covparms, likparms = default_lh_params, prior_mean = XB)
t_end = Sys.time();  time_dur = as.double(difftime(t_end, t_start, units = "mins"))
paste("Finished LR posterior in", time_dur, "minutes")


## Create and save comparison plot
pdf("server/MODIS_analysis/saved_data/MODIS_Compare_post_VL_LR.pdf", width = 12, height = 5)
par(mfrow=c(1,3))
quilt.plot(sub_locs, z_sub, zlim=c(.01, 2.9), main = "Obs", nx = 256, ny=256)
quilt.plot(sub_locs, exp(post_zy$mean), zlim=c(.01,2.9), main = "Posterior exp(VL)", nx = 256, ny=256)
quilt.plot(sub_locs, exp(post_lr$mean), zlim=c(.01,2.9),main = "Posterior exp(LR)", nx = 256, ny=256)
dev.off()

# save objects for future plotting. May be >100MB
save(sub_locs, z_sub, post_zy, post_lr, time_dur,
     file = "server/MODIS_analysis/saved_data/MODIS_posteriors.RData")
paste("Finished posteriors in", time_dur, "minutes")


## Make predictions on a zoomed region, to see small scale structure
print("Running predictions on zoomed observations!"); t_start = Sys.time()

# Subset locations
zoom_idx = intersect(which(locs[,1] >1000 & locs[,1]<1100),  which(locs[,2] <100 & locs[,2] > 0))
locs_zoom = locs[zoom_idx,]
z_zoom =z[zoom_idx]
#quilt.plot(locs_zoom, z_zoom, nx=100, ny=100)

# Calculate trend over prediction region
trend_zoom = locs_zoom
trend_zoom[,1]=1
XB_zoom = trend_zoom%*%beta

# Predict on zoomed region using VL
t_start = Sys.time() # Time VL zoom
z_VLpseudo = post_zy$t - XB
nuggets_VLpseudo = post_zy$D
vecchia.approx.zoom = vecchia_specify(sub_locs, m=global_m, locs.pred=locs_zoom)
preds_zoom=vecchia_prediction(z_VLpseudo, vecchia.approx.zoom, covparms, nuggets_VLpseudo)
t_end = Sys.time();  time_dur = as.double(difftime(t_end, t_start, units = "mins"))
paste("Finished ZY zoom pred in", time_dur, "minutes")

# predictions on zoomed region, Lowrank
t_start = Sys.time() # Time LR
z_VLRpseudo = post_lr$t-XB
nuggets_VLRpseudo = post_lr$D
vecchia.approx.pred.zoom.lr = vecchia_specify(sub_locs, m=global_m_lr, locs.pred=locs_zoom, conditioning = "firstm")
predsLR_zoom=vecchia_prediction(z_VLRpseudo, vecchia.approx.pred.zoom.lr, covparms, nuggets_VLRpseudo)
t_end = Sys.time();  time_dur = as.double(difftime(t_end, t_start, units = "mins"))
paste("Finished LR zoom pred in", time_dur, "minutes")



## Check truth at zoomed locations
zoom_idx_vis = intersect(which(sub_locs[,1] >1000 & sub_locs[,1]<1100),  which(sub_locs[,2] <100 & sub_locs[,2] > 0))
sub_locs_zoom = as.matrix(sub_locs[zoom_idx_vis,])
z_sub_zoom = as.matrix(z_sub[zoom_idx_vis])


## Zoomed plots
pdf("server/MODIS_analysis/paper_plots/MODIS_compare_zoom_VL_LR.pdf", width = 9, height = 8)
par(mfrow=c(2,2), mar = rep(2,4))
quilt.plot(sub_locs_zoom, z_sub_zoom, zlim=c(.01, .9), main = "Observed Data", nx = 100, ny=100,
           add.legend = FALSE, xaxt='n', yaxt = "n")
quilt.plot(locs_zoom, z_zoom, zlim=c(.01, .9), main = "Full Data", nx = 100, ny=100)
par(mar = c(2,1.5,2,1.5))
quilt.plot(offset_locs, exp(predsLR_grid$mu.pred+XB_pred), zlim=c(.01,.9),main = "Predicted, LowRank", nx = 128, ny=128,
           add.legend = FALSE,xaxt='n', yaxt = "n")
quilt.plot(offset_locs, exp(preds_grid$mu.pred+XB_pred), zlim=c(.01,.9), main = "Predicted, VL", nx = 128, ny=128,
           add.legend = FALSE, xaxt='n', yaxt = "n")
dev.off()

# save prediction objects for future plotting
save(locs_zoom, z_zoom, preds_zoom, predsLR_zoom, XB_zoom, file = "server/MODIS_analysis/saved_data/MODIS_zoomed_preds.RData")









##### prediction testing:  use prediction locations to measure mse  ####
print("Running predictions on observed locations!"); t_start = Sys.time()

# get locations and z for preds
pred_locs = locs[pred_idx,]
z_pred = z[pred_idx]

# calculate trend for prediction locations
trend_pred = pred_locs
trend_pred[,1]=1
XB_pred = trend_pred%*%beta

# make predictions conditional on posterior using VL
t_start = Sys.time() # Time VL full pred
z_VLpseudo = post_zy$t - XB
nuggets_VLpseudo = post_zy$D
vecchia.approx.pred = vecchia_specify(sub_locs, m=global_m, locs.pred=pred_locs)
preds=vecchia_prediction(z_VLpseudo, vecchia.approx.pred, covparms, nuggets_VLpseudo)
t_end = Sys.time();  time_dur = as.double(difftime(t_end, t_start, units = "mins"))
paste("Finished ZY full pred in", time_dur, "minutes")

# make predictions conditional on posterior using LowRank
t_start = Sys.time() # Time LR
z_VLRpseudo = post_lr$t-XB
nuggets_VLRpseudo = post_lr$D
vecchia.approx.pred.lr = vecchia_specify(sub_locs, m=global_m_lr, locs.pred=pred_locs, conditioning = "firstm")
predsLR=vecchia_prediction(z_VLRpseudo, vecchia.approx.pred.lr, covparms, nuggets_VLRpseudo)
t_end = Sys.time();  time_dur = as.double(difftime(t_end, t_start, units = "mins"))
paste("Finished LR full predictions in", time_dur, "minutes")


#### calculate comparison metrics ####
# Interval score from Gneiting and Raftery
interval_score = function(post,preds, XB_pred,  z_pred ){
  upper_quant = post$data_link(qnorm(p=.95,
                                     mean = preds$mu.pred+XB_pred,
                                     sd = sqrt(preds$var.pred)))
  lower_quant = post$data_link(qnorm(p=.05,
                                     mean = preds$mu.pred+XB_pred,
                                     sd = sqrt(preds$var.pred)))
  interval_score_val = upper_quant - lower_quant+
    2/.1*(lower_quant-z_pred)*( lower_quant > z_pred)+
    2/.1*(z_pred-upper_quant)*(upper_quant < z_pred)
  #return(mean((upper_quant > z_pred)*(lower_quant < z_pred)))
  return(interval_score_val)
}

## Assuming lognormal distribution for mean: compare truth to pred mean (lognormal median) and lognormal mean
mean_square_median_error_ZY = mean((exp(preds$mu.pred+XB_pred)-z_pred)^2)
mean_square_median_error_LR = mean((exp(predsLR$mu.pred+XB_pred)-z_pred)^2)

mean_square_error_ZY =  mean((exp(preds$mu.pred+XB_pred + 1/2*preds$var.pred)-z_pred)^2)
mean_square_error_LR =  mean((exp(predsLR$mu.pred+XB_pred + 1/2*preds$var.pred)-z_pred)^2)

interval_score_ZY =  mean(interval_score(post_zy, preds, XB_pred, z_pred))
interval_score_LR = mean(interval_score(post_lr, predsLR, XB_pred, z_pred))

print(c(mean_square_median_error_ZY,
        mean_square_median_error_LR,
        mean_square_error_ZY,
        mean_square_error_LR,
        interval_score_ZY,
        interval_score_LR))

# save predictions and resulting scores
save(pred_locs, z_pred, preds, predsLR, XB_pred,
     mean_square_mode_error_ZY, mean_square_mode_error_LR,
     mean_square_error_ZY, mean_square_error_LR,
     interval_score_ZY, interval_score_LR,
     file = "server/MODIS_analysis/saved_data/MODIS_full_preds_errors.RData")


# plots at prediction locations
pdf("server/MODIS_analysis/paper_plots/MODIS_Compare_preds_VL_LR.pdf", width = 8, height = 4)
n_pts = 128
par(mfrow=c(1,3), mar=rep(2,4), oma=rep(1.4,4))
quilt.plot(pred_locs, exp(predsLR$mu.pred+XB_pred), zlim=c(.01, 2.9), main = "Predictions LR (m=89)", nx = n_pts, ny=n_pts, add.legend = FALSE)
quilt.plot(pred_locs, z_pred, zlim=c(.01, 2.9), main = "Data, Unobserved", nx = n_pts, ny=n_pts, add.legend = FALSE)
rect(1000, 0, 1100, 100, border = "white")  # for zoomed predictions, calculated earlier
quilt.plot(pred_locs,  exp(preds$mu.pred+XB_pred), zlim=c(.01,2.9), main = "Predictions VL (m=20)", nx = n_pts, ny=n_pts, add.legend = TRUE)
dev.off()






# Calculated coarse gridded predictions, not in paper
print("Running predictions on grid!"); t_start = Sys.time()
grid.oneside=seq(1,1353,length=256)
grid.secside=seq(1,2029, length=256)
offset_locs=as.matrix(expand.grid(grid.oneside,grid.secside)) # grid of pred.locs

trend_pred = offset_locs
trend_pred[,1]=1
XB_pred = trend_pred%*%beta

t_start = Sys.time() # Time ZY grid pred
z_VLpseudo = post_zy$t - XB
nuggets_VLpseudo = post_zy$D
vecchia.approx.pred = vecchia_specify(sub_locs, m=global_m, locs.pred=offset_locs)
preds_grid=vecchia_prediction(z_VLpseudo, vecchia.approx.pred, covparms, nuggets_VLpseudo)
t_end = Sys.time();  time_dur = as.double(difftime(t_end, t_start, units = "mins"))
paste("Finished ZY grid pred pred in", time_dur, "minutes")



t_start = Sys.time() # Time LR
z_VLRpseudo = post_lr$t-XB
nuggets_VLRpseudo = post_lr$D
vecchia.approx.pred.lr = vecchia_specify(sub_locs, m=global_m_lr, locs.pred=offset_locs, conditioning = "firstm")
predsLR_grid=vecchia_prediction(z_VLRpseudo, vecchia.approx.pred.lr, covparms, nuggets_VLRpseudo)
t_end = Sys.time();  time_dur = as.double(difftime(t_end, t_start, units = "mins"))
paste("Finished LR gridded predictions in", time_dur, "minutes")


# plots for grid

pdf("server/MODIS_analysis/saved_data/MODIS_Compare_preds_grid_VL_LR.pdf", width = 12, height = 5)
par(mfrow=c(1,3))
quilt.plot(sub_locs, z_sub, zlim=c(.01, 2.9), main = "Obs", nx = 256, ny=256)
quilt.plot(offset_locs, exp(preds_grid$mu.pred+XB_pred), zlim=c(.01,2.9), main = "exp(VL)", nx = 256, ny=256)
quilt.plot(offset_locs, exp(predsLR_grid$mu.pred+XB_pred), zlim=c(.01,2.9),main = "exp(LR)", nx = 256, ny=256)
dev.off()

save(offset_locs,sub_locs, z_sub, preds_grid, predsLR_grid, XB_pred, file = "server/MODIS_analysis/saved_data/MODIS_grid_preds.RData")


