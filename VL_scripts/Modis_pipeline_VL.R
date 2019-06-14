source("VL_scripts/importer.R") #library(GPVecchia)
require(ggplot2)

load("MODIS/water_vapor_20190328.RData")
locs = as.matrix(locs)
n = length(locs[,1])

global_m = 20 #5 10 20 40
n_sample =2.5e5

# test on subset
set.seed(123)
sub_idx = sample(length(z), n_sample, replace = FALSE)
pred_idx = sample(setdiff(1:length(z), sub_idx), n_sample, replace = FALSE)

sub_locs = locs[sub_idx,]
z_sub = z[sub_idx]
#quilt.plot(sub_locs[,1], sub_locs[,2], z_sub, nx=128, ny =128)



#### trend estimation ####
X = sub_locs
X[,1]= 1



# glm: IRLS method from taking deriv of llh
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
  covparms_prev = c(.8, 40,   1.3)
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
  save(found_param, file = "VL_scripts/saved_data/covparms_mXX_nXX.RData")
}


#### Plot parameter estimation results ####
param_df = read.csv("VL_scripts/saved_data/measurements.csv")
param_df$n<-as.factor(param_df$n)
# remove the time column
param_df = param_df[,-7]
# melt data for facet plot
melted_params = melt(param_df, id.vars = c("n", "m"), measure.vars =c(3,4,5,6))
levels(melted_params$variable) <- c("a", "sigma^2", "rho", "nu")
ggplot(melted_params, aes(x = m, y = value, color = n, shape =n, linetype = n))+
  geom_line() + geom_point() +theme_bw() + theme(legend.position = "top")+ scale_shape(solid=FALSE)+
  facet_wrap(~variable,scales = "free_y", ncol=4,labeller = label_parsed)
ggsave("VL_scripts/plots/MODIS_param.pdf", width = 10, height = 3)




#### Posterior estimation and predictions  ####
# set parameters found in exploration
if(FALSE){
  a = 0.89
  covparms = c(0.25, 31,   3)
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
  #vecchia.approx.lr = vecchia_specify(sub_locs, m=global_m_lr, conditioning = "firstm")
  vecchia.approx.lr = vecchia_specify(sub_locs, m=global_m_lr, conditioning = "mra",
                                      mra.options = list(r=c(global_m_lr,1)))
  post_lr = calculate_posterior_VL(z_sub, vecchia.approx.lr, "gamma" , covparms, likparms = default_lh_params, prior_mean = XB)
  t_end = Sys.time();  time_dur = as.double(difftime(t_end, t_start, units = "mins"))
  paste("Finished LR posterior in", time_dur, "minutes")


  ## Create and save comparison plot
  pdf("VL_scripts/plots/MODIS_Compare_post_VL_LR.pdf", width = 12, height = 5)
  npixels = 128
  par(mfrow=c(1,3))
  quilt.plot(sub_locs, z_sub, zlim=c(.01, 2.9), main = "Obs", nx = npixels, ny=npixels)
  quilt.plot(sub_locs, exp(post_zy$mean), zlim=c(.01,2.9), main = "Posterior exp(VL)", nx = npixels, ny=npixels)
  quilt.plot(sub_locs, exp(post_lr$mean), zlim=c(.01,2.9),main = "Posterior exp(LR)", nx = npixels, ny=npixels)
  dev.off()

  # save objects for future plotting. May be >100MB
  save(sub_locs, z_sub, post_zy, post_lr, time_dur,
       file = "VL_scripts/saved_data/MODIS_posteriors.RData")
  paste("Finished posteriors in", time_dur, "minutes")



  #####  Zoom ####
  ## Make predictions on a zoomed region, to see small scale structure
  print("Running predictions on zoomed observations!"); t_start = Sys.time()
  # load("VL_scripts/saved_data/MODIS_posteriors.RData")
  zx1 = 800
  zx2 = 1000
  zy1 = 200
  zy2 = 400


  # Subset locations
  zoom_idx = intersect(which(locs[,1] >zx1 & locs[,1]<zx2),  which(locs[,2] >zy1 & locs[,2] < zy2))
  # remove obs locs, for comparison plot
  zoom_idx = setdiff(zoom_idx, sub_idx)
  locs_zoom = locs[zoom_idx,]
  locs_zzoom=locs_zoom
  z_zoom =z[zoom_idx]

  # some locations are missing from the data, fill in gaps for prediction
  full_set_zoom = expand.grid(seq(zx1+.5, zx2-.5), seq(zy1+.5, zy2-.5))
  colnames(full_set_zoom)=c("x","y")
  dup_array = duplicated(rbind(sub_locs,full_set_zoom))
  # take all rows that aren't duplicates of observations
  locs_zoom = full_set_zoom[which(dup_array[(length(sub_locs[,1])+1):length(dup_array)]==FALSE),]
  locs_zoom = as.matrix(locs_zoom)
  #quilt.plot(locs_zoom, z_zoom, nx=100, ny=100)

  # Calculate trend over prediction region
  trend_zoom = locs_zoom
  trend_zoom[,1]=1
  XB_zoom = trend_zoom%*%beta

  ## Predict on zoomed region using VL
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
  #vecchia.approx.pred.zoom.lr = vecchia_specify(sub_locs, m=global_m_lr, locs.pred=locs_zoom, conditioning = "firstm")
  vecchia.approx.pred.zoom.lr = vecchia_specify(sub_locs, m=global_m_lr, locs.pred=locs_zoom, conditioning = "mra",
                                                mra.options = list(r=c(global_m_lr,1)))


  predsLR_zoom=vecchia_prediction(z_VLRpseudo, vecchia.approx.pred.zoom.lr, covparms, nuggets_VLRpseudo)
  t_end = Sys.time();  time_dur = as.double(difftime(t_end, t_start, units = "mins"))
  paste("Finished LR zoom pred in", time_dur, "minutes")



  ##  obs within zoomed locations
  zoom_idx_vis = intersect(which(sub_locs[,1] >zx1 & sub_locs[,1]<zx2),  which(sub_locs[,2] >zy1 & sub_locs[,2] < zy2))
  sub_locs_zoom = as.matrix(sub_locs[zoom_idx_vis,])
  z_sub_zoom = as.matrix(z_sub[zoom_idx_vis])



  ## Zoomed plots
  pdf("VL_scripts/plots/MODIS_compare_zoom_VL_LR.pdf", width = 9, height = 8)
  npts = 200
  z_lims =c(.01, .9)
  par(mfrow=c(2,2), mar = rep(2,4))
  quilt.plot(sub_locs_zoom, z_sub_zoom, zlim=z_lims, main = "Observed Training Data", nx = npts, ny=npts,
             add.legend = FALSE, xaxt='n', yaxt = "n")
  quilt.plot(locs_zzoom, z_zoom, zlim=z_lims, main = "Test Data", nx = npts, ny=npts, nlevel=1024)
  par(mar = c(2,1.5,2,1.5))
  quilt.plot(locs_zoom, exp(predsLR_zoom$mu.pred+XB_zoom+1/2*predsLR_zoom$var.pred),
             zlim=z_lims,main = "Predictions LR (m=89)", nx = npts, ny=npts,
             add.legend = FALSE,xaxt='n', yaxt = "n", nlevel=1024)
  quilt.plot(sub_locs_zoom, exp(predsLR_zoom$mu.obs+XB+1/2*predsLR_zoom$var.obs)[zoom_idx_vis],
             zlim=z_lims, nx = npts, ny=npts, add.legend = FALSE, add = TRUE)


  quilt.plot(locs_zoom, exp(preds_zoom$mu.pred+XB_zoom+1/2*preds_zoom$var.pred), zlim=z_lims, main = "Predictions VL (m=20)",
             nx = npts, ny=npts, add.legend = FALSE, xaxt='n', yaxt = "n", nlevel=1024)

  quilt.plot(sub_locs_zoom, exp(preds_zoom$mu.obs+XB+1/2*preds_zoom$var.obs)[zoom_idx_vis], zlim=z_lims,xlim = c(zx1, zx2),
             ylim = c(zy1, zy2), nx = npts, ny=npts,
             add.legend = FALSE, add = TRUE, nlevel=1024)

  dev.off()

  # save prediction objects for future plotting
  save(locs_zoom,locs_zzoom, z_zoom, sub_locs_zoom, z_sub_zoom,zoom_idx_vis,
       preds_zoom, predsLR_zoom,XB, XB_zoom, file = "VL_scripts/saved_data/MODIS_zoomed_preds.RData")





  ##### Prediction Test####
  ##  use prediction locations to measure mse, CRPS
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


  #### Prediction Scoring  ####
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


  get_crps = function(z, preds, a, XB_pred, samples_to_gen = 10){
    # get mean from prediction method
    z = z
    pred_mu = (preds$mu.pred+XB_pred)
    pred_sd = sqrt(pmax(preds$var.pred, 0))

    sample_mat = matrix(0, ncol=samples_to_gen, nrow = length(z))
    for(smpl_idx in 1:samples_to_gen){
      sample_latent_mean = rnorm(n = length(z), mean =  pred_mu, sd =  pred_sd)
      sample_pred_mean = exp(sample_latent_mean)
      # use CRPS from package to score the prediction
      sample_mat[,smpl_idx]= rgamma(length(z), shape =a, rate =   a/sample_pred_mean )
    }

    return( scoringRules::crps_sample(y=z, dat=sample_mat) )
  }

  crps_vl = mean(get_crps(z_pred, preds, a, XB_pred, 200))
  crps_lr = mean(get_crps(z_pred, predsLR, a, XB_pred, 200))

  ## Marginal mean: compare truth to marginal predictive median and mean
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
          crps_vl,
          crps_lr))

  # save predictions and resulting scores
  save(pred_locs, z_pred, preds, predsLR, XB_pred,
       mean_square_median_error_ZY, mean_square_median_error_LR,
       mean_square_error_ZY, mean_square_error_LR,
       crps_vl,crps_lr,
       file = "VL_scripts/saved_data/MODIS_full_preds_errors.RData")


  # plots at prediction locations
  pdf("VL_scripts/plots/MODIS_Compare_preds_VL_LR.pdf", width = 9, height = 4)
  n_pts = 128
  par(mfrow=c(1,3), mar=rep(2,4), oma=rep(1.4,4))
  quilt.plot(pred_locs, exp(predsLR$mu.pred+XB_pred+1/2*predsLR$var.pred), zlim=c(.01, 2.9),  main = "Predictions LR (m=89)", nx = n_pts, ny=n_pts, add.legend = FALSE, nlevel=1024)
  quilt.plot(pred_locs, z_pred, zlim=c(.01, 2.9), main = "Test Data", nx = n_pts, ny=n_pts, add.legend = FALSE)
  rect(1000, 0, 1200, 200, border = "white")  # for zoomed predictions, calculated earlier
  quilt.plot(pred_locs,  exp(preds$mu.pred+XB_pred+1/2*preds$var.pred), zlim=c(.01,2.9), main = "Predictions VL (m=20)", nx = n_pts, ny=n_pts, add.legend = TRUE)
  dev.off()


  ##### Likelihood comparison ####
  # compare approximate likelihoods to justify m=20

  llh2 = matrix(0, nrow = 6, ncol=3)
  covparms = c(0.25, 31,  3)
  mvals = c(1, 2, 5,10,20,40)
  default_lh_params = list("alpha"=.89, "sigma"=sqrt(.1))

  for(m_idx in 1:6){
    m = mvals[m_idx]
    vecchia.approx = vecchia_specify(sub_locs, m=m, cond.yz = "zy")
    vecchia.approx.IW = vecchia_specify(sub_locs, m=m)

    vll = vecchia_laplace_likelihood(z_sub,
                                     vecchia.approx,
                                     likelihood_model="gamma",
                                     covparms = covparms,
                                     return_all = FALSE,
                                     likparms = default_lh_params,
                                     prior_mean = XB,
                                     vecchia.approx.IW=vecchia.approx.IW,
                                     y_init = NA )

    vecchia.approx.lr = vecchia_specify(sub_locs, m=m, conditioning = "mra", mra.options = list(r=c(m,1)))

    lrll = vecchia_laplace_likelihood(z_sub,
                                      vecchia.approx.lr,
                                      likelihood_model="gamma",
                                      covparms = covparms,
                                      return_all = FALSE,
                                      likparms = default_lh_params,
                                      prior_mean = XB,
                                      vecchia.approx.IW=NA,
                                      y_init = NA )
    print(c(m, vll, lrll))
    llh2[m_idx,] = c(m, vll, lrll)
  }

  vecchia.approx.lr = vecchia_specify(sub_locs, m=89, conditioning = "mra", mra.options = list(r=c(89,1)))
# special case for low rank:  m=89
  lrll_89 = vecchia_laplace_likelihood(z_sub,
                                    vecchia.approx.lr,
                                    likelihood_model="gamma",
                                    covparms = covparms,
                                    return_all = FALSE,
                                    likparms = default_lh_params,
                                    prior_mean = XB,
                                    vecchia.approx.IW=NA,
                                    y_init = NA )

  # columns of llh2
  m = c(1,2,5,10,20, 40)
  lr = c(-93916.26, -93721.05, -93509.21, -93014.81,  -91873.01, -87465.02)
  vl = c(-53210.20,  -53517.53,  -53486.39,  2.5e5/1e5*-22426.45 ,  -54214.64,  -54215.30)

  #n=1e5: 2.5e5/1e5*-22426.45 vs -53532.99 (different seed)
  pdf("VL_scripts/plots/MODIS_vll.pdf", width=5, height = 5)
  par(mar=c(4,4,2,1))
  plot(m, -lr, type ="b", ylim = c(min(-vl), max(-lr)),
       main = "Approximate Log Likelihood",
       xlab = "m", ylab = "Neg. LLH", col=2)
  points(m, -vl, type = "b", col=1)
  points(20, 77988.55, type = "p", col=2, pch = 10)
  points(20, 54214.64, type = "p", col=1, pch = 10)
  legend("right", legend = c("VL", "LowRank","LowRank: m=89"), col = c(1,2,2), pch = c(1,1,10), bty = "o")
  dev.off()
}
