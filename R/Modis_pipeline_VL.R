source("server/importer.R")

load("server/MODIS_analysis/water_vapor_20190328.RData")
locs = as.matrix(locs)
n = length(locs[,1])

global_m = 40 # 20 40
n_sample = 1e4

# test on subset or not
set.seed(123)
sub_idx = sample(length(z), n_sample, replace = FALSE)
pred_idx = setdiff(1:length(z), sub_idx)

sub_locs = locs[sub_idx,]
z_sub = z[sub_idx]
quilt.plot(sub_locs[,1], sub_locs[,2], z_sub, nx=256, ny =256)



#### trend estimation ####
X = sub_locs
X[,1]= 1



# glm: IRLS method from taking deriv of llh
#IRLS
beta = c(-1.5,0.001)
for(i in 1:10){
  XB = X%*% beta
  W = -diag.spam(array(exp(-XB)*z_sub))
  A  = exp(-XB)*z_sub-1
  U  = W%*% XB - A
  beta = solve( t(X) %*% W %*% X , t(X) %*% U)
}
print(beta)
XB = X%*% beta


### Covparm and shape estimation ####
# shape estimation: maximize conditional likelihood.  Integrated likelihood diverges
update_a = function(a_init, covparms,  vecchia.approx, vecchia.approx.IW, XB ){
  a_prev =a_init
  for(i in 1:7){
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
    sprintf("Evaluating covparms = (%.4f %.4f %.4f)",
            covparms[1], covparms[2],covparms[3])
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
  res = optim(x0, vl_likelihood, method = "Nelder-Mead", control = list("maxit" = 500, "reltol" = 1e-5))
  print(res$convergence); print(exp(res$par))
  return(exp(res$par))
}


send_result = function(a ,covparms,m, time_dur,iter_count, subject=NA){

  # emailing stuff:  for fun
  require(gmailr)
  require(mime)
  email_subject = sprintf("MODIS param estimation for m= %i, n=%i in %.1f minutes",m,n_sample,time_dur)
  if(!is.na(subject)) email_subject=subject
  email_text = sprintf("Parameters \nm= %i  \na = %.3f \ncovparms = (%.4f %.4f %.4f) \niter: %i",
                       m, a, covparms[1], covparms[2],covparms[3], iter_count)
  #file = "server/MODIS_analysis/saved_data/covparms_m.RData"
  mime() %>%
    to("dzilber@tamu.edu") %>%
    from("dzilber@tamu.edu") %>%
    subject(email_subject) %>%
    text_body(email_text) -> text_msg
  #attach_file(filename = file)
  send_message(text_msg)

  return()
}




## Do parameter estimation for multiple m values

  #param_table = matrix(0, nrow =length(m_vals), ncol =5)
  m_rep = global_m
  print(paste("Estimating parameters for m=", m_rep))

  print("Step 1, generating vecchia approximations")
  vecchia.approx = vecchia_specify(sub_locs, m=m_rep, cond.yz = "zy")
  vecchia.approx.IW = vecchia_specify(sub_locs, m=m_rep)

  ## Iterative method:  estimate a, then covparms, then a again
  print("Step 2, optimizing parameters")
  a_prev = 0.9
  covparms_prev = c(1, 1000,   1.3)
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
      iter = i
      break
    }
    a_prev = a
    covparms_prev = covparms
  }
  t_end = Sys.time()
  time_dur = as.double(difftime(t_end, t_start, units = "mins"))
  send_result(a, covparms, m_rep, time_dur, iter_count)
  cat(a, covparms, m_rep, time_dur, iter_count)
  #param_table[m_idx,] = c(a,covparms, iter_count)



  #param_df = data.frame(cbind(m_vals,param_table))
  #colnames(param_df) <- c("m", "a", "sig","rho","nu", "iters" )
  #save(param_df, file = "server/MODIS_analysis/saved_data/covparms_m.RData")




 ### Do predictions  ####
if(FALSE){
  vecchia.approx = vecchia_specify(sub_locs, m=global_m, cond.yz = "zy")
  vecchia.approx.lr = vecchia_specify(sub_locs, m=global_m, conditioning = "firstm")
  a = 0.9
  covparms = c(1.04, 1000,   1.432)
  default_lh_params = list("alpha"=a, "sigma"=sqrt(.1))
  post_zy = calculate_posterior_VL(z_sub, vecchia.approx, "gamma" , covparms, likparms = default_lh_params, prior_mean = XB)
  post_lr = calculate_posterior_VL(z_sub, vecchia.approx.lr, "gamma" , covparms, likparms = default_lh_params, prior_mean = XB)

  pdf("server/MODIS_analysis/saved_data/MODIS_Compare_post_VL_LR.pdf", width = 12, height = 5)
  par(mfrow=c(1,3))
  quilt.plot(sub_locs, z_sub, zlim=c(.01, 2.9), main = "Obs", nx = 256, ny=256)
  quilt.plot(sub_locs, exp(post_zy$mean), zlim=c(.01,2.9), main = "Posterior exp(VL)", nx = 256, ny=256)
  quilt.plot(sub_locs, exp(post_lr$mean), zlim=c(.01,2.9),main = "Posterior exp(LR)", nx = 256, ny=256)
  dev.off()
  save(post_zy, post_lr, file = "server/MODIS_analysis/saved_data/MODIS_posteriors.RData")

  # predictions v1:  grid

  grid.oneside=seq(1,1353,length=256)
  grid.secside=seq(1,2029, length=256)
  offset_locs=as.matrix(expand.grid(grid.oneside,grid.secside)) # grid of pred.locs

  trend_pred = offset_locs
  trend_pred[,1]=1
  XB_pred = trend_pred%*%beta

  z_VLpseudo = post_zy$t - XB
  nuggets_VLpseudo = post_zy$D
  vecchia.approx.pred = vecchia_specify(sub_locs, m=global_m, locs.pred=offset_locs)
  preds_grid=vecchia_prediction(z_VLpseudo, vecchia.approx.pred, covparms, nuggets_VLpseudo)

  z_VLRpseudo = post_lr$t-XB
  nuggets_VLRpseudo = post_lr$D
  vecchia.approx.pred.lr = vecchia_specify(sub_locs, m=global_m, locs.pred=offset_locs, conditioning = "firstm")
  predsLR_grid=vecchia_prediction(z_VLRpseudo, vecchia.approx.pred.lr, covparms, nuggets_VLRpseudo)

  # plots v1

  pdf("server/MODIS_analysis/saved_data/MODIS_Compare_preds_grid_VL_LR.pdf", width = 12, height = 5)
  par(mfrow=c(1,3))
  quilt.plot(sub_locs, z_sub, zlim=c(.01, 2.9), main = "Obs", nx = 256, ny=256)
  quilt.plot(offset_locs, exp(preds_grid$mu.pred+XB_pred), zlim=c(.01,2.9), main = "exp(VL)", nx = 256, ny=256)
  quilt.plot(offset_locs, exp(predsLR_grid$mu.pred+XB_pred), zlim=c(.01,2.9),main = "exp(LR)", nx = 256, ny=256)
  dev.off()

  save(preds_grid, predsLR_grid, XB_pred, file = "server/MODIS_analysis/saved_data/MODIS_grid_preds.RData")







  # predictions v2:  use prediction locations to measure mse
  pred_locs = as.matrix(locs[pred_idx,])
  z_pred = as.matrix(z[pred_idx])
  trend_pred = pred_locs
  trend_pred[,1]=1
  XB_pred = trend_pred%*%beta

  z_VLpseudo = post_zy$t - XB
  nuggets_VLpseudo = post_zy$D
  vecchia.approx.pred = vecchia_specify(sub_locs, m=global_m, locs.pred=pred_locs)
  preds=vecchia_prediction(z_VLpseudo, vecchia.approx.pred, covparms, nuggets_VLpseudo)

  z_VLRpseudo = post_lr$t-XB
  nuggets_VLRpseudo = post_lr$D
  vecchia.approx.pred.lr = vecchia_specify(sub_locs, m=global_m, locs.pred=pred_locs, conditioning = "firstm")
  predsLR=vecchia_prediction(z_VLRpseudo, vecchia.approx.pred.lr, covparms, nuggets_VLRpseudo)


  ## Calculate metrics:
  interval_score = function(post, preds, z_pred ){
    upper_quant = post$data_link(qnorm(p=.95,
                                       mean = preds$mu.pred,
                                       sd = sqrt(preds$var.pred)))
    lower_quant = post$data_link(qnorm(p=.05,
                                       mean = preds$mu.pred,
                                       sd = sqrt(preds$var.pred)))
    interval_score_val = upper_quant - lower_quant+
      2/.1*(lower_quant-z_pred)*( lower_quant > z_pred)+
      2/.1*(z_pred-upper_quant)*(upper_quant < z_pred)
    #return(mean((upper_quant > z_pred)*(lower_quant < z_pred)))
    return(interval_score_val)
  }




  mean_square_median_error_ZY = mean((exp(preds$mu.pred+XB_pred)-z_pred)^2)
  mean_square_median_error_LR = mean((exp(predsLR$mu.pred+XB_pred)-z_pred)^2)

  mean_square_error_ZY =  mean((exp(preds$mu.pred+XB_pred + 1/2*preds$var.pred)-z_pred)^2)
  mean_square_error_LR =  mean((exp(predsLR$mu.pred+XB_pred + 1/2*preds$var.pred)-z_pred)^2)

  interval_score_ZY =  mean(interval_score(post_zy, preds, z_pred))
  interval_score_LR = mean(interval_score(post_lr, predsLR, z_pred))

  print(c(mean_square_mode_error_ZY,
          mean_square_mode_error_LR,
          mean_square_error_ZY,
          mean_square_error_LR,
          interval_score_ZY,
          interval_score_LR))

  save(preds, predsLR, XB_pred,
       mean_square_mode_error_ZY, mean_square_mode_error_LR,
       mean_square_error_ZY, mean_square_error_LR,
       interval_score_ZY, interval_score_LR,
       file = "server/MODIS_analysis/saved_data/MODIS_preds_errors.RData")


  send_result(a,covparms, 10, time_dur = -1, iter_count = -1, subject = "Finished all predictions")


  # plots v2
  pdf("server/MODIS_analysis/saved_data/MODIS_Compare_preds_VL_LR.pdf", width = 12, height = 5)
  par(mfrow=c(1,3))
  quilt.plot(pred_locs, z_pred)#, zlim=c(.01, 2.9), main = "Obs", nx = 256, ny=256)
  quilt.plot(pred_locs, pred_locs,  exp(preds$mu.pred+XB_pred), zlim=c(.01,2.9), main = "exp(VL)", nx = 256, ny=256)
  quilt.plot(pred_locs, exp(predsLR$mu.pred+XB_pred), zlim=c(.01,2.9),main = "exp(LR)", nx = 256, ny=256)
  dev.off()

}


