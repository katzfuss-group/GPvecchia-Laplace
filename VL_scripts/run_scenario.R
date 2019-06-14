# This code takes a set of parameters defining some simulation scenario,
# such as data model, covariance parameters, and sample size, then simulates
# data and compares the resulting estimates of VL to other methods.


##########################################################################
##################### Compare models: computation time  ########################
##########################################################################


dposterior = function(y, pred){
  mu = pred$mean
  W= pred$W
  # VL calculates V, Laplace calculates W
  if("V" %in% names(pred)){
    det_term = sum(log(diag(pred$V)))
  }else{
    V=t(chol(forceSymmetric(pred$W)))
    #V = t(chol((W + t(W))/2))
    det_term = sum(log(diag(V)))
    #det_term =log(det(W))/2
  }
  quad_term = -t(y-mu) %*% W %*% (y-mu)/2
  pi_term = -length(y)/2*log(2*pi)
  # value is of class "dgeMatrix"
  return((quad_term+det_term+pi_term)[1,1])
}


create_scenario_tester = function(header, log_file_name){

  num_cols= length(header)

  # Write to file rather than memory, avoid losing data with crash
  filename = log_file_name
  write_to_file = !is.na(log_file_name)
  if (write_to_file) write(header, file =filename, ncolumns=num_cols, sep=",")
  #filename= paste("results_",as.Date(Sys.time()),".csv", sep="")


  run_scenario = function(seed_r, domn, dimen, samp_size, neighbors, smth, mod_type, rnge, show_output = FALSE, run_algo= TRUE){
    max_tries =  ifelse (seed_r<0, 3, 1) # in case of error, try different seed or exit without writing
    failed_convergence = FALSE
    for(i in 1:max_tries) {

      covparms=c(1, rnge, smth)
      covfun <- function(locs) Matern(rdist(locs), range = rnge, smoothness = smth)
      covmodel <- RMmatern(smth,var = 1,scale = rnge)

      if(seed_r<0 || i>1) seed_r = sample(1e6,1)

      if(mod_type==1)
        ls=gauss_sample(samp_size,covmodel, seed = seed_r , dom = domn, dimen=dimen)
      if(mod_type==2)
        ls=logistic_sample(samp_size,covmodel, seed = seed_r, dom = domn, dimen=dimen)
      if(mod_type==3)
        ls=pois_sample(samp_size,covmodel, seed = seed_r, dom = domn, dimen = dimen)
      if(mod_type==4)
        ls=gamma_sample(samp_size,covmodel, seed = seed_r , dom = domn, dimen=dimen)

      vecchia.approx=vecchia_specify(ls$locs, neighbors)
      vecchia.approx_z=vecchia_specify(ls$locs, neighbors, cond.yz = "zy")
      vecchia.approx_LR=vecchia_specify(ls$locs, neighbors, conditioning="firstm")

      default_out = list("W"=1, "mean" = 0, "runtime" = -1, "iter"=-1, "cnvgd" = TRUE)
      pred_y_lv =  default_out # SGVecchia Laplace
      pred_y_lv_z =  default_out # Observed Vecchia Laplace
      pred_y_l = default_out # Laplace
      pred_y_lr = default_out # low rank Laplace

      pred_y_lv = calculate_posterior_VL(ls$z, vecchia.approx, ls$type, covparms, return_all = TRUE)
      pred_y_lr =  calculate_posterior_VL(ls$z, vecchia.approx_LR, ls$type, covparms,  return_all = TRUE)
      pred_y_lv_z =  calculate_posterior_VL(ls$z, vecchia.approx_z, ls$type, covparms,  return_all = TRUE)
      if (run_algo) pred_y_l= calculate_posterior_laplace(ls$z, ls$type, C =covfun(ls$locs), return_all = TRUE)

      failed_convergence  = any(!c(pred_y_lv$cnvgd, pred_y_l$cnvgd, pred_y_lv_z$cnvgd, pred_y_lr$cnvgd))
      if(failed_convergence) {
        message("Failed convergence, check seed")
        if (show_output){
          mod_name = ifelse(mod_type==2, "Logistic", ifelse(mod_type==3, "Poisson", ifelse(mod_type==1, "Gaussian", "Gamma")))
          results_list = list(mod_name, domn, dimen, samp_size, smth, rnge, neighbors, seed_r)
          output_data = c("model", "domain", "dimen", "sample", "smth", "range", "neighbors", "seed_r")
          named_output = setNames( results_list, c(output_data))
          message(str(named_output))
        }
        next # try again or exit without writing result
      }
      sim_times = c(pred_y_l$runtime, pred_y_lv$runtime, pred_y_lv_z$runtime, pred_y_lr$runtime)

      #print(sim_times)
      lv_score = 0#dposterior(ls$y, pred_y_lv)
      lv_z_score = 0#dposterior(ls$y, pred_y_lv_z)
      lr_score = 0#dposterior(ls$y, pred_y_lr)
      laplace_score = 0#dposterior(ls$y, pred_y_l)
      log_score_results <-c(laplace_score, lv_score, lv_z_score, lr_score)

      mse  =  c(mean((ls$y-pred_y_l$mean)^2), mean((ls$y-pred_y_lv$mean)^2),
                mean((ls$y-pred_y_lv_z$mean)^2), mean((ls$y-pred_y_lr$mean)^2))
      NR_iters = c(pred_y_l$iter, pred_y_lv$iter, pred_y_lv_z$iter, pred_y_lr$iter)




      results_i = c(mod_type, domn, dimen, samp_size, smth, rnge, neighbors, seed_r, mse, log_score_results, sim_times, NR_iters)

      if (show_output){
        mod_name = ifelse(mod_type==2, "Logistic", ifelse(mod_type==3, "Poisson", ifelse(mod_type==1, "Gaussian", "Gamma")))
        results_list = list(mod_name, domn, dimen, samp_size, smth, rnge, neighbors, seed_r, mse, log_score_results, sim_times,NR_iters)
        output_data = c("model", "domain", "dimen", "sample", "smth", "range", "neighbors", "seed_r", "mse", "log_score_results", "sim_times", "NR_iter")
        named_output = setNames( results_list, c(output_data))
        message(str(named_output))
      }


      if(write_to_file) write(results_i, file = filename, ncolumns=num_cols, append=TRUE, sep=",")
      # successful iteration,return will break loop
      return(results_i) # for secondary output; avoids issues with buffered writing?
    }
  }
  return(run_scenario)
}


# Sample Code:
# 1.  Create scenario runner by specifying a name for the batch of scenarios
# scenario_runner = create_scenario_tester("testing12")
#2.  For every scenario, run the simulationa nd record results to the log file with the given name
# scenario_runner(seed_r, domn, samp_size, neighbors, smth, mod_type, covfun, show_output = FALSE)
# writes the output to a log, so no need to save any object
