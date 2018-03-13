#setwd("/home/grad/dzilber/Documents/dz_VL_simulation/")
#source("VL_method_cpp.R")


##########################################################################
##################### Compare models: computation time  ########################
##########################################################################


dposterior = function(y, pred){
  mu = pred$mean
  W = pred$W
  quad_term = -t(y-mu) %*% W %*% (y-mu)/2
  det_term = sum(log(diag(t(chol(W)))))# log(det(W))/2
  pi_term = -length(y)/2*log(2*pi)
  #print(paste("[Quad", quad_term, "] [Det", det_term,"]"))
  # value is of class "dgeMatrix"
  return((quad_term+det_term+pi_term)[1,1])
}


create_scenario_tester = function(log_file_name=NA){
  header = c("Mod", "Domain", "Dimen", "Sample", "C_Smoothness", "C_Range","Neighbors","Seed_off", 
             "MSE_Laplace", "MSE_VL",  "MSE_LowRank", "LS_Laplace", "LS_VL",  "LS_LowRank", 
             "Time_Laplace", "Time_VL",  "Time_LowRank", 
             "Iter_Laplace", "Iter_VL",  "Iter_LowRank")
  
  num_cols= length(header)
  
  filename = log_file_name
  if (is.na(log_file_name)) filename= paste("results_",as.Date(Sys.time()),".csv", sep="")
  
  # Write to file rather than memory, avoid interruption issues
  write(header, file =filename, ncolumns=num_cols, sep=",")
  
  run_scenario = function(seed_r, domn, dimen, samp_size, neighbors, smth, mod_type, rnge, show_output = FALSE){
    
    covparms=c(1, rnge, smth)
    covfun <- function(locs) Matern(rdist(locs), range = rnge, smoothness = smth)
    
    if(mod_type==1)
      ls=gauss_sample(samp_size,covfun, seed = samp_size+seed_r+mod_type , dom = domn, dimen=dimen)
    if(mod_type==2)
      ls=logistic_sample(samp_size,covfun, seed = samp_size+mod_type+seed_r, dom = domn, dimen=dimen)
    if(mod_type==3)
      ls=pois_sample(samp_size,covfun, seed = samp_size+seed_r+mod_type, dom = domn, dimen = dimen)
    if(mod_type==4)
      ls=gamma_sample(samp_size,covfun, seed = samp_size+seed_r+mod_type , dom = domn, dimen=dimen)
    
    
    
    t_start = Sys.time()
    pred_y_lv = estimate_laplace_vec_cpp(ls,covparms, m=neighbors)
    t_start_2 = Sys.time()
    pred_y_lr =  estimate_laplace_vec_cpp(ls, covparms, m=neighbors, use_low_rank = TRUE)
    t_start_3 = Sys.time()
    pred_y_l= estimate_laplace(ls, C =covfun(ls$locs))
    t_fin = Sys.time()
    
    LV_time = as.double(difftime(t_start_2, t_start, units = "mins"))
    LR_time = as.double(difftime(t_start_3, t_start_2, units = "mins"))
    Lap_time = as.double(difftime(t_fin, t_start_3, units = "mins"))
    sim_times = c(Lap_time, LV_time, LR_time)
    #print(sim_times)
    lv_score = dposterior(ls$y, pred_y_lv)
    lr_score = dposterior(ls$y, pred_y_lr)
    laplace_score = dposterior(ls$y, pred_y_l)
    log_score_results <-c(laplace_score, lv_score, lr_score)
    
    mse  =  c(mean((ls$y-pred_y_l$mean)^2), mean((ls$y-pred_y_lv$mean)^2), mean((ls$y-pred_y_lr$mean)^2))
    NR_iters = c(pred_y_l$iter, pred_y_lv$iter, pred_y_lr$iter)
    
    
    results_i = c(mod_type, domn, dimen, samp_size, smth, rnge, neighbors, seed_r, mse, log_score_results, sim_times, NR_iters)
    
    if (show_output){
      mod_name = ifelse(mod_type==1, "logistic", ifelse(mod_type==2, "Poisson", ifelse(mod_type==3, "Gaussian", "Gamma")))
      results_list = list(mod_name, domn, dimen, samp_size, smth, rnge, neighbors, seed_r, mse, log_score_results, sim_times,NR_iters)
      output_data = c("model", "domain", "dimen", "sample", "smth", "range", "neighbors", "seed_r", "mse", "log_score_results", "sim_times", "NR_iter")
      named_output = setNames( results_list, c(output_data))
      message(str(named_output))
    }  
    
    write(results_i, file = filename, ncolumns=num_cols, append=TRUE, sep=",")
  }
  return(run_scenario)
}



# Sample Code:
# 1.  Create scenario runner by specifying a name for the batch of scenarios
# scenario_runner = create_scenario_tester("testing12")
#2.  For every scenario, run the simulationa nd record results to the log file with the given name
# scenario_runner(seed_r, domn, samp_size, neighbors, smth, mod_type, covfun, show_output = FALSE)
# writes the output to a log, so no need to save any object
