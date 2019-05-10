
# import VL code somehow
source("server/importer.R")
t_start = Sys.time()

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



get_runner = function(loc, loc_chol){
  ## Begin MSE loop
  run_LGCP = function(m, range, smooth, run_laplace =TRUE, seed_val = NA, ds_factor = 2, locs = loc, loc_cov_chol = loc_chol){


    data_distr = "poisson"
    spatial.dim =2 # number of spatial dimensions
    dom = 1 # [0, domn]
    n=100^2 # high res locs



    # covariance parameters
    sig2=1;
    covparms = c(sig2, range, smooth)
    covfun <- function(locs) sig2*Matern(fields::rdist(locs),range=range,smoothness=smooth)


    # simulate location
    if(is.na(seed_val)){
      seed_val = sample(1e6,1)
    }
    set.seed(seed_val)



    # simulate latent process
    y= as.numeric(loc_cov_chol%*%rnorm(n))

    ###  Simulate poisson directly, convert to random points
    pt_density = n^(1/spatial.dim)
    c = dom / pt_density
    z = rpois(n, c*exp(y))
    locs_1 = locs[which(z==1),]
    locs_2 = locs[which(z==2),]
    locs_3 = locs[which(z==3),]
    locs_11 = locs_1+rnorm(length(locs_1))*c
    locs_22 = rbind(locs_2, locs_2)+rnorm(length(locs_2)*2)*c
    locs_33 = rbind(locs_3, locs_3, locs_3)+rnorm(length(locs_3)*3)*c
    pp_locs = rbind(locs_11, locs_22, locs_33)

    ### Downsample the obs and pool into larger pixels
    coarse_num = pt_density/ds_factor
    down_sampled = as.image(z, locs, nx = coarse_num, ny = coarse_num, FUN = sum)
    coarse_z = as.vector(t(down_sampled$z))
    coarse_idx = unique(down_sampled$ind)
    coarse_locs = cbind(down_sampled$x[coarse_idx[,1]], down_sampled$y[coarse_idx[,2]])


    #####################   specify Vecchia approx    #######################
    # (this only has to be run once)
    vecchia.approx = vecchia_specify(coarse_z, coarse_locs, m, cond.yz = "zy")
    #vecchia.approx.sgv = vecchia_specify(coarse_z, coarse_locs, m)
    vecchia.approx.ll=vecchia_specify(coarse_z, coarse_locs, m, conditioning = "firstm")

    #####################   prediction at observed locations    ######################
    # Perform inference on latent mean with Vecchia Laplace approximation


    ## Check:  W  = inv(Cov)
    posterior_VL = calculate_posterior_VL(vecchia.approx, likelihood_model=data_distr,
                                          covparms, max.iter = 50, return_all = TRUE)
    posterior_LL = calculate_posterior_VL(vecchia.approx.ll, likelihood_model=data_distr,
                                          covparms, max.iter = 50, return_all = TRUE)

    if(run_laplace) posterior_Lap= calculate_posterior_laplace(coarse_z, data_distr,
                                                               C =covfun(coarse_locs),
                                                               return_all = TRUE)


    # loop:  various m, LL, VL, Laplace + obs-lat

    #MSE/LS:  check against coarser grid?

    down_sampled = as.image(matrix(y), locs, nx = coarse_num,
                            ny = coarse_num, FUN = mean)
    coarse_y = as.vector(t(down_sampled$z))
    coarse_idx = unique(down_sampled$ind)
    coarse_locs = cbind(down_sampled$x[coarse_idx[,1]],
                        down_sampled$y[coarse_idx[,2]])

    lv_score = dposterior(coarse_y, posterior_VL)
    ll_score = dposterior(coarse_y, posterior_LL)
    laplace_score = 0
    laplace_mse = 0
    if(run_laplace){
      laplace_mse = mean((coarse_y - posterior_Lap$mean)^2)
      laplace_score = dposterior(coarse_y, posterior_Lap)
    }

    mse_vals = c(mean((coarse_y - posterior_VL$mean)^2),
                 mean((coarse_y - posterior_LL$mean)^2),
                 laplace_mse)
    ls_vals = c(lv_score, ll_score, laplace_score)

    results_i = c(smooth, range, m, seed_val, mse_vals, ls_vals, ds_factor)
    return(results_i)
  }
  return(run_LGCP)
}

smooth=.5
#### Sim 1; fixed cov params, for single evaluation of
run_scenarios  =function(smooth, use_parallel = FALSE){
  #nbrs = c(1, 5, 10, 20, 40, 60)
  nbrs = c(1, 5, 10, 20, 40, 60)
  range=  .1
  covfun <- function(locs) Matern(fields::rdist(locs),range=range,smoothness=smooth)
  n=100^2
  loc = create_locs(n, dimen = 2, dom = 1)
  Om0 <- covfun(loc)
  loc_chol = t(chol(Om0))
  LGCP_runner = get_runner(loc, loc_chol)

  scenario_table = c()
  for(iter in 1:100){
    seed_val = sample(1e6,1)
    for (m in nbrs){
      if(m==nbrs[1]){
        scenario_table = rbind(scenario_table,c(m, range, smooth, TRUE,seed_val ) )
      }else{
        scenario_table = rbind(scenario_table,c(m, range, smooth, FALSE,seed_val ) )
      }
    }
  }

  scenario_matrix = as.matrix(scenario_table)
  ########  Setup parallel
  aggregated_data = c()
  if(use_parallel){
    library(parallel)
    no_cores <- max(min(detectCores() - 1, 10), 1)
    cl <- makeCluster(no_cores)
    clusterEvalQ(cl, {source("server/importer.R")})
    clusterExport(cl, varlist = c("dposterior", "LGCP_runner","get_runner", "loc_chol", "loc"))
    aggregated_data=parApply(cl, scenario_matrix, 1, function(x) do.call(LGCP_runner, as.list(x)))
    stopCluster(cl)
  }else{
    for (i in 1:length(scenario_table[,1])){
      aggregated_data=rbind(aggregated_data, do.call(LGCP_runner, as.list(scenario_matrix[i,])))
    }
  }
  return(aggregated_data)


}




sm1_run = run_scenarios(.5, FALSE)
sm2_run = run_scenarios(1.5, FALSE)

sim_results = rbind(sm1_run, sm2_run)
header = c( "C_Smoothness", "C_Range","Neighbors","Seed",
            "MSE_VL_z", "MSE_LowRank","MSE_Laplace",
            "LS_VL_z", "LS_LowRank", "LS_Laplace", "Coarsening_Factor_from_100")

colnames(sim_results)<-header
write.csv(sim_results, file ="LGCP_sim_v3.csv", row.names = F)

t_end = Sys.time()
# emailing stuff:  for fun
require(gmailr)
require(mime)
runtime = as.double(difftime(t_end, t_start, units = "mins"))
email_subject = paste("LGCP Simulation completed")
email_text = paste("Time (minutes): ",runtime)
mime() %>%
  to("dzilber@tamu.edu") %>%
  from("dzilber@tamu.edu") %>%
  subject(email_subject) %>%
  text_body(email_text) -> text_msg
send_message(text_msg)

