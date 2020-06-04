#!/usr/bin/env Rscript
source("VL_scripts/importer.R")
input_settings = commandArgs(trailingOnly=TRUE)

# This code is called from terminal to generate a collection
# of scenarios to be tested with simulation.
# Input:  type of scenario, such as 2D gaussian with various sample sizes
# Output:  simulation results for the generated scenarios


print(length(input_settings))
print(input_settings)
# test if there is at least one argument: if not, return an error
if (length(input_settings)==0) {
  stop("Specify an output file name, parallelization bool, and test.\n", call.=FALSE)
} else if (length(input_settings)==1) {
  # default output file
  input_settings[2] = "F"  # "F" or "T"
  input_settings[3] = "none"  # "1D" or "2D" or "time" or anything else (runs short test)
}

filename = input_settings[1] #filename  = "test"
use_parallel = input_settings[2]
test_type = input_settings[3]# "1D" "2D" "time", misc


paste("Running simulation with output to",filename,"and parallel set to",use_parallel, "and test type of",test_type)

##########################################################################
##################### Compare models:  MSE and Log score  ########################
##########################################################################


########  Setup parallel
if(use_parallel=="T"){
  library(parallel)
  no_cores <- max(min(detectCores() - 1, 10), 1)
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {source("server/importer.R")})
}

########  Setup Simulations


#loop over different domains and sample size

# default: for fast test
d_vals = c(1)  # domain, [0,1]
s_vals = c(100)
seed_vals = rep(-1,10)# seq(141, 166)
smoothness_vals = .5
range_vals = .05
nbrs = c(5)
dimen = 2  # dimension, 1 or 2
models_tested = c(1,2,3,4)
run_laplace = TRUE



#vary smoothness and m=1,2,3:  how does approximation accuracy change with m
## y: approximation , x =m (SGV )  compare to pure lapalce
# for 2d:  smoothness 1.5 and .5
if(test_type == "1D" | test_type == "2D"){
  d_vals = c(1)  # domain, [0,1]
  s_vals = c(2500)
  seed_vals = rep(-1, 100)#
  smoothness_vals = c(.5, 1.5) #seq(2.5, 2.8, length.out = 25)
  nbrs = c(1, 3, 5, 10, 20, 40, 60)
  dimen = ifelse(test_type == "1D",1, 2)
}
if(test_type == "3D"){
  s_vals = c(13^3)
  range_vals = .1
  seed_vals = rep(-1, 50)#
  smoothness_vals = c(.5, 1.5) #seq(2.5, 2.8, length.out = 25)
  nbrs = c(1, 3, 5, 10, 20, 40)
  dimen = 3
}
if( test_type == "4D" ){
  s_vals = c(7^4)
  range_vals = .2
  seed_vals = rep(-1, 50)#
  smoothness_vals = c(.5, 1.5) #seq(2.5, 2.8, length.out = 25)
  nbrs = c(1, 3, 5, 10, 20, 40)
  dimen = 4
}

if(test_type == "time"){
  d_vals = c(10)  # domain, [0,1]
  s_vals = c(250, 500, 1000, 2000, 4000, 8000, 16000)
  models_tested = 1#c(2,3,4)
  nbrs = c(10)
  dimen=1
}
if(test_type == "time2"){
  d_vals = c(10)  # domain, [0,1]
  s_vals = c(16, 25, 33, 45, 64, 90, 130, 180, 273, 387, 550)^2
  models_tested = 1#c(2,3,4)
  range_vals = .3
  nbrs = c(10)
  dimen=2
  run_laplace = FALSE
}
if(test_type == "2D_sm"){
  d_vals = c(1)  # domain, [0,1]
  s_vals = c(2500)
  seed_vals = rep(-1, 100)#
  smoothness_vals = c(.5, 1.5) #seq(2.5, 2.8, length.out = 25)
  range_vals = .2
  nbrs = c(1, 3, 5, 10, 20, 40, 60)
  dimen = 2
}






scenario_table = c()

#2d Scenarios
#nbrs = c(1, 2, 3, 4, 5, 10, 20, 50, 100, 250, 500)
#dimen = 2


# Compare log likelihood of true y for posterior density of approx p(y | z)
for (seed_r in seed_vals){
  for (domn in d_vals){
    for (samp_size in s_vals){
      #model_names = c("logistic", "pois", "gauss", "gamma")[mod_type]
      for( neighbors in nbrs ){
        for(smth in smoothness_vals){
          for(rnge in range_vals){
            for(mod_type in models_tested){
              scenario_table= rbind(scenario_table, c(seed_r, domn, dimen, samp_size, neighbors, smth, mod_type, rnge, TRUE, run_laplace))
            }
          }
        }
      }
    }
  }
}

aggregated_data = c()



header = c("Mod", "Domain", "Dimen", "Sample", "C_Smoothness", "C_Range","Neighbors","Seed_off",
           "MSE_Laplace", "MSE_VL", "MSE_VL_z",  "MSE_LowRank",
           "LS_Laplace", "LS_VL",  "LS_VL_z", "LS_LowRank",
           "Time_Laplace", "Time_VL", "Time_VL_z",  "Time_LowRank",
           "Iter_Laplace", "Iter_VL", "Iter_VL_z",  "Iter_LowRank")

scen_params =c("seed_r", "domn", "dimen", "samp_size", "neighbors", "smth", "mod_type", "rnge", "show_output", "run_algo")


source("VL_scripts/create_scenario_tester.R")

t_start = Sys.time()

##  Non-parallel
if(use_parallel == "F"){
  scenario_runner = create_scenario_tester(header, filename)# filename= "delete_me.csv"
  for (i in 1:length(scenario_table[,1])){
    params = setNames( as.list(scenario_table[i,]), c(scen_params))
    aggregated_data=rbind(aggregated_data, do.call(scenario_runner, params))
  }
}
## Parallel
if(use_parallel=="T"){
  scenario_runner = create_scenario_tester(header, NA) # cant write during parallel
  clusterExport(cl, varlist = c("scenario_runner"))
  aggregated_data=parApply(cl, scenario_table, 1 , function(x) do.call(scenario_runner, as.list(x)))
  stopCluster(cl)
  aggregated_data = t(aggregated_data)
}

t_end = Sys.time()

colnames(aggregated_data)<-header
write.csv(aggregated_data, file = paste("alt_",filename, sep = ""), row.names = F)


# emailing stuff:  for fun
require(gmailr)
require(mime)
parallel_time = as.double(difftime(t_end, t_start, units = "mins"))
email_subject = paste("Simulation",filename,"completed")
email_text = paste("Time (minutes): ",parallel_time)
mime() %>%
  to("your_email@school.edu") %>%
  from("your_email@school.edu") %>%
  subject(email_subject) %>%
  text_body(email_text) -> text_msg
send_message(text_msg)


