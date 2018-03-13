setwd("/Users/danielzilber/Desktop/Katzfuss/Transfer_CPP/")
source("VL_method_cpp.R")
source("run_scenario_cpp.R")
source("data_generating_functions_cpp.R")

##########################################################################
##################### Compare models:  MSE and Log score  ########################
##########################################################################




#loop over different domains and sample size


#vary smoothness and m=1,2,3:  how does approximation accuracy change with m
## y: approximation , x =m (SGV )  compare to pure lapalce
# for 2d:  smoothness 1.5 and .5
d_vals = c(1)  # domain, [0,1]
s_vals = c(625)
seed_vals = c(120:221)# seq(141, 166)
smoothness_vals = c(.5, 1.5)
range_vals = .05
#nbrs = c(1, 2, 3, 4, 5, 10, 20, 50)#, 100, 250, 500)
# dimen = 1  # dimension, 1 or 2
models_tested = c(1,2,3,4)


#2d Scenarios
nbrs = c(1, 2, 3, 4, 5, 10, 20, 50, 100, 250, 500)
dimen = 2

scenario_table = c()

# Compare log likelihood of true y for posterior density of approx p(y | z)
for (seed_r in seed_vals){
  for (domn in d_vals){
    for (samp_size in s_vals){
      #model_names = c("logistic", "pois", "gauss", "gamma")[mod_type]
      for( neighbors in nbrs ){
        for(smth in smoothness_vals){
          for(rnge in range_vals){
            for(mod_type in models_tested){
              scenario_table= rbind(scenario_table, c(seed_r, domn, dimen, samp_size, neighbors, smth, mod_type, rnge, TRUE))
            }
          }
        }
      }
    }
  }
}


filename = "2D_sim.csv"
scenario_runner = create_scenario_tester(filename)
scen_params =c("seed_r", "domn", "dimen", "samp_size", "neighbors", "smth", "mod_type", "rnge", "show_output")
#for (i in 1:length(scenario_table[,1])){
#  params = setNames( as.list(scenario_table[i,]), c(scen_params))
#  do.call(scenario_runner, params)
#}


library(parallel)

no_cores <- min(detectCores() - 1, 10)
cl <- makeCluster(no_cores, type="FORK")
start = Sys.time()
parApply(cl, scenario_table, 1 , function(x) do.call(scenario_runner, as.list(x)))
time_para = Sys.time()-start
stopCluster(cl)

#start = Sys.time()
#apply(scenario_table, 1 , function(x) do.call(scenario_runner, as.list(x)) )
#time_serial = Sys.time()-start
