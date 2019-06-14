rm(list = ls())


setwd("../GPVecchia/")

###  load GPvecchia package manually, if not using library(GPVecchia)
for (nm in list.files("R",pattern = "\\.[RrSsQq]$")) {
  cat(nm,":"); source(file.path("R",nm)); cat("\n")
}


Rcpp::sourceCpp('src/U_NZentries.cpp')
Rcpp::sourceCpp('src/MaxMin.cpp')
Rcpp::sourceCpp('src/fastTree.cpp')

library(GpGp); library(Matrix); library(RcppParallel)
library(parallel); library(sparseinv); library(fields)

# for running the scenario simulations, uncomment
#source("server/run_scenario.R")
#source("server/data_generating_functions2.R")
