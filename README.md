# GPvecchia_laplace
This repository contains code related to Vecchia-Laplace Approximations for Generalized GP's paper, excluding the package GPVecchia. 

The MODIS folder contains the MODIS satellite data sample, along with a short R script to extract the raw data and create an RData file for performing analysis.

The VL-scripts folder contains the code used to for simulation and application sections of the paper.  
  - HMCMC.R contains all code and plots that involve Hamiltonian MC posterior estimation.
  - LGCP_plot_generation.R simulates realizations of a Log Gaussian Cox process and downsamples to create gridded Poisson observations.
  - Modis_pipeline_VL.R uses the data file from the MODIS folder to perform parameter estimation via Nelder-Mead and posterior predictions using the GPVecchia package.
  - VL_manuscript_plots.R contains scripts for reproducing all remaining plots from the paper, given the saved data.
  - data_generating_functions.R and data_generating_functions2.R are both helper code files that generate samples for the simulation-based experiments.  The second version is scalable.
  - generate_scenarios.R and run_scenarios.R are the scripts that ran the 1D to 4D simulations for various sample sizes, paremeter settings, and likelihood models.  
  -parameter_estimation_grid_example.R demonstrates grid-search based parameter estimation for Poisson data.
  
The importer.R script can be run to manually load the GPVecchia package, in case the package cannot be installed via devtools::install_github("katzfuss-group/GPvecchia").  Scripts run directly if the GPVecchia library is loaded, or else should be copied into the GPVecchia folder so that GPVecchia/MODIS/ and GPVecchia/VL_scripts/ are valid subdirectories.  
