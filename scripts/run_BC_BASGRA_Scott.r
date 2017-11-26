## 1. INITIALISE MCMC ##
  cat(file=stderr(), 'Starting BC_BASGRA_Scott.r', "\n")

  # unload packages and remove variables
  pkgs = names(sessionInfo()$otherPkgs)
  if (!is.null(pkgs)){
    pkgs = paste('package:', pkgs, sep = "")
    lapply(pkgs, detach, character.only=TRUE, unload=TRUE)
  }
  rm(list=ls()) # kills breakpoints!
  
  # point to parameters
  file_params    <- 'model_inputs/parameters_Scott.txt' # can contain multiple columns
  parcol       <- 1 # which one are we going to use? (row names are ignored)
  
  # initialise BC
  source('scripts/BC_BASGRA_MCMC_init_Scott.r')

## 2. RUNNING THE MCMC ##
  cat(file=stderr(), 'Calling  BC_BASGRA_MCMC.r', "\n")
  source('scripts/BC_BASGRA_MCMC.R')

## 3. PRINTING AND PLOTTING ##
  cat(file=stderr(), 'Calling BC output scripts', "\n")
  source('scripts/BC_export_parModes.R')
  source('scripts/BC_plot_parameters_traceplots.R') # can be slow!
  source('scripts/BC_plot_parameters_priorbeta_histograms.R')
  source('scripts/BC_plot_outputs_data.R')
  
  # run with MAP parameters
  cat(file=stderr(), 'Running BASGRA_Scott with BC MAP parameters', "\n")
  file_params    <- 'model_outputs/BASGRA_parModes.txt' 
  parcol         <- 2 
  source('scripts/initialise_BASGRA_Scott.r')
  cat(file=stderr(), 'Calling run_model()', "\n")
  output <- run_model()
  cat(file=stderr(), 'Calling export_output()', "\n")
  export_output()

  #
  cat(file=stderr(), 'Finished BC_BASGRA_Scott.r', "\n")
  