## 1. INITIALISE MCMC ##
  cat(file=stderr(), 'Starting BC_BASGRA_Scott.r', "\n")
  #rm(list=ls()) # kills breakpoints!
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
  
  #
  cat(file=stderr(), 'Finished BC_BASGRA_Scott.r', "\n")
  