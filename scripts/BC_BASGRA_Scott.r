## 1. INITIALISE MCMC ##
   source('scripts/BC_BASGRA_MCMC_init_Scott.r')

## 2. RUNNING THE MCMC ##
   source('scripts/BC_BASGRA_MCMC.R')

## 3. PRINTING AND PLOTTING ##
   source('scripts/BC_export_parModes.R')
#   source('scripts/BC_plot_parameters_traceplots.R') # can be slow!
   source('scripts/BC_plot_parameters_priorbeta_histograms.R')
   source('scripts/BC_plot_outputs_data.R')