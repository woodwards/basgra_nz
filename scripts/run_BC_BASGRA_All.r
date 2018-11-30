# run all 9 sites

cat(file=stderr(), 'Starting BC_BASGRA_All.r', "\n")

# unload packages and remove variables
pkgs = names(sessionInfo()$otherPkgs)
if (FALSE){ # detach all packages
  pkgs = paste('package:', pkgs, sep = "")
  suppressMessages({
    lapply(pkgs, detach, character.only=TRUE, unload=TRUE)
  })
}
rm(list=ls()) # kills breakpoints! frees memory.
graphics.off() # closes all graphics

# load packages
suppressMessages({
  library(tidyverse)
})
# options(warn=2) # trap warnings

#### 1. INITIALISE MCMC ####
  
  # initialise BC
  source('scripts/BC_BASGRA_MCMC_init_All.r')
  # which in turn calls:
  # source('scripts/fLogL_Sivia.R')
  # source('scripts/fLogL_mm_Beta.R')
  # source('scripts/BC_BASGRA_MCMC_init_general.R')
  # which in turn calls:
  # source( sitesettings_filenames[s] ) # for each site
  # which in turn calls:
  # source('scripts/initialise_BASGRA_general.R')
  
#### 2. RUNNING THE MCMC ####
  
  cat(file=stderr(), 'Calling BC_BASGRA_BT.r', "\n")
  source('scripts/BC_BASGRA_BT.R') # new solver

#### 2b. PLOT BT RESULTS ####
  
  cat(file=stderr(), 'Calling BC_BASGRA_BT_results.r', "\n")
  source('scripts/BC_BASGRA_BT_results.R') # plot new solver results
  # which in turn calls:
  # source('scripts/plotResiduals_BT.r') # replacement functions
  
#### 3. OTHER OUTPUTS AND PLOTTING ####

  cat(file=stderr(), 'Calling BC_BASGRA_BT_other_results.r', "\n")
  source("scripts/BC_BASGRA_BT_other_results.r")
  
#### 4. CLEAN UP ####
  
  # save workspace since it takes a long time to generate
  cat(file=stderr(), 'Saving checkpoint_finished.RData', "\n")
  file_save <- 'model_outputs/checkpoint_finished.RData' 
  save.image(file_save)
  
  # finish
  cat(file=stderr(), 'Finished BC_BASGRA_All.r', "\n")
  # dyn.unload(BASGRA_DLL) 
  # installr::kill_all_Rscript_s() # kills processes left from BT parallel
  closeAllConnections()
  graphics.off()
  
  # reload workspace
  if (FALSE){
    suppressMessages(library(tidyverse))
    file_save <- 'model_outputs/checkpoint_finished.RData' 
    load(file_save)  
    BASGRA_DLL <- "model/BASGRA_WG.DLL" # use _WG version to use own PET
    dyn.load(BASGRA_DLL) 
  }

  