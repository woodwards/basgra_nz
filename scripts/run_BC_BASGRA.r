# BC_BASGRA run 
cat(file=stderr(), "Starting BC_BASGRA.r", "\n")

# unload packages and remove variables
pkgs = names(sessionInfo()$otherPkgs)
if (FALSE){ # detach all packages
  pkgs = paste("package:", pkgs, sep = "")
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

#### point to scenario directory ####
scenarios <- c("run_mean")
scenario <- scenarios[[1]]
for (scenario in scenarios){
  
#### 1. INITIALISE MCMC ####

  set.seed(1234)
  
  # initialise BC
  file_name <- paste(scenario, "/BC_BASGRA_MCMC_init.r", sep="")
  cat(file=stderr(), "Calling", file_name, "\n")
  source(file_name)
  # which in turn calls:
  # source("scripts/fLogL_Sivia.R")
  # source("scripts/fLogL_mm_Beta.R")
  # source("scripts/BC_BASGRA_MCMC_init_general.R")
  # which in turn calls:
  # source( sitesettings_filenames[s] ) # for each site
  # which in turn calls:
  # source("scripts/initialise_BASGRA_general.R")
  
#### 2. RUNNING THE MCMC ####
  
  file_name <- "scripts/BC_BASGRA_BT.R"
  cat(file=stderr(), "Calling", file_name, "\n")
  source(file_name)

#### 2b. PLOT BT RESULTS ####
  
  file_name <- "scripts/BC_BASGRA_BT_results.R"
  cat(file=stderr(), "Calling", file_name, "\n")
  source(file_name)
  # which in turn calls:
  # source("scripts/plotResiduals_BT.r") # replacement functions
  
#### 3. OTHER OUTPUTS AND PLOTTING ####

  file_name <- "scripts/BC_BASGRA_BT_other_results.R"
  cat(file=stderr(), "Calling", file_name, "\n")
  source(file_name)

#### 4. CLEAN UP ####
  
  # save workspace since it takes a long time to generate
  cat(file=stderr(), "Saving checkpoint_finished.RData", "\n")
  file_save <- paste(scenario, "/checkpoint_finished.RData", sep="") 
  save.image(file_save)
  
  # reload workspace
  if (FALSE){
    suppressMessages(library(tidyverse))
    file_save <- paste(scenario, "/checkpoint_finished.RData", sep="") 
    load(file_save)  
    BASGRA_DLL <- "model/BASGRA_WG.DLL" # use _WG version to use own PET
    dyn.load(BASGRA_DLL) 
  }

} # next scenario

cat(file=stderr(), "Finished BC_BASGRA.r", "\n")
if (FALSE){
  dyn.unload(BASGRA_DLL)
  installr::kill_all_Rscript_s() # kills processes left from BT parallel
  closeAllConnections()
  graphics.off()
}

# open silly browser window to tell me it's finished :)
browseURL('https://www.youtube.com/watch?v=QH2-TGUlwu4')

