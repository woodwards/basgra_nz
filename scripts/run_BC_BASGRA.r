# run_BC_BASGRA
cat(file = stderr(), "Starting run_BC_BASGRA.r", "\n")

# closes all graphics
graphics.off()

# load packages
suppressMessages({
  library(tidyverse)
})

# Load DLL here (parallel version uses BASGRA package instead)
BASGRA_DLL <- "model/BASGRA_WG.DLL" # use _WG version to use own PET
if (.Machine$sizeof.pointer == 4) {
  print("32-bit R detected")
  print("Warning: Might have memory problems")
  dyn.load(BASGRA_DLL)
} else if (.Machine$sizeof.pointer == 8) {
  print("64-bit R detected")
  # print(paste("Can't load", BASGRA_DLL))
}

#### user options ####
# now able to calibrate and run scenarios in same batch
run_parallel <- TRUE
fixpars <- FALSE
fixpar <- c()
single_year <- 2013
base_s <- 1
scenarios <- c("run_mean", "run_mean/scenario_northland", "run_mean/scenario_scott")
parameters <- c("run_mean", "run_mean", "run_mean")
fit_mcmcs <- c(FALSE, FALSE, FALSE)
plot_varsets <- c(1, 2, 2)

#### point to scenario directory ####
is <- 1:3
for (i in is) {

  #
  scenario <- scenarios[i]
  fit_mcmc <- fit_mcmcs[i]
  plot_varset <- plot_varsets[i]
  parameter_location <- parameters[i]
  cat(file = stderr(), "\nStarting scenario", scenario, "\n")
  if (!fit_mcmc){
    cat(file = stderr(), "Using parameters from", parameter_location, "\n")
  }
  
  #### 1. INITIALISE MCMC ####

  # random seed
  set.seed(123)

  # initialise BC
  file_name <- paste(scenario, "/BC_BASGRA_MCMC_init.r", sep = "")
  cat(file = stderr(), "Calling", file_name, "\n")
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
  cat(file = stderr(), "Calling", file_name, "\n")
  source(file_name)

  #### 2b. PLOT BT RESULTS ####

  file_name <- "scripts/BC_BASGRA_BT_results.R"
  cat(file = stderr(), "Calling", file_name, "\n")
  source(file_name)
  # which in turn calls:
  # source("scripts/plotResiduals_BT.r") # replacement functions

  #### 3. OTHER OUTPUTS AND PLOTTING ####

  if (FALSE) {
    file_name <- "scripts/BC_BASGRA_BT_other_results.R"
    cat(file = stderr(), "Calling", file_name, "\n")
    source(file_name)
  }

  #### 4. CLEAN UP ####

  # save workspace since it takes a long time to generate
  cat(file = stderr(), "Saving checkpoint_finished.RData", "\n")
  file_save <- paste(scenario, "/checkpoint_finished.RData", sep = "")
  save.image(file_save)
  
  #
  cat(file = stderr(), "Finished scenario", scenario, "\n")
  
  # reload workspace
  if (FALSE) {
    suppressMessages(library(tidyverse))
    file_save <- paste(scenario, "/checkpoint_finished.RData", sep = "")
    load(file_save)
    # BASGRA_DLL <- "model/BASGRA_WG.DLL" # use _WG version to use own PET
    # dyn.load(BASGRA_DLL)
  }
} # next scenario

cat(file = stderr(), "Finished BC_BASGRA.r", "\n")
if (FALSE) {
  # dyn.unload(BASGRA_DLL)
  installr::kill_all_Rscript_s() # kills processes left from BT parallel
  closeAllConnections()
  graphics.off()
}

# open amusing browser window to tell me it's finished :)
browseURL("https://www.youtube.com/watch?v=QH2-TGUlwu4")
