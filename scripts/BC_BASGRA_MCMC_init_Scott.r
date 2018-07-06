## BC_BASGRA_MCMC_init_xxx.R ##
   cat(file=stderr(), 'Starting BC_BASGRA_MCMC_init_Scott.r', "\n")

## MCMC chain length (total number of iterations across all chains)
   nChains       <- 3 # chains for repeatability testing
   nBurnin       <- as.integer(12000 * nChains) # need to check convergence
   nSampling     <- as.integer(12000 * nChains)
   nChain        <- nBurnin + nSampling

## FILE FOR PRIOR PARAMETER DISTRIBUTION
   # this file provides parameters for all the sites
   # shape <- 0 # shape parameter (0=noninformative, 4=previous method) = now included in par file
   file_prior    <- "model_inputs/parameters_BC_Scott.txt"
   
   # parameters file has 5 columns
   # name  minimum  mode  maximum  (beta distribution)  unknown
   # require
   # .	DLMXGE > DAYLB 
   # .	TOPTGE > TBASE 
   # .	FSMAX has a theoretical upper limit < 1. 
   # .	HAGERE <= 1 
   # .	SHAPE <= 1 
   # .	SLAMAX > SLAMIN 
   # .	TRANCO may have physical limits [a,b] where a>0 and b<infinity. 
   # .	YG < 1 because it is the Growth Yield, the fraction of C allocated to growth that 
   #      actually ends up in new biomass, with the remainder being lost to growth respiration. 
   
## LIKELIHOOD FUNCTION ##
   source('scripts/fLogL_Sivia.R')
   source('scripts/fLogL_mm_Beta.R')
   
## SETTINGS FOR THE DIFFERENT CALIBRATION SITES (at least one site)
   sitesettings_filenames <- c("scripts/initialise_BASGRA_Scott_2.r",
                               "scripts/initialise_BASGRA_Scott_3.r",
                               "scripts/initialise_BASGRA_Scott_4.r"
                               )
   sitedata_filenames     <- c("model_inputs/data_calibration_Scott_2.txt",
                               "model_inputs/data_calibration_Scott_3.txt",
                               "model_inputs/data_calibration_Scott_4.txt"
                               )
   nSites                 <- length(sitedata_filenames)
   sitelist               <- list() ; length(sitelist) <- nSites

   # Specify data uncertainties (the max of: cv for relative uncertainty, sd for absolute)   
   # These are used in BC_BASGRA_MCMC_init_general.r to set the data uncertainites
   # Data uncertainties are now specified in the data file
   # cv_default    <- rep( 0.2 , nSites ) # Simon change from 0.5 to 0.2
   # cv_DM         <- rep( 0.05, nSites ) ; sd_DM_min     <- rep(   0, nSites )
   # cv_LAI        <- rep( 0.1 , nSites ) ; sd_LAI_min    <- rep(   0, nSites )
   # cv_TILTOT     <- rep( 0.0 , nSites ) ; sd_TILTOT_min <- rep( 200, nSites )
   # cv_YIELD      <- rep( 0.05, nSites ) ; sd_YIELD_min  <- rep(   0, nSites )
   # sd_LT50       <- rep( 5   , nSites )
   # sd_CST        <- rep( 5   , nSites ) # need to add this because has values of 0
   # sd_CLV        <- rep( 5  , nSites ) 
   # sd_WCL        <- rep(  2.5  , nSites ) 
   # cv_mm_default <- rep( 0.2 , nSites ) # “minimum and maximum” variables
   # cv_mm_FRTILG  <- rep( 0.2 , nSites )
   
## PROPOSAL TUNING FACTOR  
   # fPropTuning   <- 0.05 # This factor is used to modify Gelman's suggested average step length
   #                       # (2.38^2 / np_BC) which seems too big
   fPropTuning   <- 0.1  # Simon (not used in BT)
   
## GENERAL INITIALISATION FOR MCMC
   source('scripts/BC_BASGRA_MCMC_init_general.R')

   #
   cat(file=stderr(), 'Finished BC_BASGRA_MCMC_init_Scott.r', "\n")
   