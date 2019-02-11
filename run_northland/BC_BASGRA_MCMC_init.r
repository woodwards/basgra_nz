## BC_BASGRA_MCMC_init.R ##
   cat(file=stderr(), "Starting BC_BASGRA_MCMC_init.r", "\n")

## MCMC chain length (total number of iterations across all chains)
   # https://stats.stackexchange.com/questions/266749/mcmc-convergence
   # suggest nBurnin up to 50% of nChain
   nChains       <- 3 # independent chains for repeatability testing
   nInternal     <- 3 # internal chains in DREAMzs
   nBurnin       <- as.integer(00000 * nChains * nInternal) 
   nSampling     <- as.integer(10000 * nChains * nInternal)
   nChain        <- nBurnin + nSampling # total samples 

## FILE FOR PRIOR PARAMETER DISTRIBUTION
   # this file provides parameters for all the sites
   # shape <- 0 # shape parameter (0=noninformative, 4=previous method) = now included in par file
   
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
   source("scripts/fLogL_Sivia.R")
   source("scripts/fLogL_mm_Beta.R")
   
## SETTINGS FOR THE DIFFERENT CALIBRATION SITES (at least one site)
   sitesettings_filenames <- c()
   sitedata_filenames     <- c()
   file_prior    <- paste(scenario, "/parameters_BC.txt", sep="")

   sitesettings_filenames <- c(sitesettings_filenames,
                               paste(scenario, "/initialise_BASGRA_Northland_0.r", sep=""))
   sitedata_filenames     <- c(sitedata_filenames,
                               paste(scenario, "/data_calibration_Northland_0.txt", sep=""))
   # 
   # sitesettings_filenames <- c(sitesettings_filenames,
   #                             paste(scenario, "/initialise_BASGRA_Scott_2.r", sep=""),
   #                             paste(scenario, "/initialise_BASGRA_Scott_3.r", sep=""),
   #                             paste(scenario, "/initialise_BASGRA_Scott_4.r", sep=""))
   # sitedata_filenames     <- c(sitedata_filenames,
   #                             paste(scenario, "/data_calibration_Scott_2.txt", sep=""),
   #                             paste(scenario, "/data_calibration_Scott_3.txt", sep=""),
   #                             paste(scenario, "/data_calibration_Scott_4.txt", sep=""))
   
   # sitesettings_filenames <- c(sitesettings_filenames,
   #                             paste(scenario, "/initialise_BASGRA_Lincoln_2.r", sep=""),
   #                             paste(scenario, "/initialise_BASGRA_Lincoln_3.r", sep=""),
   #                             paste(scenario, "/initialise_BASGRA_Lincoln_5.r", sep=""))
   # sitedata_filenames     <- c(sitedata_filenames,
   #                             paste(scenario, "/data_calibration_Lincoln_2.txt", sep=""),
   #                             paste(scenario, "/data_calibration_Lincoln_3.txt", sep=""),
   #                             paste(scenario, "/data_calibration_Lincoln_5.txt", sep=""))
   
   nSites                 <- length(sitedata_filenames)
   sitelist               <- list() ; length(sitelist) <- nSites
   
   # additional outputs to plot
   extraOutputs <- c("LAI", "TSIZE", "CLVD", "CRT", "DM", "RES", "TRANRF")
   
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
   # fPropTuning   <- 0.05 # This factor is used to modify Gelman"s suggested average step length
   #                       # (2.38^2 / np_BC) which seems too big
   fPropTuning   <- 0.1  # Simon (not used in BT)
   
## GENERAL INITIALISATION FOR MCMC
   source("scripts/BC_BASGRA_MCMC_init_general.R")

   #
   cat(file=stderr(), "Finished BC_BASGRA_MCMC_init.r", "\n")
   