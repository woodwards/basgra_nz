# Bayesian calibration of BASGRA model using BayesianTools package
# https://cran.r-project.org/web/packages/BayesianTools/vignettes/BayesianTools.html
# https://www.rdocumentation.org/packages/BayesianTools/versions/0.1.3/topics/VSEM

#
suppressMessages({
  library(BayesianTools)
})

#
cat(file=stderr(), 'Calibrating BASGRA using BayesianTools package', "\n")

# parameter names
bt_names <- dplyr::if_else(parsites_BC=="1:nSites",
                    as.character(parname_BC),
                    paste(parname_BC, "(", parsites_BC, ")", sep=""))
cat(file=stderr(), paste('Adjustable parameters =', length(bt_names)), "\n")

# likelihood function (work with scaled parameters)
par <- rep(1, length(bt_names)) # for testing
s <- 1 # for testing
bt_likelihood <- function(par){
  # profvis::profvis({for (i in 1:1000){
    # use loop from BC_BASGRA_MCMC.R  
    candidatepValues_BC   <- par * sc
    for (s in 1:nSites) {
      params         <- list_params        [[s]] # get site parameters initial values (in parameters.txt)
      matrix_weather <- list_matrix_weather[[s]] # get site weather
      days_harvest   <- list_days_harvest  [[s]] # get site harvest
      NDAYS          <- list_NDAYS         [[s]] # get site NDAYS
      # ip_BC_site[[s]] = indicies of model parameters being changed (in parameters.txt)
      # icol_pChain_site[[s]] = indices of calibration parameters being used (in parameters_BC.txt)
      params[ ip_BC_site[[s]] ] <- candidatepValues_BC[ icol_pChain_site[[s]] ]
      list_output[[s]] <- run_model(params,matrix_weather,days_harvest,NDAYS,NOUT,list_output[[s]])
    }
    # use functions from BC_BASGRA_init_general.R
    logL1 <- calc_sum_logL( list_output ) # likelihood function in BC_BASGRA_init_general.R
  # }}) # profvis
  # return likelihood
  return(logL1)
}
bt_likelihood(rep(1, np_BC))

# construct priors (scaled parameter space)
bt_prior <- createBetaPrior(aa, bb, scparmin_BC[1:np_BC], scparmax_BC[1:np_BC])

# construct setup
bt_setup <- createBayesianSetup(likelihood=bt_likelihood, 
                                prior=bt_prior, 
                                parallel=FALSE, # number of cores
                                names=bt_names)

# construct settings (note: DREAMzs has startValue=3 internal chains by default)
# Possibly DREAMzs has limited capability to use parallel cores
nInternal   <- 3 # internal chains for DREAMzs
bt_settings <- list(iterations=nChain/nChains,
                    startValue=nInternal, 
                    # startValue=matrix(1,nInternal,np_BC), 
                    nrChains=nChains, 
                    # burnin=0, # because can't analyse convergence if we discard burnin
                    burnin=nBurnin/nChains+nChains, # to give correct number of samples
                    parallel=FALSE, # can overrule parallel=FALSE in BayesianSetup
                    consoleUpdates=1000,
                    message=TRUE) 

# run BT until stopping conditions met (these can be changed in the file BC_BASGRA_BT_stop.csv)
bt_chains <- nInternal * nChains 
bt_pars <- length(bt_names)
bt_conv <- NA
bt_time <- 0
# profile <- profvis({
# source("scripts/mcmcRun.R") # for profiling
# source("scripts/mcmcDREAMzs.R") # for profiling
# source("scripts/mcmcDREAM_helperFunctions.R") # for profiling
repeat{
  
  # run Bayesian Tools
  if (is.na(bt_conv)){ 
    # first run
    print(elapsed <- system.time({
      bt_out <- runMCMC(bayesianSetup = bt_setup, 
                        sampler = "DREAMzs", 
                        settings = bt_settings)
      cat(file=stderr(), " ", "\n")
    }))
  } else {
    # continuation
    if (nBurnin==0){
      cat(file=stderr(), paste("Greater than", conv_target, "and less than", conv_minutes, "so continuing..."), "\n")
    } else {
      cat(file=stderr(), "Restart doesn't work with nBurnin>0 so stopping...\n")
      # stop()
    }
    print(elapsed <- system.time({
      bt_out <- runMCMC(bayesianSetup = bt_out)
      cat(file=stderr(), " ", "\n")
    }))
  }
  
  # assess convergence  
  bt_length <- dim(bt_out[[1]]$chain[[1]])[[1]]
  cat(file=stderr(), paste("Total chains =", bt_chains), "\n")
  cat(file=stderr(), paste("Total samples per chain =", bt_length), "\n")
  # bt_conv <- gelmanDiagnostics(bt_out)$mpsrf 
  # cat(file=stderr(), paste("Overall convergence (mpsrf) =", round(bt_conv,3)), "\n")
  bt_conv <- max(gelmanDiagnostics(bt_out)$psrf[,1])
  cat(file=stderr(), paste("Convergence max(psf) =", round(bt_conv,3)), "\n")
  # bt_time <- sum(sapply(bt_out, function(x) x$settings$runtime[3]/60, simplify=TRUE))
  # cat(file=stderr(), paste("Total time =", round(bt_time,2), "minutes"), "\n")
  bt_time <- bt_time + elapsed[[3]]/60
  cat(file=stderr(), paste("Elapsed time =", round(bt_time,2), "minutes"), "\n")
  
  # read stopping criterion from file (allows us to change it on the fly)
  conv_criteria <- read.csv("scripts/BC_BASGRA_BT_stop.csv", header=FALSE) 
  conv_target <- conv_criteria[1,1]
  conv_minutes <- conv_criteria[2,1]
  if ((bt_conv <= conv_target) || (bt_time >= conv_minutes)){
    break
  }
}
# }) # end profvis
# htmlwidgets::saveWidget(profile, "profile.html")
# browseURL("profile.html")

# stop parallel if needed
stopParallel(bayesianSetup=bt_setup)

# report convergence
cat(file=stderr(), paste("Convergence of individual parameters (psf)"), "\n")
psf <- gelmanDiagnostics(bt_out)$psrf[,1]
print(round(psf,3))
cat(file=stderr(), paste("Convergence of worst parameter (psf)"), "\n")
print(round(psf[which.max(psf)],3))
cat(file=stderr(), paste("Multivariate convergence (mpsrf) ="), "\n")
print(round(bt_conv,3))

# memory management
cat(file=stderr(), 'Saving checkpoint after BASGRA calibration', "\n")
file_save <- paste(scenario, "/checkpoint_after_calibration.RData", sep="")
save.image(file=file_save)
rm(list=setdiff(ls(), c("scenario", "scenarios"))) # avoid memory overflow

