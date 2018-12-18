# Bayesian calibration of BASGRA model using BayesianTools package
# https://cran.r-project.org/web/packages/BayesianTools/vignettes/BayesianTools.html
# https://www.rdocumentation.org/packages/BayesianTools/versions/0.1.3/topics/VSEM

#
suppressMessages({
  library(BayesianTools)
  library(coda)
})

#
cat(file=stderr(), 'Calibrating BASGRA using BayesianTools package', "\n")

# parameter names
bt_names <- if_else(parsites_BC=="1:nSites",
                    as.character(parname_BC),
                    paste(parname_BC, "(", parsites_BC, ")", sep=""))
cat(file=stderr(), paste('Adjustable parameters =', length(bt_names)), "\n")

# likelihood function (work with scaled parameters)
bt_likelihood <- function(par){
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
    output                    <- run_model(params,matrix_weather,days_harvest,NDAYS)
    list_output[[s]]          <- output
  }
  # use functions from BC_BASGRA_init_general.R
  logL1 <- calc_sum_logL( list_output ) # likelihood function in BC_BASGRA_init_general.R
  # return likelihood
  return(logL1)
}

# construct priors (scaled parameter space)
bt_prior <- createBetaPrior(aa, bb, scparmin_BC[1:np_BC], scparmax_BC[1:np_BC])

# construct setup
bt_setup <- createBayesianSetup(likelihood=bt_likelihood, 
                                prior=bt_prior, 
                                parallel=TRUE,
                                parallelOptions=list(dlls=list(BASGRA_DLL)),
                                names=bt_names)

# construct settings (note: DREAMzs has startValue=3 internal chains by default)
# Possibly DREAMzs has limited capability to use parallel cores
nInternal   <- 3 # internal chains for DREAMzs
bt_settings <- list(startValue=nInternal, 
                    iterations=nChain/nChains, 
                    nrChains=nChains, 
                    # burnin=0, # because can't analyse convergence if we discard burnin
                    burnin=nBurnin/nChains+nChains, # to give correct number of samples
                    parallel=TRUE,
                    message=TRUE) 

# run BT
bt_out <- runMCMC(bayesianSetup = bt_setup, 
                  sampler = "DREAMzs", 
                  settings = bt_settings)
cat(file=stderr(), " ", "\n")
bt_chains <- nInternal * nChains 
bt_length <- dim(bt_out[[1]]$chain[[1]])[[1]]
bt_pars <- length(bt_names)
bt_conv <- gelmanDiagnostics(bt_out)$mpsrf
cat(file=stderr(), paste("Total chains =", bt_chains), "\n")
cat(file=stderr(), paste("Total samples per chain =", bt_length), "\n")
cat(file=stderr(), paste("Overall convergence (mpsrf) =", round(bt_conv,3)), "\n")

# rerun BT until target mpsrf achieved (provided run time remains reasonable)
# DOESN'T REALLY SEEM TO CONVERGE, AND OFTEN CAUSES A CRASH, ALSO CAN'T HANDLE BURNIN
# target_mpsrf <- 2.0
# max_minutes <- 60
# if (FALSE){ 
#   while ((gelmanDiagnostics(bt_out)$mpsrf > target_mpsrf)&
#          (bt_out[[1]]$settings$runtime[3] < max_minutes*60)){
#     # restart
#   	if (nBurnin==0){
#   	  cat(file=stderr(), paste("Greater than", target_mpsrf, "so continuing..."), "\n")
#   	} else {
#   		stop("restart doesn't work with nBurnin>0 due to a bug")
#   	}
#     bt_out <- runMCMC(bayesianSetup = bt_out) 
#     cat(file=stderr(), " ", "\n")
#     # bt_chains <- nInternal * nChains
#     bt_length <- dim(bt_out[[1]]$chain[[1]])[[1]]
#     # bt_pars <- length(bt_names)
#     bt_conv <- gelmanDiagnostics(bt_out)$mpsrf
#     # cat(file=stderr(), paste("Total chains =", bt_chains), "\n")
#     cat(file=stderr(), paste("Total samples per chain =", bt_length), "\n")
#     cat(file=stderr(), paste("Overall convergence (mpsrf) =", round(bt_conv,3)), "\n")
#   }
# }

# stop parallel
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


