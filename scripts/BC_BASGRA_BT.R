# Bayesian calibration of BASGRA model using BayesianTools package
# https://cran.r-project.org/web/packages/BayesianTools/vignettes/BayesianTools.html
library(BayesianTools)

#
cat(file=stderr(), 'Calibrating BASGRA using BayesianTools package', "\n")

# construct names
bt_names <- parname_BC

# construct likelihood
bt_likelihood <- function(par, sum=TRUE){
  # use loop from BC_BASGRA_MCMC.R  
  candidatepValues_BC   <- par
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
  if (sum==TRUE){
    logL1 <- calc_sum_logL( list_output ) # likelihood function in BC_BASGRA_init_general.R
    return(logL1)
  } else {
    stop() # FIXME currently don't have this information
  }
}

# construct priors
bt_prior <- createBetaPrior(aa, bb, scparmin_BC[1:np_BC], scparmax_BC[1:np_BC])

# construct setup
bt_setup <- createBayesianSetup(likelihood=bt_likelihood, 
                                prior=bt_prior, 
                                names=bt_names, 
                                parallel=TRUE, 
                                parallelOptions=list(dlls=list(BASGRA_DLL))
                                )

# construct settings (note: DREAMzs has startValue=3 internal chains by default)
bt_settings <- list(startValue=nChains, iterations=nChain, nrChains=1, burnin=0)

# run BT
bt_out <- runMCMC(bayesianSetup=bt_setup, 
                  sampler = "DREAMzs", 
                  settings=bt_settings)

# get samples
pChain <- getSample(bt_out)

# display results
summary(bt_out)
# plot(bt_out)
# marginalPlot(bt_out)
# # correlationPlot(bt_out) # too many parameters
# gelmanDiagnostics(bt_out, plot=TRUE)


