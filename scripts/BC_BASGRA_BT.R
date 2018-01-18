# Bayesian calibration of BASGRA model using BayesianTools package
# https://cran.r-project.org/web/packages/BayesianTools/vignettes/BayesianTools.html
# https://www.rdocumentation.org/packages/BayesianTools/versions/0.1.3/topics/VSEM

#
library(BayesianTools)

#
cat(file=stderr(), 'Calibrating BASGRA using BayesianTools package', "\n")

# parameter names
bt_names <- parname_BC

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
  # update MaxL
  if (logL1 > logMaxL) {
    logMaxL      <<- logL1
    scparMaxL_BC <<- as.numeric(par)
  }
  # return likelihood
  return(logL1)
}

# construct priors (scaled parameter space)
bt_prior <- createBetaPrior(aa, bb, scparmin_BC[1:np_BC], scparmax_BC[1:np_BC])

# construct setup
bt_setup <- createBayesianSetup(likelihood=bt_likelihood, 
                                prior=bt_prior, 
                                names=bt_names, 
                                parallel=TRUE, 
                                parallelOptions=list(dlls=list(BASGRA_DLL))
                                )

# construct settings (note: DREAMzs has startValue=3 internal chains by default)
nInternal   <- 3 # internal chains for DREAMzs
bt_settings <- list(startValue=nInternal, 
                    iterations=nChain/nChains, 
                    nrChains=nChains, 
                    burnin=nBurnin/nChains*0) # burnin gets discarded but this causes a crash

# run BT
bt_out <- runMCMC(bayesianSetup=bt_setup, 
                  sampler = "DREAMzs", 
                  settings=bt_settings)
cat(file=stderr(), " ", "\n")

# return samples (scaled parameter space)
pChain       <- getSample(bt_out) 
scparMAP_BC  <- MAP(bt_out)$parametersMAP 
cat(file=stderr(), paste("Stored pChain =", dim(pChain)[1], "Iterations with", dim(pChain)[2], "Parameters"), "\n")

# results are generated in BC_BASGRA_BT_results.R
