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
bt_names <- parname_BC
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
                                parallel=T,
                                parallelOptions=list(dlls=list(BASGRA_DLL)),
                                names=bt_names)

# construct settings (note: DREAMzs has startValue=3 internal chains by default)
nInternal   <- 3 # internal chains for DREAMzs
bt_settings <- list(startValue=nInternal, 
                    iterations=nChain/nChains, 
                    nrChains=nChains, 
                    burnin=nBurnin/nChains+nChains, # burnin gets discarded but this causes a crash
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
# browser()

# rerun BT until target mpsrf achieved (provided run time remains reasonable)
target_mpsrf <- 2.0
max_seconds <- 1800
if (FALSE){
  while ((gelmanDiagnostics(bt_out)$mpsrf > target_mpsrf)&
         (bt_out[[1]]$settings$runtime[3] < max_seconds)){
    cat(file=stderr(), paste("Greater than", target_mpsrf, "so continuing..."), "\n")
    stopifnot(nBurnin==0) # restart doesn't currently work with burnin due to bug
    bt_out <- runMCMC(bayesianSetup = bt_out) 
    cat(file=stderr(), " ", "\n")
    # bt_chains <- nInternal * nChains
    bt_length <- dim(bt_out[[1]]$chain[[1]])[[1]]
    # bt_pars <- length(bt_names)
    bt_conv <- gelmanDiagnostics(bt_out)$mpsrf
    # cat(file=stderr(), paste("Total chains =", bt_chains), "\n")
    cat(file=stderr(), paste("Total samples per chain =", bt_length), "\n")
    cat(file=stderr(), paste("Overall convergence (mpsrf) =", round(bt_conv,3)), "\n")
    # browser()
  }
}

# stop parallel
stopParallel(bayesianSetup=bt_setup)

# report convergence
cat(file=stderr(), paste("Convergence of individual parameters (psf)"), "\n")
psf <- gelmanDiagnostics(bt_out)$psrf[,1]
print(round(psf,3))

# memory management
save.image(file="temp.RData")
rm(list=ls())
load(file="temp.RData")

# return samples (scaled parameter space) 
# pChain       <- getSample(bt_out) #### Uses too much memory !
post_num <- nSampling / bt_chains # posterior samples
post_end <- bt_length  
post_start <- bt_length - post_num + 1 
pChain       <- getSample(bt_out, start=post_start, end=post_end) # extract sampling period only
cat(file=stderr(), paste("Stored pChain =", dim(pChain)[1], "Iterations with", dim(pChain)[2], "Parameters"), "\n")

# define ML function
ML <- function(bayesianOutput, ...){
  samples = getSample(bayesianOutput, parametersOnly = F, ...)
  if("mcmcSamplerList" %in% class(bayesianOutput)) nPars <- bayesianOutput[[1]]$setup$numPars
  else nPars = bayesianOutput$setup$numPars
  best = which.max(samples[,nPars + 2])
  return(list(parametersML = samples[best, 1:nPars], valuesML = samples[best, (nPars + 1):(nPars + 3)] ))
}

# store best par
# all <- getSample(bt_out)
scparMAP_BC  <- MAP(bt_out)$parametersMAP 
scparMAP_BC_values  <- MAP(bt_out)$valuesMAP 
logMAP_final <- scparMAP_BC_values[[1]]
scparMaxL_BC <- ML(bt_out)$parametersML 
scparMaxL_BC_values  <- ML(bt_out)$valuesML 
logMaxL_final <- scparMaxL_BC_values[[1]]
params_BC_MAP <- scparMAP_BC * sc
names(params_BC_MAP) <- parname_BC
cat(file=stderr(), paste("MAP parameter values"), "\n")
print(round(params_BC_MAP,4))

# results are generated in BC_BASGRA_BT_results.R

