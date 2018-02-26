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
                    burnin=nBurnin/nChains*0, # burnin gets discarded but this causes a crash
                    parallel=TRUE,
                    message=TRUE) 

# run BT
bt_out <- runMCMC(bayesianSetup=bt_setup, 
                  sampler = "DREAMzs", 
                  settings=bt_settings)
cat(file=stderr(), " ", "\n")

# rerun BT
if (FALSE){
  while ((gelmanDiagnostics(bt_out)$mpsrf>1.1)&(bt_out[[1]]$settings$runtime[3]<1800)){
    conv <- gelmanDiagnostics(bt_out)$mpsrf
    cat(file=stderr(), paste("Overall convergence (mpsrf) =", round(conv,3)), "\n")
    cat(file=stderr(), paste("Continuing..."), "\n")
    bt_out <- runMCMC(bayesianSetup=bt_out)
    cat(file=stderr(), " ", "\n")
  }
}

# stop parallel
stopParallel(bayesianSetup=bt_setup)

# report convergence
cat(file=stderr(), paste("Convergence of individual parameters (psf)"), "\n")
psf <- gelmanDiagnostics(bt_out)$psrf[,1]
print(round(psf,3))
conv <- gelmanDiagnostics(bt_out)$mpsrf
cat(file=stderr(), paste("Overall convergence (mpsrf) =", round(conv,3)), "\n")

# memory management
save.image(file="temp.RData")
rm(list=ls())
load(file="temp.RData")

# return samples (scaled parameter space)
pChain       <- getSample(bt_out) 
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

