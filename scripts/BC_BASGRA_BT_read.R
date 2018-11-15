
# reload BT results
cat(file=stderr(), 'Reading BASGRA calibration from temp.RData', "\n")
load(file="temp.RData")
suppressMessages({
  library(BayesianTools)
  library(coda)
})
dyn.load(BASGRA_DLL) 

# memory management        
library(pryr)
mem_used() # 1.18Gb
mem_objects <- as.data.frame(sort(sapply(ls(), function(x){
  object.size(get(x))
}), decreasing=TRUE)) 

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

# plots are generated in BC_BASGRA_BT_plots.R
