# restart R and then run this to complete BT run if memory error occurs

load(file="temp.RData")
library(BayesianTools)
library(coda)
dyn.load(BASGRA_DLL) 

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

extraOutputs <- c("TSIZE", "CRT", "LAI", "DM", "RES")
source('scripts/BC_BASGRA_BT_results.R') # new solver results
# which in turn calls:
# source('scripts/plotResiduals_BT.r') # replacement functions

## 3. PRINTING AND PLOTTING ##

cat(file=stderr(), 'Calling BC output scripts', "\n")
source('scripts/BC_export_parModes.R')
# source('scripts/BC_plot_parameters_traceplots.R') # very slow! but useful for convergence check
source('scripts/BC_plot_parameters_priorbeta_histograms.R')
# outputMax[which(outputNames=="CLV")] <- 200 # it's pretty awkward to set these programmatically
# outputMax[which(outputNames=="CST")] <- 100
# outputMax[which(outputNames=="TILTOT")] <- 10000
# outputMax[which(outputNames=="WCL")] <- 60
# source('scripts/BC_plot_outputs_data.R')

# run with MAP parameters (can run shorter time window)
for (s in 1:nSites){
  cat(file=stderr(), paste('Running BASGRA with BC MAP parameters, site',s), "\n")
  source(sitesettings_filenames[[s]])
  # modify sitesettings
  # year_start     <- as.integer(2011)
  # doy_start      <- as.integer(152) # 1 June
  year_stop      <- as.integer(2017)
  doy_stop       <- as.integer(151) # 31 May
  NDAYS_all <- NDAYS # will need to restore this for other functions
  NDAYS          <- as.integer(difftime(
    as.Date(doy_stop -1, origin=paste(year_stop ,1,1,sep='-')), 
    as.Date(doy_start-1, origin=paste(year_start,1,1,sep='-')),
    units="days")) + 1L
  y              <- matrix(0,NDAYS,NOUT)
  # overwrite parameters with MAP
  file_params    <- 'model_outputs/BASGRA_parModes.txt' 
  # temp <- read.csv(file_params, sep="\t")
  parcol         <- 2 + 8*(s-1)
  df_params      <- read.table(file_params,header=T,sep="\t",row.names=1)
  params         <- df_params[,parcol]
  # now run
  output <- run_model()
  file_table  = paste("model_outputs/basgra_trace_table_",s,".txt", sep="")
  file_plot   = paste("model_outputs/basgra_trace_plots_",s,".png", sep="")
  export_output(file_table=file_table, file_plot=file_plot)
  NDAYS <- NDAYS_all
  y     <- matrix(0,NDAYS,NOUT)
}

# save workspace since it takes a long time to generate
cat(file=stderr(), 'Saving BASGRA_Workspace.RData', "\n")
file_save <- 'model_outputs/BASGRA_Workspace.RData' 
save.image(file_save)

# latest output
output <- as.data.frame(output)
names(output) <- outputNames
x <- output$Time
# plot(x,output$DEBUG,main="DEBUG")

# finish
cat(file=stderr(), 'Finished BC_BASGRA.r', "\n")
dyn.unload(BASGRA_DLL) 
# installr::kill_all_Rscript_s() # kills processes left from BT parallel

# reload workspace
if (FALSE){
  file_save <- 'model_outputs/BASGRA_Workspace.RData' 
  load(file_save)  
  BASGRA_DLL <- "model/BASGRA_WG.DLL" # use _WG version to use own PET
  dyn.load(BASGRA_DLL) 
}



