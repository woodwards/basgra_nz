# Show results of Bayesian calibration of BASGRA model using BayesianTools package

#
cat(file=stderr(), 'Results of calibrating BASGRA using BayesianTools package', "\n")

# diagnostic results
summary(bt_out) # summary, posteriors, convergence, correlation
# plot(bt_out) # parameter traces
# marginalPlot(bt_out) # prior and posterior histograms
# correlationPlot(bt_out, whichParameters=1:10) # parameter correlation plot 
# gelmanDiagnostics(bt_out, plot=TRUE, start=nBurnin/(nChains*nInternal)) # Rhat for each parameter

# prediction function
# https://github.com/florianhartig/BayesianTools/blob/master/Examples/PlotTimeSeriesResults.Rmd
s          <- (1:nSites)[1] # example
data_col   <- unique(data_index[[s]])[1] # example
bt_predict <- function(par){
  # use loop from BC_BASGRA_MCMC.R  
  candidatepValues_BC   <- par
  # for (s in 1:nSites) {
    params         <- list_params        [[s]] # get site parameters initial values (in parameters.txt)
    matrix_weather <- list_matrix_weather[[s]] # get site weather
    days_harvest   <- list_days_harvest  [[s]] # get site harvest
    NDAYS          <- list_NDAYS         [[s]] # get site NDAYS
    # ip_BC_site[[s]] = indicies of model parameters being changed (in parameters.txt)
    # icol_pChain_site[[s]] = indices of calibration parameters being used (in parameters_BC.txt)
    params[ ip_BC_site[[s]] ] <- candidatepValues_BC[ icol_pChain_site[[s]] ]
    output                    <- run_model(params,matrix_weather,days_harvest,NDAYS)
    # list_output[[s]]          <- output
  # }
  this_output                 <- output[,data_col] 
  return(this_output)
}
bt_samples <- pChain[(nBurnin+1):nChain,]

# plot predictive results
for (s in 1:nsites){ 
  for (data_col in unique(data_index[[s]])){ 
    
    p <- data_col
    datap     <- which( data_name[[s]] == as.character(outputNames[p]) ) # which data points are this variable?
    bt_obs_rows <- list_output_calibr_rows[[s]]
    bt_obs_vals <- rep( as.double(NA), NDAYS )
    bt_obs_vals[bt_obs_rows[datap]] <- data_value[[s]][datap]
    bt_obs_errs <- rep( as.double(NA), NDAYS )
    bt_obs_errs[bt_obs_rows[datap]] <- data_sd[[s]][datap] # note: errors are constant
    bt_error_constant <- data_sd[[s]][datap][1] 
    bt_obs_time <- data_year[[s]][datap]+(data_doy[[s]][datap]-0.5)/366
    # error function
    bt_error <- function(mean, par){
      return(bt_error_constant)
    }
    plotTimeSeriesResults(sampler=bt_samples, 
                          model=bt_predict, 
                          observed=bt_obs_vals,
                          error=bt_error,
                          main=paste("Site", s, "Variable", data_name[[s]][datap][1]),
                          na.rm=TRUE)
    
  }
}