# Show results of Bayesian calibration of BASGRA model using BayesianTools package

# 
cat(file=stderr(), 'Results of calibrating BASGRA using BayesianTools package', "\n")
suppressMessages({
  library(tidyverse)
  library(coda)
})

# summary and correlation matrix
if (TRUE){
  invisible(capture.output(cmatrix <- summary(bt_out), # this gets correlation matrix into cmatrix
                           file=paste("model_outputs/BC_summary_BT.txt",sep=""))) # and export text summary
  # correlation plot (for most correlated parameters)
  nmatrix <- dim(cmatrix)[1]
  flat <- tibble(row=rep(1:nmatrix, times=nmatrix),
                     col=rep(1:nmatrix,  each=nmatrix),
                     val=as.vector(cmatrix)) %>%
            filter(row>col) %>%
            mutate(absval=abs(val)) %>%
            arrange(desc(absval))
  whichc <- 1:8
  whichp <- unique(c(flat$row[whichc], flat$col[whichc]))
  png( paste("model_outputs/BC_parameters_correlations_BT.png",sep=""),
       width=11*3, height=8*3, units="in", type="windows", res=300)  
  correlationPlot(bt_out, whichParameters=whichp) # parameter correlation plot, very slow and big!
  dev.off()
}

# traceplots (multiple plots)
if (FALSE){
  tracePlot(bt_out) # parameter traces (don't know how to combine onto one sheet)
}

# prior and posterior histograms
if (TRUE){
  png( paste("model_outputs/BC_parameters_histograms_BT.png",sep=""),
       width=11*3, height=8*3, units="in", type="windows", res=300)  
  marginalPlot(bt_out) # prior and posterior histograms (scaled parameters)
  dev.off()
}

# gelman convergence plots (multiple plots)
if (FALSE){
  # gelmanDiagnostics(bt_out, plot=TRUE, 
  #                   start=1, 
  #                   end=nChain/(nInternal*nChains)) # Rhat for each parameter
  gelman.plot(getSample(bt_out, coda=TRUE), ylim=c(1.0,1.5))
}

# prediction function
# https://github.com/florianhartig/BayesianTools/blob/master/Examples/PlotTimeSeriesResults.Rmd
if (FALSE){
  
  bt_samples <- pChain[(nBurnin+1):nChain, ]
  bt_predict <- function(par){
    # use loop from BC_BASGRA_MCMC.R  
    candidatepValues_BC   <- par * sc
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
  s <- 1
  data_col <- 1
  bt_pred_times <- bt_predict(scparMAP_BC)

  # error function
  bt_error <- function(mean, par){
    return(rnorm(length(mean), mean=mean, sd=bt_error_constant)) # copied from VSEM vignette, weird
  }
  
  # plot predictive results for each site
  for (s in 1:nSites){ 
    
    # pdf( paste('model_outputs/BC_calibration_fits_BT_', s, '.pdf',sep=""),
    #      width=pagew, height=pageh)
    png( paste('model_outputs/BC_calibration_fits_BT_', s, '.png',sep=""),
         width=pagew, height=pageh, units="in", type="windows", res=300)

    # set up plot grid
    noutputsMeasured     <- length(unique(data_index[[s]]))
    nrowsPlots           <- ceiling(sqrt(noutputsMeasured+1))
    ncolsPlots           <- ceiling((noutputsMeasured+1)/nrowsPlots)
    oldpar <- par(mfrow=c(nrowsPlots,ncolsPlots),omi=c(0,0,0.5,0), mar=c(2, 2, 2, 1) )
    
    # loop through variables
    for (data_col in unique(data_index[[s]])){ 
      p <- data_col
      datap     <- which( data_name[[s]] == as.character(outputNames[p]) ) # which data points are this variable?
      bt_obs_rows <- list_output_calibr_rows[[s]]
      bt_obs_vals <- rep( as.double(NA), NDAYS )
      bt_obs_vals[bt_obs_rows[datap]] <- data_value[[s]][datap]
      bt_obs_wts <- rep( 0, NDAYS )
      bt_obs_wts[bt_obs_rows[datap]] <- data_weight[[s]][datap]
      bt_obs_vals[bt_obs_wts==0] <- NA # remove unweighted data
      bt_obs_errs <- rep( as.double(NA), NDAYS )
      bt_obs_errs[bt_obs_rows[datap]] <- data_sd[[s]][datap] # note: errors are constant
      bt_error_constant <- data_sd[[s]][datap][1] 
      # bt_obs_times <- data_year[[s]][datap]+(data_doy[[s]][datap]-0.5)/366
      bt_pred_MAP <- bt_predict(scparMAP_BC)
      if (TRUE){
        pred <- getPredictiveIntervals(parMatrix=bt_samples,
                                       model=bt_predict,
                                       numSamples=1000,
                                       quantiles=c(0.025, 0.5, 0.975),
                                       error=bt_error)
        plotTimeSeries(       observed=bt_obs_vals, 
                              predicted = pred$posteriorPredictivePredictionInterval[2,],
                              confidenceBand = pred$posteriorPredictiveCredibleInterval[c(1,3),],
                              predictionBand = pred$posteriorPredictivePredictionInterval[c(1,3),],
                              x=bt_pred_times,
                              main=paste("Site", s, "Variable", data_name[[s]][datap][1])
        )
      }
      if (FALSE){ # this doesn't work with par(mfrow) but gives analysis of residuals
        try(plotTimeSeriesResults(sampler=bt_samples,
                                  model=bt_predict,
                                  observed=bt_obs_vals,
                                  error=bt_error,
                                  main=paste("Site", s, "Variable", data_name[[s]][datap][1])
                                  ),
            silent=TRUE
        )
      }

    } # next data_col
    
    dev.off() # close figure
    par(oldpar)
    
  } # next site
  
}
