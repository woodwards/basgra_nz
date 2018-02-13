plotTimeSeriesResults <- function(sampler, model, observed, error = NULL, plotResiduals = TRUE, start = 1, prior = FALSE, ...){
  oldPar = par(no.readonly = TRUE)
  
  if (plotResiduals == TRUE && is.null(error)) {
    warning("Can not plot residuals without an error function.")
  }
  
  if (plotResiduals == TRUE && !is.null(error)) {
    layout(matrix(c(1, 1, 1, 2, 3, 4), 2, 3, byrow = TRUE))
    par(mar = c(3, 3, 3, 3), oma = c(2, 2, 2, 2))
  }
  
  # ... can we pass on to both getSample and plot? 
  
  if(prior == FALSE){
    if(inherits(sampler,"bayesianOutput")) parMatrix = getSample(sampler, start = start)
    else if (class(sampler) == "matrix") parMatrix = sampler
    else if ("mcmc.list" %in% class(sampler) || "mcmc" %in% class(sampler)) parMatrix = getSample(sampler, start = start)
    else stop("wrong type given to variable sampler")    
  }else if (prior == TRUE){
    if(inherits(sampler,"bayesianOutput")) {
      if(class(sampler)[1] == "mcmcSamplerList") parMatrix = sampler[[1]]$setup$prior$sampler(1000)
      else parMatrix <- sampler$setup$prior$sampler(1000)
    } else {
      stop("prior==TRUE is only available for sampler of type bayesianOutput") 
    }
  }else stop("BayesianTools::plotTimeSeriesResults - wrong argument to prior")
  
  numSamples = min(1000, nrow(parMatrix))
  
  pred <- getPredictiveIntervals(parMatrix = parMatrix,
                                 model = model,
                                 numSamples = numSamples,
                                 quantiles = c(0.025, 0.5, 0.975),
                                 error = error)
  
  if(!is.null(error)) plotTimeSeries(observed = observed,
                                     predicted = pred$posteriorPredictivePredictionInterval[2,],
                                     confidenceBand = pred$posteriorPredictiveCredibleInterval[c(1,3),],
                                     predictionBand = pred$posteriorPredictivePredictionInterval[c(1,3),],
                                     ...)
  else plotTimeSeries(observed = observed, predicted = pred$posteriorPredictiveSimulations,
                      confidenceBand = pred$posteriorPredictiveSimulations[c(1,3),],
                      ...)
  
  if (plotResiduals && !is.null(error)) {
    dh = getDharmaResiduals(model = model,
                            parMatrix = parMatrix,
                            numSamples = numSamples,
                            observed = observed,
                            error = error,
                            plot = FALSE)
    # qq-plot
    gap::qqunif(dh$scaledResiduals, pch=2, bty="n", logscale = F, col = "black", cex = 0.6, main = "QQ plot residuals", cex.main = 1)
    
    # missing data causes crash
    # browser()
    keeps <- which(!is.na(dh$scaledResiduals))
    if (length(keeps)>1){
      
    # residuals vs fitted
    DHARMa::plotResiduals(dh$fittedPredictedResponse[keeps], dh$scaledResiduals[keeps], main = "Residual vs. predicted\n quantile lines should be\n horizontal lines at 0.25, 0.5, 0.75", cex.main = 1, xlab = "Predicted value", ylab = "Standardized residual")
    
    # residuals vs time
    t <- 1:length(dh$fittedPredictedResponse)
    DHARMa::plotResiduals(t[keeps], dh$scaledResiduals[keeps], xlab = "Time", ylab = "Standardized residual", main = "Residual vs. time\n quantile lines should be\n horizontal lines at 0.25, 0.5, 0.75", cex.main = 1)
    
    message("DHARMa::plotTimeSeriesResults called with posterior predictive (residual) diagnostics. Type vignette(\"DHARMa\", package=\"DHARMa\") for a guide on how to interpret these plots")
    
    } # end check keeps
    
  }
  
  
  par(oldPar)
}

getDharmaResiduals <- function(model, parMatrix, numSamples, observed, error, plot = TRUE){
  
  predDistr <- getPredictiveDistribution(parMatrix = parMatrix,
                                         model = model,
                                         numSamples = numSamples)
  # apply error to predictions
  for (i in 1:nrow(predDistr)){
    predDistr[i,] = error(mean = predDistr[i,], par = parMatrix[i,]) 
  }
  
  fittedPars = apply(parMatrix, 2, median)
  fittedPredictedResponse = model(fittedPars)
  
  dh = DHARMa::createDHARMa(simulatedResponse = t(predDistr),
                            observedResponse = observed,
                            fittedPredictedResponse = fittedPredictedResponse)
  if (plot == TRUE) {
    DHARMa::plotSimulatedResiduals(dh)
  }
  
  return(dh)
}
