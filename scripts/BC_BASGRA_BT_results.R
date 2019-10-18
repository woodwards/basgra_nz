# Show results of Bayesian calibration of BASGRA model using BayesianTools package

# requires scenario folder
cat(file = stderr(), "Plotting results of calibrating BASGRA using BayesianTools package", "\n")

# reload BT results
cat(file = stderr(), "Reading checkpoint after BASGRA calibration", "\n")
load(file = paste(scenario, "/checkpoint_after_calibration.RData", sep = ""))
suppressMessages({
  library(tidyverse)
  library(ggthemes) # for ggplot
  library(lemon) # helps ggplot
  library(BayesianTools)
  library(BASGRA)
  library(scales)
  library(coda) # helps BayesianTools
})
# dyn.load(BASGRA_DLL)
graphics.off()

# weather summary
cat(file = stderr(), "Weather summary by site\n")
lapply(1:nSites, function(i){
  c(
    list_NDAYS[[i]],
    mean(list_matrix_weather[[i]][1:NDAYS,3]), # rad
    mean(list_matrix_weather[[i]][1:NDAYS,4]), # tmin
    mean(list_matrix_weather[[i]][1:NDAYS,5]), # tmax
    sum(list_matrix_weather[[i]][1:NDAYS,6])/6, # rain
    sum(list_matrix_weather[[i]][1:NDAYS,7])/6 # pet
  )
}) %>% print()

# additional outputs to plot
if (fit_mcmc){
  calibOutputs <- c("BASAL", "CLV", "CST", "CLVD", "CRT", "TILTOT", "WCLM")
  otherOutputs <- c("LAI", "TSIZE", "RGRTV", "RDRTIL", "VERN", "YIELD", "DEBUG")
  otherOutputs <- c("LAI", "TSIZE", "GTILV", "DTILV", "TRANRF", "YIELD", "VERN")
} else {
  otherOutputs <- c("LAI", "TSIZE", "BASAL", "TILTOT", "VERN", "CRT", "YIELD")
}
# otherOutputs <- c("RLEAF", "DEBUG", "RDRL", "VERN", "YIELD", "RES", "HARVFR")

# fix parameters
params
parnames <- rownames(df_params)
parname_BC
fixpari <- match(names(fixpar), parnames)

# memory management
library(pryr)
mem_used() # 1.18Gb
mem_objects <- as.data.frame(sort(sapply(ls(), function(x) {
  object.size(get(x))
}), decreasing = TRUE))

# return samples (scaled parameter space)
# pChain       <- getSample(bt_out) #### Uses too much memory !
post_num <- nSampling / bt_chains # posterior samples
post_end <- bt_length
post_start <- bt_length - post_num + 1
pChain <- getSample(bt_out, start = post_start, end = post_end) # extract sampling period only
cat(file = stderr(), paste("Stored pChain =", dim(pChain)[1], "Iterations with", dim(pChain)[2], "Parameters"), "\n")

# define ML function
ML <- function(bayesianOutput, ...) {
  samples <- getSample(bayesianOutput, parametersOnly = F, ...)
  if ("mcmcSamplerList" %in% class(bayesianOutput)) {
    nPars <- bayesianOutput[[1]]$setup$numPars
  } else {
    nPars <- bayesianOutput$setup$numPars
  }
  best <- which.max(samples[, nPars + 2])
  return(list(parametersML = samples[best, 1:nPars], valuesML = samples[best, (nPars + 1):(nPars + 3)]))
}

# store best par
# all <- getSample(bt_out, start=post_start, end=post_end)
scparMAP_BC <- MAP(bt_out, start = post_start, end = post_end)$parametersMAP
scparMAP_BC_values <- MAP(bt_out, start = post_start, end = post_end)$valuesMAP
logMAP_final <- scparMAP_BC_values[[1]]
scparMaxL_BC <- ML(bt_out, start = post_start, end = post_end)$parametersML
scparMaxL_BC_values <- ML(bt_out, start = post_start, end = post_end)$valuesML
logMaxL_final <- scparMaxL_BC_values[[1]]
params_BC_MAP <- scparMAP_BC * sc
names(params_BC_MAP) <- parname_BC
cat(file = stderr(), paste("MAP parameter values"), "\n")
print(round(params_BC_MAP, 4))

# correlation matrix ####
if (TRUE && fit_mcmc) {
  cat(file = stderr(), "Plot largest cells of posterior correlation grid", "\n")
  cmatrix <- round(cor(pChain), 3)
  nmatrix <- dim(cmatrix)[1]
  flat <- tibble(
    row = rep(1:nmatrix, times = nmatrix),
    col = rep(1:nmatrix, each = nmatrix),
    val = as.vector(cmatrix)
  ) %>%
    filter(row > col) %>%
    mutate(absval = abs(val)) %>%
    arrange(desc(absval))
  top <- 10 # most correlated parameters, retaining order
  whichp <- c()
  whichx <- c()
  # whichx <- match(c("RUBISC", "SIMAX1T"), parname_BC) # additional params to check
  for (i in 1:top) {
    whichp <- c(whichp, flat$row[i], flat$col[i])
  }
  whichp <- unique(whichp, whichx)
  png(paste(scenario, "/BC_parameter_correlations_BT.png", sep = ""),
    width = 11 * 3, height = 8 * 3, units = "in", type = "windows", res = 300
  )
  correlationPlot(bt_out, start = post_start, end = post_end, whichParameters = whichp) # parameter correlation plot, very slow and big!
  dev.off()
}

# traceplots (multiple plots) ####
if (TRUE && fit_mcmc) {
  # modified coda::traceplot code : edit title, include sc
  traceplot <-
    function(x, smooth = FALSE, col = 1:6, type = "l", xlab = "Iterations",
                 ylab = "", ...) {
      x <- mcmc.list(x)
      args <- list(...)
      for (j in 1:nvar(x)) {
        xp <- as.vector(time(x))
        yp <- if (nvar(x) > 1) {
          x[, j, drop = TRUE]
        } else {
          x
        }
        yp <- do.call("cbind", yp)
        matplot(xp, yp * sc[j], xlab = xlab, ylab = ylab, type = type, col = col, ...)
        if (!is.null(varnames(x)) && is.null(list(...)$main)) {
          title(paste(varnames(x)[j]))
        }
        if (smooth) {
          scol <- rep(col, length = nchain(x))
          for (k in 1:nchain(x)) {
            lines(lowess(xp, yp[, k]),
              col = scol[k]
            )
          }
        }
      }
    }
  cat(file = stderr(), "Plot parameter traces", "\n")
  pagew <- 11
  pageh <- 8
  png(paste(scenario, "/BC_parameter_traces.png", sep = ""),
    width = pagew, height = pageh, units = "in", type = "windows", res = 300
  )
  nrowsPlots <- ceiling(sqrt((np_BC + 1) * pageh / pagew))
  ncolsPlots <- ceiling((np_BC + 1) / nrowsPlots)
  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(nrowsPlots, ncolsPlots))
  par(mar = c(2, 2, 2, 1))
  traceplot(getSample(bt_out, start = post_start, end = post_end, coda = TRUE))
  dev.off()
  par(oldpar)
}

# gelman convergence plots (multiple plots) ####
# only works on sampling period, not burnin, so not very useful
# if (FALSE){
#   # gelmanDiagnostics(bt_out, start=post_start, end=post_end, plot=TRUE)
#   cat(file=stderr(), "Plot Gelman-Rubin statistic", "\n")
#   pagew <- 11 ; pageh <- 8
#   png( paste(scenario, "/BC_parameter_gelman.png", sep=""),
#        width=pagew, height=pageh, units="in", type="windows", res=300)
#   nrowsPlots <- ceiling( sqrt((np_BC+1)*pageh/pagew) )
#   ncolsPlots <- ceiling( (np_BC+1)/nrowsPlots )
#   oldpar <- par(no.readonly = TRUE)
#   par( mfrow = c(nrowsPlots,ncolsPlots) )
#   par(mar=c(2, 2, 2, 1))
#   gelman.plot(getSample(bt_out, start=post_start, end=post_end, coda=TRUE),
#               ylim=c(1.0,2.0),
#               ask=FALSE,
#               autoburnin=FALSE, # TRUE=ignores first half of chain
#               auto.layout=FALSE, # FALSE=one plot at a time
#               bin.width=50)
#   plot(1,type='n', axes=FALSE, xlab="", ylab="")
#   plot_colors <- c("black","red")
#   legend("bottomright", c("Gelman-R","95% C.I."),
#          bty="n", lty=c(1,2), lwd=c(1,2), col=plot_colors, title = "LEGEND:")
#   dev.off()
#   par(oldpar)
#
# }

# prediction function (old method) ####
# https://github.com/florianhartig/BayesianTools/blob/master/Examples/PlotTimeSeriesResults.Rmd
# if (FALSE){
#
#   source("scripts/plotResiduals_BT.r") # replacement functions
#
#   bt_predict <- function(par){ # needs s and data_col
#     # use loop from BC_BASGRA_MCMC.R
#     candidatepValues_BC   <- par * sc
#     # for (s in 1:nSites) {
#       params         <- list_params        [[s]] # get site parameters initial values (in parameters.txt)
#       matrix_weather <- list_matrix_weather[[s]] # get site weather
#       days_harvest   <- list_days_harvest  [[s]] # get site harvest
#       NDAYS          <- list_NDAYS         [[s]] # get site NDAYS
#       # ip_BC_site[[s]] = indicies of model parameters being changed (in parameters.txt)
#       # icol_pChain_site[[s]] = indices of calibration parameters being used (in parameters_BC.txt)
#       params[ ip_BC_site[[s]] ] <- candidatepValues_BC[ icol_pChain_site[[s]] ]
#       output                    <- run_model(params,matrix_weather,days_harvest,NDAYS,NOUT,matrix(0,NDAYS,NOUT))
#       # list_output[[s]]          <- output
#     # }
#     this_output                 <- output[,data_col]
#     this_output[is.na(this_output)] <- -999 # catch NA
#     return(this_output)
#   }
#
#   # error function
#   bt_error <- function(mean, par){
#     return(rnorm(length(mean), mean=mean, sd=bt_error_constant)) # copied from VSEM vignette, weird
#   }
#
#   # statistics calculations
#   calc_rmse <- function(m,d){
#     if (length(m)==0 && length(d)==0){
#       NA_real_
#     } else {
#       sqrt(mean((m-d)^2, na.rm=TRUE))
#     }
#   }
#   calc_rsq <- function(m,d){
#     if (length(m)==0 && length(d)==0){
#       NA_real_
#     } else {
#       d[is.na(m)] <- NA
#       1-mean((m-d)^2, na.rm=TRUE)/var(d, na.rm=TRUE)
#     }
#   }
#   # https://stackoverflow.com/questions/17549762/is-there-such-colsd-in-r
#   colSds <- function(x, na.rm=TRUE) {
#     if (is.null(dim(x))){ # vector
#       return(sd(x, na.rm=na.rm))
#     } else if (na.rm) {
#       n <- colSums(!is.na(x)) # thanks @flodel
#     } else {
#       n <- nrow(x)
#     }
#     colVar <- colMeans(x*x, na.rm=na.rm) - (colMeans(x, na.rm=na.rm))^2
#     return(sqrt(colVar * n/(n-1)))
#   }
#
#   # plot predictive results for each site (and collect residuals into a dataframe) ####
#   residual_df <- vector("list", nSites * nBCvar)
#   sitenames <- gsub( "\\.R", "", sub(".*BASGRA_","",sitesettings_filenames), ignore.case=TRUE )
#   site_name <- case_when(
#     str_detect(sitenames, "Northland") ~ "Northland",
#     str_detect(sitenames, "Scott") ~ "Waikato",
#     str_detect(sitenames, "Lincoln") ~ "Canterbury")
#   sitenames <- factor(sitenames, sitenames) # set in BC_BASGRA_MCMC_init.R, preserve order
#   s <- 1
#   for (s in 1:nSites){
#
#     # predictins against data
#     cat(file=stderr(), "Plot model calibration fits against data, site", s, "\n")
#     # pdf( paste("model_outputs/BC_calibration_fits_BT_", s, ".pdf",sep=""),
#     #      width=pagew, height=pageh)
#     png( paste(scenario, "/BC_calibration_fits_BT_", s, ".png", sep=""),
#          width=11, height=8, units="in", type="windows", res=300)
#
#     # set up plot grid
#     noutputsMeasured     <- length(unique(data_index[[s]]))
#     nrowsPlots           <- ceiling(sqrt(noutputsMeasured+1))
#     ncolsPlots           <- ceiling((noutputsMeasured+1)/nrowsPlots)
#     oldpar <- par(mfrow=c(nrowsPlots,ncolsPlots),omi=c(0,0,0.5,0), mar=c(2, 2, 2, 1) )
#
#     # loop through calibration variables
#     data_col <- 1
#     bt_pred_times <- bt_predict(scparMAP_BC)
#     i <- 1
#     for (i in 1:nBCvar){
#       data_col <- unique(data_index[[s]])[[i]] # corresponding column in output
#       datap     <- which( data_name[[s]] == as.character(outputNames[data_col]) ) # which data points are this variable?
#       bt_obs_rows <- list_output_calibr_rows[[s]][datap] # corresponding rows in output
#       bt_obs_vals <- data_value[[s]][datap]
#       bt_obs_wts <- data_weight[[s]][datap]
#       bt_obs_errs <- data_sd[[s]][datap]
#       bt_error_constant <- data_sd[[s]][datap][1] # used in bt_error
#       bt_obs_times <- data_year[[s]][datap]+(data_doy[[s]][datap]-0.5)/366
#       bt_pred_MAP <- bt_predict(scparMAP_BC)
#       bt_pred_MAP_obs <- bt_pred_MAP[bt_obs_rows]
#       keeps <- bt_obs_wts>0
#       rmse <- calc_rmse(bt_pred_MAP_obs[keeps], bt_obs_vals[keeps])
#       rsq <- calc_rsq(bt_pred_MAP_obs[keeps], bt_obs_vals[keeps])
#       keeps2 <- bt_obs_wts==0 | bt_obs_wts>0
#       rmse2 <- calc_rmse(bt_pred_MAP_obs[keeps2], bt_obs_vals[keeps2])
#       rsq2 <- calc_rsq(bt_pred_MAP_obs[keeps2], bt_obs_vals[keeps2])
#       bt_pred_ML <- bt_predict(scparMaxL_BC)
#       scparMode_BC <- parmode_BC / sc
#       bt_pred_Mode <- bt_predict(scparMode_BC)
#       if (TRUE){
#         pred <- getPredictiveIntervals(parMatrix=pChain,
#                                        model=bt_predict,
#                                        numSamples=1000,
#                                        quantiles=c(0.05, 0.5, 0.95),
#                                        error=bt_error)
#         # error function gets sampled on top of parameter variation
#         confidenceBand <- pred$posteriorPredictiveCredibleInterval[c(1,3),]
#         confidenceBand <- pmax(confidenceBand, 0.0)
#         predictedMedian <- pred$posteriorPredictiveCredibleInterval[2,]
#         predictedMedian_obs <- predictedMedian[bt_obs_rows]
#         predicted <- pred$posteriorPredictivePredictionInterval[2,]
#         predictionBand <- pred$posteriorPredictivePredictionInterval[c(1,3),]
#         predictionBand <- pmax(predictionBand, 0.0)
#         x <- pred$posteriorPredictiveSimulations[,bt_obs_rows]
#         if (length(bt_obs_rows)>1){
#           bt_pred_mean <- colMeans(x)
#           bt_pred_sd <- colSds(x)
#         } else if (length(bt_obs_rows)==1) {
#           bt_pred_mean <- mean(x)
#           bt_pred_sd <- sd(x)
#         } else {
#           bt_pred_mean <- NA
#           bt_pred_sd <- NA
#         }
#         plotTimeSeries <- function(observed = NULL, predicted = NULL, x = NULL, xlim = NULL,
#                                    confidenceBand = NULL, predictionBand = NULL,
#                                    xlab = "Time", ylab = "Observed / predicted values", ...){
#           ylim = range(observed, predicted, confidenceBand, predictionBand,na.rm=TRUE)
#           # ylim = range(observed, predicted, na.rm=TRUE)
#           if (is.null(x)){
#             if(!is.null(observed)) x = 1:length(observed)
#             else if(!is.null(predicted)) x = 1:length(predicted)
#             else stop("either observed or predicted must be supplied")
#           }
#           len = length(x)
#           plot(x, y=rep(0,len), xlim = xlim, ylim = ylim, type = "n", xlab = xlab, ylab = ylab, ...)
#           if(!is.null(predictionBand))
#             polygon(c(x,rev(x)),c(predictionBand[1,],predictionBand[2,len:1]),col="moccasin",border=NA)
#           # polygon(c(1:len,len:1),c(predictionBand[1,],predictionBand[2,len:1]),col="moccasin",border=NA)
#           if(!is.null(confidenceBand))
#             polygon(c(x,rev(x)),c(confidenceBand[1,],confidenceBand[2,len:1]),col="#99333380",border=NA)
#           # polygon(c(1:len,len:1),c(confidenceBand[1,],confidenceBand[2,len:1]),col="#99333380",border=NA)
#           if(!is.null(predicted)) lines(x, predicted, col = "red")
#           if(!is.null(observed)) points(x, observed, col = NA, pch = 3, cex = 0.6)
#         }
#         plotTimeSeries(       observed = rep_len(bt_obs_vals, length(predicted)), # only used to set the range
#                               predicted = predicted,
#                               confidenceBand = confidenceBand,
#                               predictionBand = predictionBand,
#                               x=bt_pred_times,
#                               # xlim=c(2012,2015), # show only a subset of time line (else = NULL)
#                               # main=paste(easyNames[data_col], outputUnits[data_col],
#                               #            "RSME_MAP =", signif(rmse,3), "/", signif(rmse2,3),
#                               #            "RSQ_MAP =", signif(rsq,3), "/", signif(rsq2,3))
#                               main=paste(easyNames[data_col], outputUnits[data_col])
#         )
#         # plot key prediction lines
#         lines(x=bt_pred_times, y=bt_pred_Mode, col=NA)
#         lines(x=bt_pred_times, y=bt_pred_ML, col="lightblue")
#         lines(x=bt_pred_times, y=bt_pred_MAP, col="blue")
#         abline(h=0, col="darkgrey")
#         # plot all data
#         keeps <- bt_obs_wts==0 | bt_obs_wts>0
#         x_obs <- bt_obs_times[keeps]
#         suppressWarnings({
#           # arrows(x0=x_obs, y0=bt_obs_vals[keeps],
#           #        x1=x_obs, y1=bt_pred_MAP_obs[keeps],
#           #        col="black", lwd=1.5, angle=45, length=0.05) # residual arrow
#           # arrows(x0=x_obs, y0=bt_obs_vals[keeps]-bt_obs_errs[keeps]*1.96,
#           #        x1=x_obs, y1=bt_obs_vals[keeps]+bt_obs_errs[keeps]*1.96,
#           #        col="grey", lwd=1.5, angle=90, code=3, length=0.05) # error bars
#         })
#         points( x=x_obs, y=bt_obs_vals[keeps],
#                 pch=16, col="grey", cex=1.5)
#         # plot weighted data
#         keeps <- bt_obs_wts>0
#         x_obs <- bt_obs_times[keeps]
#         suppressWarnings({
#           arrows(x0=x_obs, y0=bt_obs_vals[keeps],
#                  x1=x_obs, y1=predictedMedian_obs[keeps],
#                  col="black", lwd=1.5, angle=45, length=0.05) # residual arrow
#           # arrows(x0=x_obs, y0=bt_obs_vals[keeps]-bt_obs_errs[keeps]*1.96,
#           #        x1=x_obs, y1=bt_obs_vals[keeps]+bt_obs_errs[keeps]*1.96,
#           #        col="darkblue", lwd=1.5, angle=90, code=3, length=0.05) # error bars
#         })
#         points( x=x_obs, y=bt_obs_vals[keeps],
#                 pch=16, col="darkblue", cex=1.5)
#         # collect residual_df
#         temp <- tibble(
#           site_num=s,
#           site_name=sitenames[s],
#           var_name=easyNames[data_col],
#           var_units=outputUnits[data_col],
#           times=bt_pred_times[bt_obs_rows],
#           obs_vals=data_value[[s]][datap],
#           obs_errs=data_sd[[s]][datap],
#           obs_wts=data_weight[[s]][datap],
#           pred_map=bt_pred_MAP[bt_obs_rows],
#           pred_med=predictedMedian[bt_obs_rows],
#           pred_min=confidenceBand[1,bt_obs_rows],
#           pred_max=confidenceBand[2,bt_obs_rows],
#           pred_min2=predictionBand[1,bt_obs_rows],
#           pred_max2=predictionBand[2,bt_obs_rows],
#           pred_mean=bt_pred_mean,
#           pred_sd=bt_pred_sd,
#           resid_mean=pred_mean-obs_vals,
#           resid_sd=pred_sd,
#           resid_rmse=sqrt(pred_sd^2+2*resid_mean*pred_mean+(obs_vals^2-pred_mean^2)),
#           rmse=rmse,
#           rmse2=rmse2
#         )
#         residual_df[[(s-1)*nBCvar+i]] <- temp
#       }
#
#
#     } # next data_col
#
#     # legend and title
#     plot(1, type="n", axes=FALSE, xlab="", ylab="") # empty plot with legend
#     legend( "bottomright", title="Predictions",
#             legend=c("Prior Mode", "Median", "Max L",      "MAP",      "Calib Data", "Other Data", "Residuals"),
#             # col   =c(NA,  "red",    "lightblue",  "blue",     "darkblue",   "grey",       "black"),
#             col   =c(NA,  "red",    "lightblue",  "blue",     "darkblue",   "grey",       "black"),
#             lty=1, lwd=1)
#     mtext( paste("SITE ",s," (",site_name[s],")",sep=""),
#            side=3, line=1, outer=TRUE, cex=1, font=2)
#
#     # close figure
#     dev.off()
#     par(oldpar)
#
#     # plot residual analysis (slow) ####
#     # FIXME doesn't work if more than one obs one a row
#     # if (FALSE){
#     #   cat(file=stderr(), "Plot residual analysis, site", s, "\n")
#     #   i <- 1
#     #   for (i in 1:nBCvar){
#     #     data_col <- unique(data_index[[s]])[[i]] # corresponding column in output
#     #     datap     <- which( data_name[[s]] == as.character(outputNames[data_col]) ) # which data points are this variable?
#     #     bt_obs_rows <- list_output_calibr_rows[[s]][datap] # corresponding rows in output
#     #     bt_obs_vals <- data_value[[s]][datap]
#     #     bt_obs_wts <- data_weight[[s]][datap]
#     #     bt_obs_errs <- data_sd[[s]][datap]
#     #     bt_error_constant <- data_sd[[s]][datap][1] # used in bt_error
#     #     bt_obs_times <- data_year[[s]][datap]+(data_doy[[s]][datap]-0.5)/366
#     #     if (TRUE){
#     #       # this doesn't work with par(mfrow) but gives analysis of residuals
#     #       # try({ # ignore errors thrown some subplots
#     #       # debug(plotTimeSeriesResults)
#     #       suppressMessages({
#     #         plotTimeSeriesResults(sampler=pChain,
#     #                               model=bt_predict,
#     #                               observed=bt_obs_vals,
#     #                               error=bt_error,
#     #                               main=paste("Site", s, easyNames[data_col]," ",outputUnits[data_col])
#     #                               )
#     #       })
#     #       # }, silent=TRUE)
#     #       # save
#     #       dev.copy(png, paste("model_outputs/Residuals_Site_", s, "_", easyNames[data_col], ".png", sep=""),
#     #                width = 480*2, height = 480*2)
#     #       dev.off()
#     #     }
#     #   } # next data_col
#     # }
#
#     # plot other model outputs for each site ####
#     if (TRUE){
#       cat(file=stderr(), "Plot model calibration other outputs, site", s, "\n")
#       png( paste(scenario, "/BC_calibration_other_BT_", s, ".png", sep=""),
#            width=11, height=8, units="in", type="windows", res=300)
#
#       # set up plot grid
#       data_cols <- match(otherOutputs, outputNames)
#       noutputsMeasured     <- length(data_cols)
#       nrowsPlots           <- ceiling(sqrt(noutputsMeasured+1))
#       ncolsPlots           <- ceiling((noutputsMeasured+1)/nrowsPlots)
#       oldpar <- par(mfrow=c(nrowsPlots,ncolsPlots),omi=c(0,0,0.5,0), mar=c(2, 2, 2, 1) )
#
#       # loop through other selected variables (there are no observations for these variables)
#       for (data_col in data_cols){
#         bt_error_constant <- 0 # used by bt_error
#         bt_pred_MAP <- bt_predict(scparMAP_BC)
#         bt_pred_ML <- bt_predict(scparMaxL_BC)
#         scparMode_BC <- parmode_BC / sc
#         bt_pred_Mode <- bt_predict(scparMode_BC)
#         if (TRUE){
#           pred <- getPredictiveIntervals(parMatrix=pChain,
#                                          model=bt_predict,
#                                          numSamples=1000,
#                                          quantiles=c(0.05, 0.5, 0.95),
#                                          error=bt_error)
#           confidenceBand <- pred$posteriorPredictiveCredibleInterval[c(1,3),]
#           confidenceBand <- pmax(confidenceBand, 0.0)
#           predictedMedian <- pred$posteriorPredictiveCredibleInterval[2,]
#           predicted <- pred$posteriorPredictivePredictionInterval[2,]
#           predictionBand <- pred$posteriorPredictivePredictionInterval[c(1,3),]
#           predictionBand <- pmax(predictionBand, 0.0)
#           plotTimeSeries(       predicted = predicted,
#                                 confidenceBand = confidenceBand,
#                                 predictionBand = predictionBand,
#                                 x=bt_pred_times,
#                                 # xlim=c(2012,2015), # show only a subset of time line (else = NULL)
#                                 main=paste(easyNames[data_col], outputUnits[data_col])
#           )
#           # plot key prediction lines
#           lines(x=bt_pred_times, y=bt_pred_Mode, col=NA)
#           lines(x=bt_pred_times, y=bt_pred_ML, col="lightblue")
#           lines(x=bt_pred_times, y=bt_pred_MAP, col="blue")
#           abline(h=0, col="darkgrey")
#         }
#
#       } # next data_col
#
#       # legend and title
#       plot(1, type="n", axes=FALSE, xlab="", ylab="") # empty plot with legend
#       legend( "bottomright", title="Predictions",
#               legend=c("Prior Mode", "Median", "Max L",      "MAP",      "Calib Data", "Other Data", "Residuals"),
#               # col   =c(NA,  "red",    "lightblue",  "blue",     "darkblue",   "grey",       "black"),
#               col   =c(NA,  "red",    "lightblue",  "blue",     "darkblue",   "grey",       "black"),
#               lty=1, lwd=1)
#       mtext( paste("SITE ",s," (",site_name[s],")",sep=""),
#              side=3, line=1, outer=TRUE, cex=1, font=2)
#
#       # close figure
#       dev.off()
#       par(oldpar)
#     }
#
#   } # next site
#
#   # tidy up residual_df
#   residual_df <- bind_rows(residual_df) %>%
#     mutate(
#       site_name=case_when(
#         str_detect(site_name, "Northland") ~ "Northland",
#         str_detect(site_name, "Scott") ~ "Waikato",
#         str_detect(site_name, "Lincoln") ~ "Canterbury"),
#       site_name=factor(site_name, levels=c("Northland", "Waikato", "Canterbury")),
#       obs_min=obs_vals-obs_errs*1.96,
#       obs_max=obs_vals+obs_errs*1.96
#     ) %>%
#     group_by(var_name) %>%
#     mutate(
#       all_min=min(obs_min, pred_min, pred_min2),
#       all_max=max(obs_max, pred_max, pred_max2),
#       score=case_when(
#         obs_vals==pred_med ~ 0,
#         obs_vals==pred_max ~ 1,
#         obs_vals==pred_min ~ -1,
#         obs_vals==pred_max2 ~ 2,
#         obs_vals==pred_min2 ~ -2,
#         obs_vals>pred_med & obs_vals<pred_max ~ (obs_vals-pred_med)/(pred_max-pred_med),
#         obs_vals<pred_med & obs_vals>pred_min ~ (obs_vals-pred_med)/(pred_med-pred_min),
#         obs_vals>pred_max & obs_vals<pred_max2 ~ (obs_vals-pred_max)/(pred_max2-pred_max)+1,
#         obs_vals<pred_min & obs_vals>pred_min2 ~ (obs_vals-pred_min)/(pred_min-pred_min2)-1,
#         obs_vals>pred_max2 ~ (obs_vals-pred_max2)/(obs_errs*1.96)+2,
#         obs_vals<pred_min2 ~ (obs_vals-pred_min2)/(obs_errs*1.96)-2)
#     ) %>%
#     ungroup()
#
# }

# prediction function (for ggplot) ####
if (TRUE) {
  
  bt_predict <- function(par) { # needs s and data_col
    # use loop from BC_BASGRA_MCMC.R
    candidatepValues_BC <- par * sc
    # for (s in 1:nSites) {
    params <- list_params        [[s]] # get site parameters initial values (in parameters.txt)
    matrix_weather <- list_matrix_weather[[s]] # get site weather
    days_harvest <- list_days_harvest  [[s]] # get site harvest
    NDAYS <- list_NDAYS         [[s]] # get site NDAYS
    # ip_BC_site[[s]] = indicies of model parameters being changed (in parameters.txt)
    # icol_pChain_site[[s]] = indices of calibration parameters being used (in parameters_BC.txt)
    params[ ip_BC_site[[s]] ] <- candidatepValues_BC[ icol_pChain_site[[s]] ]
    if (fixpars){
      params[fixpari] <- fixpar
    }
    output <- run_model(params, matrix_weather, days_harvest, NDAYS, NOUT, matrix(0, NDAYS, NOUT))
    # list_output[[s]]          <- output
    # }
    this_output <- output[, data_col]
    this_output[is.na(this_output)] <- -999 # catch NA
    return(this_output)
  }

  # error function
  bt_error <- function(mean, par) {
    return(rnorm(length(mean), mean = mean, sd = bt_error_constant)) # copied from VSEM vignette, weird
  }

  # statistics calculations
  calc_rmse <- function(m, d) {
    if (length(m) == 0 && length(d) == 0) {
      NA_real_
    } else {
      sqrt(mean((m - d)^2, na.rm = TRUE))
    }
  }
  calc_rsq <- function(m, d) {
    if (length(m) == 0 && length(d) == 0) {
      NA_real_
    } else {
      d[is.na(m)] <- NA
      1 - mean((m - d)^2, na.rm = TRUE) / var(d, na.rm = TRUE)
    }
  }
  # https://stackoverflow.com/questions/17549762/is-there-such-colsd-in-r
  colSds <- function(x, na.rm = TRUE) {
    if (is.null(dim(x))) { # vector
      return(sd(x, na.rm = na.rm))
    } else if (na.rm) {
      n <- colSums(!is.na(x)) # thanks @flodel
    } else {
      n <- nrow(x)
    }
    colVar <- colMeans(x * x, na.rm = na.rm) - (colMeans(x, na.rm = na.rm))^2
    return(sqrt(colVar * n / (n - 1)))
  }

  # predictive results for each site (collect results and residuals into dataframes) ####
  nOvar <- length(otherOutputs)
  result_df <- vector("list", nSites * nBCvar)
  other_df <- vector("list", nSites * nOvar)
  residual_df <- vector("list", nSites * nBCvar)
  sitenames <- factor(sitenames, sitenames) # set in BC_BASGRA_MCMC_init.R, preserve order
  s <- 1
  for (s in 1:nSites) {

    # predictions against data
    cat(file = stderr(), "Get model calibration fits against data, site", s, "\n")

    # loop through calibration variables
    data_col <- 1
    bt_pred_times <- bt_predict(scparMAP_BC)
    i <- 1
    for (i in 1:nBCvar) {
      data_col <- unique(data_index[[s]])[[i]] # corresponding column in output
      datap <- which(data_name[[s]] == as.character(outputNames[data_col])) # which data points are this variable?
      bt_obs_rows <- list_output_calibr_rows[[s]][datap] # corresponding rows in output
      bt_obs_vals <- data_value[[s]][datap]
      bt_obs_wts <- data_weight[[s]][datap]
      bt_obs_errs <- data_sd[[s]][datap]
      bt_error_constant <- data_sd[[s]][datap][1] # used in bt_error
      bt_obs_times <- data_year[[s]][datap] + (data_doy[[s]][datap] - 0.5) / 366
      bt_pred_MAP <- bt_predict(scparMAP_BC)
      bt_pred_MAP_obs <- bt_pred_MAP[bt_obs_rows]
      keeps <- bt_obs_wts > 0
      rmse <- calc_rmse(bt_pred_MAP_obs[keeps], bt_obs_vals[keeps])
      rsq <- calc_rsq(bt_pred_MAP_obs[keeps], bt_obs_vals[keeps])
      keeps2 <- bt_obs_wts == 0 | bt_obs_wts > 0
      rmse2 <- calc_rmse(bt_pred_MAP_obs[keeps2], bt_obs_vals[keeps2])
      rsq2 <- calc_rsq(bt_pred_MAP_obs[keeps2], bt_obs_vals[keeps2])
      bt_pred_ML <- bt_predict(scparMaxL_BC)
      scparMode_BC <- parmode_BC / sc
      bt_pred_Mode <- bt_predict(scparMode_BC)
      if (TRUE) {
        pred <- getPredictiveIntervals(
          parMatrix = pChain,
          model = bt_predict,
          numSamples = 1000,
          quantiles = c(0.05, 0.5, 0.95),
          error = bt_error
        )
        # error function gets sampled on top of parameter variation
        confidenceBand <- pred$posteriorPredictiveCredibleInterval[c(1, 3), ]
        confidenceBand <- pmax(confidenceBand, 0.0)
        predictedMedian <- pred$posteriorPredictiveCredibleInterval[2, ]
        predictedMedian_obs <- predictedMedian[bt_obs_rows]
        predicted <- pred$posteriorPredictivePredictionInterval[2, ]
        predictionBand <- pred$posteriorPredictivePredictionInterval[c(1, 3), ]
        predictionBand <- pmax(predictionBand, 0.0)
        x <- pred$posteriorPredictiveSimulations[, bt_obs_rows]
        if (length(bt_obs_rows) > 1) {
          bt_pred_mean <- colMeans(x)
          bt_pred_sd <- colSds(x)
        } else if (length(bt_obs_rows) == 1) {
          bt_pred_mean <- mean(x)
          bt_pred_sd <- sd(x)
        } else {
          bt_pred_mean <- NA
          bt_pred_sd <- NA
        }
        # collect results_df
        temp <- tibble(
          site_num = s,
          site_name = sitenames[s],
          var_name0 = outputNames[data_col],
          var_name = easyNames[data_col],
          var_units = outputUnits[data_col],
          var_name2 = paste(var_name, var_units),
          times = bt_pred_times,
          pred_map = bt_pred_MAP,
          pred_med = predictedMedian,
          pred_med_mean = mean(predictedMedian), 
          pred_min = confidenceBand[1, ],
          pred_max = confidenceBand[2, ],
          pred_min2 = predictionBand[1, ],
          pred_max2 = predictionBand[2, ],
          obs_vals = data_value[[s]][datap][match(bt_pred_times, bt_obs_times)],
          obs_errs = data_sd[[s]][datap][match(bt_pred_times, bt_obs_times)],
          obs_wts = data_weight[[s]][datap][match(bt_pred_times, bt_obs_times)],
          pred_mean = bt_pred_mean[match(bt_pred_times, bt_obs_times)],
          pred_sd = bt_pred_sd[match(bt_pred_times, bt_obs_times)],
          resid_mean = pred_mean - obs_vals,
          resid_sd = pred_sd,
          resid_rmse = sqrt(pred_sd^2 + 2 * resid_mean * pred_mean + (obs_vals^2 - pred_mean^2)),
          rmse = rmse,
          rmse2 = rmse2
        )
        result_df[[(s - 1) * nBCvar + i]] <- temp
        # collect residual_df
        temp <- tibble(
          site_num = s,
          site_name = sitenames[s],
          var_name0 = outputNames[data_col],
          var_name = easyNames[data_col],
          var_units = outputUnits[data_col],
          var_name2 = paste(var_name, var_units),
          times = bt_pred_times[bt_obs_rows],
          obs_vals = data_value[[s]][datap],
          obs_errs = data_sd[[s]][datap],
          obs_wts = data_weight[[s]][datap],
          pred_map = bt_pred_MAP[bt_obs_rows],
          pred_med = predictedMedian[bt_obs_rows],
          pred_min = confidenceBand[1, bt_obs_rows],
          pred_max = confidenceBand[2, bt_obs_rows],
          pred_min2 = predictionBand[1, bt_obs_rows],
          pred_max2 = predictionBand[2, bt_obs_rows],
          pred_mean = bt_pred_mean,
          pred_sd = bt_pred_sd,
          resid_mean = pred_mean - obs_vals,
          resid_sd = pred_sd,
          resid_rmse = sqrt(pred_sd^2 + 2 * resid_mean * pred_mean + (obs_vals^2 - pred_mean^2)),
          rmse = rmse,
          rmse2 = rmse2
        )
        residual_df[[(s - 1) * nBCvar + i]] <- temp
      } # if
    } # next data_col

    # predictions against other outputs
    cat(file = stderr(), "Get model calibration other outputs, site", s, "\n")
    data_cols <- match(otherOutputs, outputNames)

    # loop through other selected variables (there are no observations for these variables)
    for (i in seq_along(data_cols)) {
      data_col <- data_cols[i]
      bt_error_constant <- 0 # used by bt_error
      bt_pred_MAP <- bt_predict(scparMAP_BC)
      bt_pred_ML <- bt_predict(scparMaxL_BC)
      scparMode_BC <- parmode_BC / sc
      bt_pred_Mode <- bt_predict(scparMode_BC)
      if (TRUE) {
        pred <- getPredictiveIntervals(
          parMatrix = pChain,
          model = bt_predict,
          numSamples = 1000,
          quantiles = c(0.05, 0.5, 0.95),
          error = bt_error
        )
        confidenceBand <- pred$posteriorPredictiveCredibleInterval[c(1, 3), ]
        confidenceBand <- pmax(confidenceBand, 0.0)
        predictedMedian <- pred$posteriorPredictiveCredibleInterval[2, ]
        predicted <- pred$posteriorPredictivePredictionInterval[2, ]
        predictionBand <- pred$posteriorPredictivePredictionInterval[c(1, 3), ]
        predictionBand <- pmax(predictionBand, 0.0)

        # collect other_df
        temp <- tibble(
          site_num = s,
          site_name = sitenames[s],
          var_name0 = outputNames[data_col],
          var_name = easyNames[data_col],
          var_units = outputUnits[data_col],
          var_name2 = paste(var_name, var_units),
          times = bt_pred_times,
          pred_map = bt_pred_MAP,
          pred_med = predictedMedian,
          pred_med_mean = mean(predictedMedian),
          pred_min = confidenceBand[1, ],
          pred_max = confidenceBand[2, ],
          pred_min2 = predictionBand[1, ],
          pred_max2 = predictionBand[2, ],
        )
        other_df[[(s - 1) * nOvar + i]] <- temp
      } # if
    } # next data_col
  } # next site

  # tidy up results_df
  result_df <- bind_rows(result_df) %>%
    mutate(
      obs_min = obs_vals - obs_errs * 1.96,
      obs_max = obs_vals + obs_errs * 1.96,
      logl = flogLi(pred_map, obs_vals, obs_errs, obs_wts) - flogLi(obs_vals, obs_vals, obs_errs, obs_wts)
    ) %>%
    group_by(var_name) %>%
    mutate(
      all_min = min(obs_min, pred_min, pred_min2),
      all_max = max(obs_max, pred_max, pred_max2),
      score = case_when(
        obs_vals == pred_med ~ 0,
        obs_vals == pred_max ~ 1,
        obs_vals == pred_min ~ -1,
        obs_vals == pred_max2 ~ 2,
        obs_vals == pred_min2 ~ -2,
        obs_vals > pred_med & obs_vals < pred_max ~ (obs_vals - pred_med) / (pred_max - pred_med),
        obs_vals < pred_med & obs_vals > pred_min ~ (obs_vals - pred_med) / (pred_med - pred_min),
        obs_vals > pred_max & obs_vals < pred_max2 ~ (obs_vals - pred_max) / (pred_max2 - pred_max) + 1,
        obs_vals < pred_min & obs_vals > pred_min2 ~ (obs_vals - pred_min) / (pred_min - pred_min2) - 1,
        obs_vals > pred_max2 ~ (obs_vals - pred_max2) / (obs_errs * 1.96) + 2,
        obs_vals < pred_min2 ~ (obs_vals - pred_min2) / (obs_errs * 1.96) - 2
      )
    ) %>%
    ungroup()

  # tidy up other_df
  other_df <- bind_rows(other_df)

  # tidy up residual_df
  residual_df <- bind_rows(residual_df) %>%
    mutate(
      obs_min = obs_vals - obs_errs * 1.96,
      obs_max = obs_vals + obs_errs * 1.96
    ) %>%
    group_by(var_name) %>%
    mutate(
      all_min = min(obs_min, pred_min, pred_min2),
      all_max = max(obs_max, pred_max, pred_max2),
      score = case_when(
        obs_vals == pred_med ~ 0,
        obs_vals == pred_max ~ 1,
        obs_vals == pred_min ~ -1,
        obs_vals == pred_max2 ~ 2,
        obs_vals == pred_min2 ~ -2,
        obs_vals > pred_med & obs_vals < pred_max ~ (obs_vals - pred_med) / (pred_max - pred_med),
        obs_vals < pred_med & obs_vals > pred_min ~ (obs_vals - pred_med) / (pred_med - pred_min),
        obs_vals > pred_max & obs_vals < pred_max2 ~ (obs_vals - pred_max) / (pred_max2 - pred_max) + 1,
        obs_vals < pred_min & obs_vals > pred_min2 ~ (obs_vals - pred_min) / (pred_min - pred_min2) - 1,
        obs_vals > pred_max2 ~ (obs_vals - pred_max2) / (obs_errs * 1.96) + 2,
        obs_vals < pred_min2 ~ (obs_vals - pred_min2) / (obs_errs * 1.96) - 2
      )
    ) %>%
    ungroup()
}

# plot results using ggplot ####
if (TRUE) {
  cat(file = stderr(), "Plot results using ggplot", "\n")

  # tans
  xpale <- "#FFF7BC" # SRON Fig 7 row 1
  xlight <- "#FEC44F" # SRON Fig 7 row 1
  xmid <- "#D95F0E" # SRON Fig 7 row 1
  xdark <- "#993404" # SRON Fig 7 row 3
  xdata <- "darkblue"
  xaxis <- "#999999"

  # greens 2 # https://www.color-hex.com/color-palette/5016
  xpale <- "#F0F7DA"
  xlight <- "#C9DF8A"
  xmid <- "#77AB59"
  xdark <- "#36802D"
  xdata <- "#043927" # https://graf1x.com/shades-of-green-color-palette-html-hex-rgb-code/
  xaxis <- "#999999"
  xgrey <- "grey"

  # greens 3 # https://graf1x.com/shades-of-green-color-palette-html-hex-rgb-code/

  # calibration
  cat(file = stderr(), "Plotting calibration variables against data", "\n")
  limits <- NULL
  plot3 <- result_df %>%
    filter(var_name0 %in% calibOutputs) %>% 
    ggplot() +
    labs(title = "", x = "", y = "") +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min2, ymax = pred_max2), fill = xpale, colour = xlight, size = 0.5) +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min, ymax = pred_max), fill = xlight, colour = xlight, size = 0.5) +
    geom_line(mapping = aes(x = times, y = pred_med), colour = xmid, size = 0.5, na.rm = TRUE) +
    geom_line(mapping = aes(x = times, y = pred_map), colour = xdark, size = 0.5, na.rm = TRUE) +
    # geom_point(mapping=aes(x=times, y=obs_vals, size=-logl), colour=xdata) +
    geom_point(mapping = aes(x = times, y = obs_vals, colour = factor(obs_wts)), size = 2, na.rm = TRUE) +
    # geom_point(mapping=aes(x=times, y=obs_vals, colour=site_name)) +
    # scale_color_manual(values=c(blue, green, red)) +
    geom_abline(mapping = aes(slope = 0, intercept = 0), colour = xaxis, size = 0.5) +
    # geom_abline(mapping = aes(slope = 0, intercept = pred_med_mean), colour = "red", size = 0.5) +
    scale_y_continuous(expand = expand_scale(0, 0)) +
    scale_x_continuous(expand = expand_scale(0, 0), limits = limits) +
    scale_colour_manual(values = c("1" = xdata, "0" = "grey")) +
    facet_grid(var_name2 ~ site_name, scales = "free_y") +
    guides(colour = FALSE, size = FALSE) +
    theme_few() +
    theme(panel.grid.major.x = element_line(colour = xgrey))
  # print(plot3)
  file_name <- paste(scenario, "/calibration_time.png", sep = "")
  png(file_name, width = 210, height = 297 / 7 * length(calibOutputs), units = "mm", type = "windows", res = 600)
  print(plot3)
  dev.off()

  # calibration diff
  cat(file = stderr(), "Plotting calibration variables differences", "\n")
  limits <- NULL
  temp <- result_df %>%
    filter(var_name0 %in% calibOutputs) %>% 
    group_by(var_name, times) %>%
    arrange(site_name) %>% 
    mutate(
      pred_min2 = pred_min2 - pred_min2[base_s],
      pred_max2 = pred_max2 - pred_max2[base_s],
      pred_min = pred_min - pred_min[base_s],
      pred_max = pred_max - pred_max[base_s],
      pred_med = pred_med - pred_med[base_s],
      pred_med_mean = pred_med_mean - pred_med_mean[base_s],
      pred_map = pred_map - pred_map[base_s],
      obs_vals = obs_vals - obs_vals[base_s],
      n = n()
    ) %>% 
    ungroup()
  plot3 <- temp %>%
    ggplot() +
    labs(title = "", x = "", y = "") +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min2, ymax = pred_max2), fill = xpale, colour = xlight, size = 0.5) +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min, ymax = pred_max), fill = xlight, colour = xlight, size = 0.5) +
    geom_line(mapping = aes(x = times, y = pred_med), colour = xmid, size = 0.5, na.rm = TRUE) +
    geom_line(mapping = aes(x = times, y = pred_map), colour = xdark, size = 0.5, na.rm = TRUE) +
    # geom_point(mapping=aes(x=times, y=obs_vals, size=-logl), colour=xdata) +
    geom_point(mapping = aes(x = times, y = obs_vals, colour = factor(obs_wts)), size = 2, na.rm = TRUE) +
    # geom_point(mapping=aes(x=times, y=obs_vals, colour=site_name)) +
    # scale_color_manual(values=c(blue, green, red)) +
    geom_abline(mapping = aes(slope = 0, intercept = 0), colour = xaxis, size = 0.5) +
    geom_abline(mapping = aes(slope = 0, intercept = pred_med_mean), colour = "red", size = 0.5) +
    scale_y_continuous(expand = expand_scale(0, 0)) +
    scale_x_continuous(expand = expand_scale(0, 0), limits = limits) +
    scale_colour_manual(values = c("1" = xdata, "0" = "grey")) +
    facet_grid(var_name2 ~ site_name, scales = "free_y") +
    guides(colour = FALSE, size = FALSE) +
    theme_few() +
    theme(panel.grid.major.x = element_line(colour = xgrey))
  # print(plot3)
  file_name <- paste(scenario, "/calibration_time_diff.png", sep = "")
  png(file_name, width = 210, height = 297 / 7 * length(calibOutputs), units = "mm", type = "windows", res = 600)
  print(plot3)
  dev.off()
  
  # calibration single year
  cat(file = stderr(), "Plotting calibration variables single year", "\n")
  limits <- single_year + c(152 / 365, 1 + 151 / 365)
  plot3 <- result_df %>%
    filter(var_name0 %in% calibOutputs) %>% 
    ggplot() +
    labs(title = "", x = "", y = "") +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min2, ymax = pred_max2), fill = xpale, colour = xlight, size = 0.5) +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min, ymax = pred_max), fill = xlight, colour = xlight, size = 0.5) +
    geom_line(mapping = aes(x = times, y = pred_med), colour = xmid, size = 0.5, na.rm = TRUE) +
    geom_line(mapping = aes(x = times, y = pred_map), colour = xdark, size = 0.5, na.rm = TRUE) +
    # geom_point(mapping=aes(x=times, y=obs_vals, size=-logl), colour=xdata) +
    geom_point(mapping = aes(x = times, y = obs_vals, colour = factor(obs_wts)), size = 2, na.rm = TRUE) +
    # geom_point(mapping=aes(x=times, y=obs_vals, colour=site_name)) +
    # scale_color_manual(values=c(blue, green, red)) +
    geom_abline(mapping = aes(slope = 0, intercept = 0), colour = xaxis, size = 0.5) +
    scale_y_continuous(expand = expand_scale(c(0, 0.05), 0)) +
    scale_x_continuous(expand = expand_scale(0, 0), limits = limits) +
    scale_colour_manual(values = c("1" = xdata, "0" = "grey")) +
    facet_grid(var_name2 ~ site_name, scales = "free_y") +
    guides(colour = FALSE, size = FALSE) +
    theme_few()
  # print(plot3)
  file_name <- paste(scenario, "/calibration_", single_year, ".png", sep = "")
  png(file_name, width = 210, height = 297 / 7 * length(calibOutputs), units = "mm", type = "windows", res = 600)
  print(plot3)
  dev.off()
  
  # calibration single year diff
  cat(file = stderr(), "Plotting calibration variables single year differences", "\n")
  limits <- single_year + c(152 / 365, 1 + 151 / 365)
  plot3 <- temp %>%
    ggplot() +
    labs(title = "", x = "", y = "") +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min2, ymax = pred_max2), fill = xpale, colour = xlight, size = 0.5) +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min, ymax = pred_max), fill = xlight, colour = xlight, size = 0.5) +
    geom_line(mapping = aes(x = times, y = pred_med), colour = xmid, size = 0.5, na.rm = TRUE) +
    geom_line(mapping = aes(x = times, y = pred_map), colour = xdark, size = 0.5, na.rm = TRUE) +
    # geom_point(mapping=aes(x=times, y=obs_vals, size=-logl), colour=xdata) +
    geom_point(mapping = aes(x = times, y = obs_vals, colour = factor(obs_wts)), size = 2, na.rm = TRUE) +
    # geom_point(mapping=aes(x=times, y=obs_vals, colour=site_name)) +
    # scale_color_manual(values=c(blue, green, red)) +
    geom_abline(mapping = aes(slope = 0, intercept = 0), colour = xaxis, size = 0.5) +
    scale_y_continuous(expand = expand_scale(c(0, 0.05), 0)) +
    scale_x_continuous(expand = expand_scale(0, 0), limits = limits) +
    scale_colour_manual(values = c("1" = xdata, "0" = "grey")) +
    facet_grid(var_name2 ~ site_name, scales = "free_y") +
    guides(colour = FALSE, size = FALSE) +
    theme_few()
  # print(plot3)
  file_name <- paste(scenario, "/calibration_", single_year, "_diff.png", sep = "")
  png(file_name, width = 210, height = 297 / 7 * length(calibOutputs), units = "mm", type = "windows", res = 600)
  print(plot3)
  dev.off()
  
  # other
  cat(file = stderr(), "Plotting other variables ", "\n")
  limits <- NULL
  plot3 <- other_df %>%
    filter(var_name0 %in% otherOutputs) %>% 
    ggplot() +
    labs(title = "", x = "", y = "") +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min2, ymax = pred_max2), fill = xpale, colour = xlight, size = 0.5) +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min, ymax = pred_max), fill = xlight, colour = xlight, size = 0.5) +
    geom_line(mapping = aes(x = times, y = pred_med), colour = xmid, size = 0.5, na.rm = TRUE) +
    geom_line(mapping = aes(x = times, y = pred_map), colour = xdark, size = 0.5, na.rm = TRUE) +
    # geom_point(mapping=aes(x=times, y=obs_vals), colour="xdata") +
    geom_abline(mapping = aes(slope = 0, intercept = 0), colour = xaxis, size = 0.5) +
    # geom_abline(mapping = aes(slope = 0, intercept = pred_med_mean), colour = "red", size = 0.5) +
    scale_y_continuous(expand = expand_scale(c(0, 0.05), 0)) +
    scale_x_continuous(expand = expand_scale(0, 0), limits = limits) +
    facet_grid(var_name2 ~ site_name, scales = "free_y") +
    theme_few() +
    theme(panel.grid.major.x = element_line(colour = xgrey))
  # print(plot3)
  file_name <- paste(scenario, "/other_time.png", sep = "")
  png(file_name, width = 210, height = 297 / 7 * length(otherOutputs), units = "mm", type = "windows", res = 600)
  print(plot3)
  dev.off()
  
  # other diff
  cat(file = stderr(), "Plotting other variables differences", "\n")
  limits <- NULL
  temp <- other_df %>%
    filter(var_name0 %in% otherOutputs) %>% 
    group_by(var_name, times) %>%
    arrange(site_name) %>% 
    mutate(
      pred_min2 = pred_min2 - pred_min2[base_s],
      pred_max2 = pred_max2 - pred_max2[base_s],
      pred_min = pred_min - pred_min[base_s],
      pred_max = pred_max - pred_max[base_s],
      pred_med = pred_med - pred_med[base_s],
      pred_med_mean = pred_med_mean - pred_med_mean[base_s],
      pred_map = pred_map - pred_map[base_s],
      # obs_vals = obs_vals - obs_vals[base_s],
      n = n()
    ) %>% 
    ungroup()
  plot3 <- temp %>%
    ggplot() +
    labs(title = "", x = "", y = "") +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min2, ymax = pred_max2), fill = xpale, colour = xlight, size = 0.5) +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min, ymax = pred_max), fill = xlight, colour = xlight, size = 0.5) +
    geom_line(mapping = aes(x = times, y = pred_med), colour = xmid, size = 0.5, na.rm = TRUE) +
    geom_line(mapping = aes(x = times, y = pred_map), colour = xdark, size = 0.5, na.rm = TRUE) +
    # geom_point(mapping=aes(x=times, y=obs_vals), colour="xdata") +
    geom_abline(mapping = aes(slope = 0, intercept = 0), colour = xaxis, size = 0.5) +
    geom_abline(mapping = aes(slope = 0, intercept = pred_med_mean), colour = "red", size = 0.5) +
    scale_y_continuous(expand = expand_scale(c(0, 0.05), 0)) +
    scale_x_continuous(expand = expand_scale(0, 0), limits = limits) +
    facet_grid(var_name2 ~ site_name, scales = "free_y") +
    theme_few() +
    theme(panel.grid.major.x = element_line(colour = xgrey))
  # print(plot3)
  file_name <- paste(scenario, "/other_time_diff.png", sep = "")
  png(file_name, width = 210, height = 297 / 7 * length(otherOutputs), units = "mm", type = "windows", res = 600)
  print(plot3)
  dev.off()
  
  # other single year
  cat(file = stderr(), "Plotting other variables single year", "\n")
  limits <- single_year + c(152 / 365, 1 + 151 / 365)
  plot3 <- other_df %>%
    filter(var_name0 %in% otherOutputs) %>% 
    ggplot() +
    labs(title = "", x = "", y = "") +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min2, ymax = pred_max2), fill = xpale, colour = xlight, size = 0.5) +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min, ymax = pred_max), fill = xlight, colour = xlight, size = 0.5) +
    geom_line(mapping = aes(x = times, y = pred_med), colour = xmid, size = 0.5, na.rm = TRUE) +
    geom_line(mapping = aes(x = times, y = pred_map), colour = xdark, size = 0.5, na.rm = TRUE) +
    # geom_point(mapping=aes(x=times, y=obs_vals), colour="xdata") +
    geom_abline(mapping = aes(slope = 0, intercept = 0), colour = xaxis, size = 0.5) +
    scale_y_continuous(expand = expand_scale(c(0, 0.05), 0)) +
    scale_x_continuous(expand = expand_scale(0, 0), limits = limits) +
    facet_grid(var_name2 ~ site_name, scales = "free_y") +
    theme_few()
  # print(plot3)
  file_name <- paste(scenario, "/other_", single_year, ".png", sep = "")
  png(file_name, width = 210, height = 297 / 7 * length(otherOutputs), units = "mm", type = "windows", res = 600)
  print(plot3)
  dev.off()
  
  # other single year diff
  cat(file = stderr(), "Plotting other variables single year differences", "\n")
  limits <- single_year + c(152 / 365, 1 + 151 / 365)
  plot3 <- temp %>%
    ggplot() +
    labs(title = "", x = "", y = "") +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min2, ymax = pred_max2), fill = xpale, colour = xlight, size = 0.5) +
    geom_ribbon(mapping = aes(x = times, ymin = pred_min, ymax = pred_max), fill = xlight, colour = xlight, size = 0.5) +
    geom_line(mapping = aes(x = times, y = pred_med), colour = xmid, size = 0.5, na.rm = TRUE) +
    geom_line(mapping = aes(x = times, y = pred_map), colour = xdark, size = 0.5, na.rm = TRUE) +
    # geom_point(mapping=aes(x=times, y=obs_vals), colour="xdata") +
    geom_abline(mapping = aes(slope = 0, intercept = 0), colour = xaxis, size = 0.5) +
    scale_y_continuous(expand = expand_scale(c(0, 0.05), 0)) +
    scale_x_continuous(expand = expand_scale(0, 0), limits = limits) +
    facet_grid(var_name2 ~ site_name, scales = "free_y") +
    theme_few()
  # print(plot3)
  file_name <- paste(scenario, "/other_", single_year, "_diff.png", sep = "")
  png(file_name, width = 210, height = 297 / 7 * length(otherOutputs), units = "mm", type = "windows", res = 600)
  print(plot3)
  dev.off()
  
  # set site colours
  palette <- hue_pal()(nSites)
  palette[1:3] <- rev(palette[1:3]) # nicer for Northland, Waikato, Canterbury
  show_col(palette)
  
  # residuals (data - model)/error
  cat(file = stderr(), "Plotting residuals density", "\n")
  temp <- filter(residual_df, obs_wts > 0, pred_map > 0 | obs_vals > 0) # avoid Stem C data and model == 0
  plot2 <- temp %>%
    ggplot() +
    labs(title = "", x = "", y = "", colour = "Run") +
    # geom_density(mapping=aes(x=(pred_map-obs_vals)/obs_errs, colour=site_name)) +
    stat_density(mapping = aes(x = -(pred_map - obs_vals) / obs_errs, colour = site_name), position = "identity", fill = NA, adjust = 1) +
    geom_vline(mapping = aes(xintercept = 0), colour = "black") +
    theme_few() +
    scale_y_continuous(expand = expand_scale(c(0, 0.05), 0)) +
    scale_x_continuous(limits = c(-5, 5), expand = expand_scale(0, 0)) +
    scale_color_manual(values=palette) +
    scale_fill_manual(values=palette) +
    facet_wrap(~var_name, scale = "free")
  # print(plot2)
  png(paste(scenario, "/residual_density.png", sep = ""), width = 11, height = 8, units = "in", type = "windows", res = 300)
  print(plot2)
  dev.off()

  # residuals bias and precision
  cat(file = stderr(), "Plotting residuals bias and precision", "\n")
  temp <- filter(
    residual_df,
    # pred_map>0 | obs_vals>0, # avoid Stem C data and model == 0
    obs_wts > 0
  ) %>%
    mutate(xjitter = runif(n()) * 0.1 - 0.05)
  plot3 <- temp %>%
    ggplot() +
    labs(title = "", x = "", y = "", colour = "Run", fill = "Run") +
    # geom_ribbon(mapping=aes(x=times+xjitter, ymin=-obs_errs*1.96, ymax=obs_errs*1.96), fill="lightgrey") +
    geom_ribbon(mapping = aes(x = times + xjitter, ymin = pred_min2 - pred_med, ymax = pred_max2 - pred_med, fill = site_name, colour = site_name), alpha = 0.1) +
    geom_ribbon(mapping = aes(x = times + xjitter, ymin = pred_min - pred_med, ymax = pred_max - pred_med, fill = site_name, colour = site_name), alpha = 0.3) +
    # geom_errorbar(mapping=aes(x=times+xjitter, ymin=pred_min2-obs_vals, ymax=pred_max2-obs_vals), colour="grey", width=0.1) +
    # geom_errorbar(mapping=aes(x=times+xjitter, ymin=pred_min2-pred_min, ymax=pred_max2-pred_max), colour="grey", width=0.1) +
    # geom_errorbar(mapping=aes(x=times+xjitter, ymin=obs_min-pred_med, ymax=obs_max-pred_med, colour=site_name)) +
    # geom_errorbarh(mapping=aes(y=obs_vals, xmin=pred_min, xmax=pred_max, colour=site_name)) +
    geom_point(mapping = aes(x = times + xjitter, y = obs_vals - pred_med, colour = site_name), na.rm = TRUE) +
    geom_abline(mapping = aes(slope = 0, intercept = 0), colour = "black") +
    # geom_smooth(mapping=aes(x=times, y=resid_mean, colour=site_name), method="lm", se=FALSE) +
    scale_color_manual(values=palette) +
    scale_fill_manual(values=palette) +
    # scale_y_continuous(expand=expand_scale(0,0)) +
    # scale_x_continuous(expand=expand_scale(0,0)) +
    # facet_wrap(~var_name, scale="free_y") +
    lemon::facet_rep_wrap(~var_name2, scales = "free_y", repeat.tick.labels = "bottom") +
    theme_few()
  # print(plot3)
  png(paste(scenario, "/residual_time.png", sep = ""), width = 11, height = 8, units = "in", type = "windows", res = 300)
  print(plot3)
  dev.off()

  # residuals bias and precision #2
  cat(file = stderr(), "Plotting residuals bias and precision 2", "\n")
  temp <- filter(
    residual_df,
    # pred_map>0 | obs_vals>0, # avoid Stem C data and model == 0
    obs_wts > 0
  ) %>%
    mutate(xjitter = runif(n()) * 0.1 - 0.05)
  plot3b <- temp %>%
    ggplot() +
    labs(title = "", x = "", y = "", colour = "Run", fill = "Run") +
    # geom_ribbon(mapping=aes(x=pred_med, ymin=-obs_errs*1.96, ymax=obs_errs*1.96), fill="lightgrey") +
    geom_ribbon(mapping = aes(x = pred_med, ymin = pred_min2, ymax = pred_max2, fill = site_name, colour = site_name), alpha = 0.1) +
    geom_ribbon(mapping = aes(x = pred_med, ymin = pred_min, ymax = pred_max, fill = site_name, colour = site_name), alpha = 0.3) +
    # geom_errorbar(mapping=aes(x=pred_med, ymin=pred_min2-obs_vals, ymax=pred_max2-obs_vals), colour="grey", width=0.1) +
    # geom_errorbar(mapping=aes(x=pred_med, ymin=pred_min2-pred_min, ymax=pred_max2-pred_max), colour="grey", width=0.1) +
    # geom_errorbar(mapping=aes(x=pred_med, ymin=obs_min-pred_med, ymax=obs_max-pred_med, colour=site_name)) +
    # geom_errorbarh(mapping=aes(y=pred_med, xmin=pred_min, xmax=pred_max, colour=site_name)) +
    geom_point(mapping = aes(x = pred_med, y = obs_vals, colour = site_name), na.rm = TRUE) +
    geom_abline(mapping = aes(slope = 1, intercept = 0), colour = "black") +
    # geom_smooth(mapping=aes(x=times, y=resid_mean, colour=site_name), method="lm", se=FALSE) +
    scale_color_manual(values=palette) +
    scale_fill_manual(values=palette) +
    scale_y_continuous(expand = expand_scale(0, 0)) +
    scale_x_continuous(expand = expand_scale(0, 0)) +
    facet_wrap(~var_name2, scale = "free") +
    theme_few()
  # print(plot3b)
  png(paste(scenario, "/residual_data2.png", sep = ""), width = 11, height = 8, units = "in", type = "windows", res = 300)
  print(plot3b)
  dev.off()
}

# prior and posterior histograms ####
if (TRUE) {
  cat(file = stderr(), "Plot prior/posterior histograms", "\n")
  post_df <- as.data.frame(pChain %*% diag(sc))
  names(post_df) <- colnames(pChain)
  post_df <- gather(post_df)
  prior_df <- vector("list", nBCvar)
  my_pretty_breaks <- function(n = 5, ...) {
    n_default <- n
    function(x, n = n_default) {
      minx <- min(x)
      maxx <- max(x)
      midx <- (minx + maxx) / 2
      x2 <- midx + (x - midx) * 0.9
      breaks <- pretty(x2, n, ...)
      names(breaks) <- attr(breaks, "labels")
      breaks
    }
  }
  i <- 1
  for (i in seq_along(parname_BC)) {
    key2 <- colnames(pChain)[i]
    x <- seq(parmin_BC[i], parmax_BC[i], (parmax_BC[i] - parmin_BC[i]) / 100)
    y <- dbeta((x - parmin_BC[i]) / (parmax_BC[i] - parmin_BC[i]), aa[i], bb[i]) / (parmax_BC[i] - parmin_BC[i])
    xmap <- scparMAP_BC[i] * sc[i]
    # ymap <- dbeta((xmap - parmin_BC[i]) / (parmax_BC[i] - parmin_BC[i]), aa[i], bb[i]) / (parmax_BC[i] - parmin_BC[i])
    temp <- ggplot() + 
      geom_histogram(data = post_df %>% filter(key == key2),
                     mapping = aes(x = value, y = ..density..), bins = 30) +
      scale_x_continuous(expand = expand_scale(0, 0), breaks = my_pretty_breaks(n = 2, min.n = 2))
    temp <- layer_data(temp)
    j <- which(abs(temp$x - xmap) == min(abs(temp$x - xmap)))
    ymap <- temp$density[j]
    prior_df[[i]] <- tibble(key = key2, x = x, y = y, xmap = xmap, ymap = ymap)
  }
  prior_df <- bind_rows(prior_df) 
  post_df <- post_df %>% filter(str_detect(key, "I\\([1-9]:[1-9]\\)") == FALSE)
  prior_df <- prior_df %>% filter(str_detect(key, "I\\([1-9]:[1-9]\\)") == FALSE)
  plot1 <- ggplot(data = prior_df) +
    labs(title = "", x = "", y = "") +
    geom_line(mapping = aes(x = x, y = y), colour = xmid, size = 1) +
    geom_histogram(data = post_df, mapping = aes(x = value, y = ..density..), fill = xlight, bins = 30) +
    geom_line(mapping = aes(x = x, y = y), colour = xmid, size = 1) +
    # geom_point(mapping = aes(x = xmap, y = ymap), colour = xdata, size = 3, na.rm = TRUE) +
    geom_segment(mapping = aes(x = xmap, xend = xmap, y = 0, yend = ymap), colour = xdata, size = 1) +
    theme_few() +
    theme(panel.spacing.x = unit(9, "mm"), plot.margin = unit(c(0, 5, 0, 0), "mm")) +
    # theme(axis.text=element_text(size=9)) +
    scale_x_continuous(expand = expand_scale(0, 0), breaks = my_pretty_breaks(n = 2, min.n = 2)) +
    scale_y_continuous(expand = expand_scale(0, 0), breaks = NULL) +
    facet_wrap(vars(key), scales = "free")
  # print(plot1)
  png(paste(scenario, "/BC_parameter_histograms_BT.png", sep = ""),
    width = 210 * 1.5, height = 210 * 1.5, units = "mm",
    type = "windows", res = 600
  )
  print(plot1)
  dev.off()
}

# memory management
cat(file = stderr(), "Saving checkpoint after BASGRA results", "\n")
save.image(file = paste(scenario, "/checkpoint_after_results.RData", sep = ""))
