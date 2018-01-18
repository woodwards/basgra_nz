params_BC_ModePrior <-   parmod_BC
params_BC_MAP       <- scparMAP_BC  * sc
params_BC_MaxL      <- scparMaxL_BC * sc

dev.set()
pagew <- 11 ; pageh <- 8
# png( paste('model_outputs/BC_outputs_data',format(Sys.time(),"_%H_%M.png"),sep=""),
#      width=pagew, height=pageh, units="in", type="windows", res=300)

s <- 1
for (s in 1:nSites) {
  
  png( paste('model_outputs/BC_calibration_fits_', s, '.png',sep=""),
       width=pagew, height=pageh, units="in", type="windows", res=300)
  
  params          <- list_params      [[s]] ; matrix_weather <- list_matrix_weather[[s]]
  days_harvest    <- list_days_harvest[[s]] ; NDAYS          <- list_NDAYS         [[s]]
  ip_BC_s         <- ip_BC_site       [[s]]
  icol_pChain_s   <- icol_pChain_site [[s]]
# Calculate model output for the prior mode
  params[ip_BC_s] <- params_BC_ModePrior[icol_pChain_s]
  outputPriorMode <- run_model(params,matrix_weather,days_harvest,NDAYS)
# Calculate model output for the MAP parameter vector
  params[ip_BC_s] <- params_BC_MAP      [icol_pChain_s]
  outputMAP       <- run_model(params,matrix_weather,days_harvest,NDAYS)
# Calculate model output for the MaxL parameter vector
  params[ip_BC_s] <- params_BC_MaxL     [icol_pChain_s]
  outputMaxL      <- run_model(params,matrix_weather,days_harvest,NDAYS)
# Calculate model output for a sample from the posterior
# Take a sample (of size nSample) from the chain generated using MCMC
  nSample         <- 1000
  nStep           <- (nChain-nBurnin) / nSample
  outputSample    <- array( 0, c(nSample,NDAYS,NOUT) )
  print(paste("Running", nSample, "posteriors of", nChain-nBurnin, ", site", s, "of", nSites))
  ii              <- 0   
  jj              <- as.integer(seq(from=nBurnin+nStep, to=nChain, length.out=nSample))
  for (j in jj) {
    ii <- ii+1
    params_j           <- pChain[j,] * sc
    params[ip_BC_s]    <- params_j[icol_pChain_s]
    outputSample[ii,,] <- run_model(params,matrix_weather,days_harvest,NDAYS)
  } # end of sample loop
# Analyse the posterior output sample: calculate quantiles 5% and 95%. Simon added na.rm=TRUE
  q5  <- sapply( 1:NOUT, function(i) sapply(1:NDAYS,function(j)quantile(outputSample[,j,i],0.05,na.rm=TRUE)) ) 
  q95 <- sapply( 1:NOUT, function(i) sapply(1:NDAYS,function(j)quantile(outputSample[,j,i],0.95,na.rm=TRUE)) )

# Plot
  list_runs <- list( outputPriorMode, outputMaxL, outputMAP, q5, q95 )
  isite       = s
  list_runs   = list_runs
  nruns        = length(list_runs)
  leg_title   = "BC"
  leg         = c("Prior","MaxL","MAP", "q5", "q95")
  cols        = c("lightgrey",  "lightpink", "firebrick3", "firebrick3", "firebrick3" )
  lwds        = c( 2, 2, 2, 1, 1 )
  ltys        = c( 1, 1, 1, 7, 7 ) # ltys==7 in 2 series draws a ploygon
  col01       = c('darkgrey', 'dodgerblue3')
  plot_outputs_data_s( isite       = s,
                       list_runs   = list_runs,
                       leg_title   = "BC",
                       leg         = c("Prior","MaxL","MAP", "q5", "q95"),
                       cols        = c("lightgrey",  "lightpink", "firebrick3", "firebrick3", "firebrick3" ),
                       lwds        = c( 2, 2, 2, 1, 1 ),
                       ltys        = c( 1, 1, 1, 7, 7 ), # ltys==7 in 2 series draws a ploygon
                       col01       = c('darkgrey', 'dodgerblue3')
                       )   
dev.off()   
}