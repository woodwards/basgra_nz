## Preparing for plotting ##
 # Read the parameter names
  titles <- ifelse(parsites_BC=="1:nSites", 
                   as.character(parname_BC), 
                   paste(parname_BC, "(", parsites_BC, ")", sep=""))

 # We will write the plots to a pdf file:
   pagew <- 11 ; pageh <- 8
   # png( paste("model_outputs/BC_parameters_priorbeta_histograms",format(Sys.time(),"_%H_%M.png"),sep=""),
   #      width=pagew, height=pageh, units="in", type="windows", res=300)
   png( paste(scenario, "/BC_parameter_histograms.png",sep=""),
        width=pagew, height=pageh, units="in", type="windows", res=300)
   
## Histogram plots ##
   nrowsPlots <- ceiling( sqrt((np_BC+1)*pageh/pagew) )
   ncolsPlots <- ceiling( (np_BC+1)/nrowsPlots )
   par( mfrow = c(nrowsPlots,ncolsPlots) )
   par( mar   =c(2, 2, 2, 1))
   nbreaks <- 20
   # keeps <- seq(nBurnin+1, nChain, length.out=min(nChain-nBurnin,100000)) # avoid overflow error from large samples
   keeps <- seq(1, nSampling, length.out=min(nSampling,100000)) # avoid overflow error from large samples
   for (i in seq(1,np_BC)){
        hist( pChain[keeps,i] * sc[i],
              xlab="", ylab="", main=titles[i], cex.main=1,
              breaks=nbreaks, freq=FALSE, xlim=c(parmin_BC[i],parmax_BC[i]) )
        parseq_BC <- seq(parmin_BC[i],parmax_BC[i],(parmax_BC[i]-parmin_BC[i])/100)
        points( parseq_BC, dbeta( (parseq_BC-parmin_BC[i])/(parmax_BC[i]-parmin_BC[i]), aa[i], bb[i] ) /
                           (parmax_BC[i]-parmin_BC[i]),
                type='l',col='red')
        points( scparMAP_BC[i]*sc[i], dbeta( (scparMAP_BC[i]*sc[i]-parmin_BC[i])/(parmax_BC[i]-parmin_BC[i]),
                                             aa[i], bb[i] ) / (parmax_BC[i]-parmin_BC[i]),
                pch=16, cex=2, col='blue') # Simon add MAP point
   }
   plot(1,type='n', axes=FALSE, xlab="", ylab="")
   plot_colors <- c("red","black","blue")
   legend("bottomright", c("Prior","Post.","MAP"),
          bty="n", lty=1, lwd=3, col=plot_colors, title = "LEGEND:")

## Closing ##
   dev.off()
