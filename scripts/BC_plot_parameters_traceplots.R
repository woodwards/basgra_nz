## Preparing for plotting ##
 # Read the parameter names
   titles     <- parname_BC
   
 # We will write the plots to a pdf file
   pagew <- 11 ; pageh <- 8
   # png( paste("model_outputs/BC_parameters_traceplots",format(Sys.time(),"_%H_%M.png"),sep=""),
   #      width=pagew, height=pageh, units="in", type="windows", res=300)
   png( paste("model_outputs/BC_parameters_convergence.png",sep=""),
        width=pagew, height=pageh, units="in", type="windows", res=300)
   
## Parameter trace plots ##
   nrowsPlots <- ceiling( sqrt((np_BC+1)*pageh/pagew) )
   ncolsPlots <- ceiling( (np_BC+1)/nrowsPlots )
   par( mfrow = c(nrowsPlots,ncolsPlots) )
   par(mar=c(2, 2, 2, 1))
   for (i in seq(1,np_BC)){
        plot( pChain[,i] * sc[i],
		      type='l', xlab="", ylab="", col=rgb(0,0,0,0.3), 
          main=paste(titles[i],parsites_BC[i]), cex.main=1 )
        abline( v=nBurnin, col='red', lwd=2, lty=2 )
   }
   plot(1,type='n', axes=FALSE, xlab="", ylab="")
   plot_colors <- c("black","red")
   legend("bottomright", c("Parameter trace","End of burn-in"),
          bty="n", lty=c(1,2), lwd=c(1,2), col=plot_colors, title = "LEGEND:")


## Closing ##
   dev.off()
