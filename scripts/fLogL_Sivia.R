flogL <- function(sims,data,data_s,data_w=1)
{ 
Ri         <- (sims - data) / data_s
#i0         <- which( abs(Ri)<1.e-08 ) # catch tiny values which cause problems
logLi      <- log(1-exp(-0.5*Ri^2)) - 2*log(abs(Ri)) - 0.5*log(2*pi) - log(data_s)
i0 <- which(!is.finite(logLi)) # Simon this is more robust I think
logLi[i0]  <- -0.5*log(2*pi) - log(2*data_s[i0])

sum(logLi*data_w)
}
