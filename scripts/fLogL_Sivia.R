flogL <- function(sims,data,data_s,data_w=1)
{ 
Ri         <- (sims - data) / data_s
#i0         <- which( abs(Ri)<1.e-08 ) # catch tiny values which cause problems
logLi      <- log(1-exp(-0.5*Ri^2)) - 2*log(abs(Ri)) - 0.5*log(2*pi) - log(data_s)
i0 <- which(!is.finite(logLi)) # Simon this is more robust I think
logLi[i0]  <- -0.5*log(2*pi) - log(2*data_s[i0])

sum(logLi*data_w)
}

#   logLi <- -0.5*z^2 - 0.5*log(2*pi) - log(data_s) # Gaussian
#   logLi <- log(1-exp(-0.5*z^2)) - 2*log(abs(z)) - 0.5*log(2*pi) - log(data_s) # Sivia
#   logLi <- - log(1+z^2) - log(pi) # Cauchy
#   logLi <- -((7+1)/2)*log(1+z^2/7) - 0.5*log(7*pi) + log(gamma((7+1)/2)/gamma(7/2)) # Students7
