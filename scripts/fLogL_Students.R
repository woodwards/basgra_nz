flogL <- function(sims,data,data_s,data_w=1,dof=7)
{ 
Ri         <- (sims - data) / data_s
logLi      <- - ((dof+1)/2)*log(1+(Ri^2)/dof) - 0.5*log(dof*pi) + log(gamma((dof+1)/2)/gamma(dof/2))
i0         <- which(!is.finite(logLi))
logLi[i0]  <- 0

sum(logLi*data_w)
}

#   logLi <- -0.5*z^2 - 0.5*log(2*pi) - log(data_s) # Gaussian
#   logLi <- log(1-exp(-0.5*z^2)) - 2*log(abs(z)) - 0.5*log(2*pi) - log(data_s) # Sivia
#   logLi <- - log(1+z^2) - log(pi) # Cauchy
#   logLi <- -((7+1)/2)*log(1+z^2/7) - 0.5*log(7*pi) + log(gamma((7+1)/2)/gamma(7/2)) # Students7
