flogL <- function(sims,data,data_s,data_w=1)
{ 
  Ri         <- (sims - data) / data_s
  logLi      <- log(1-exp(-0.5*Ri^2)) - 2*log(abs(Ri)) - 0.5*log(2*pi) - log(data_s)
  i0         <- which(!is.finite(logLi)) # Simon this is more robust I think
  logLi[i0]  <- -0.5*log(2*pi) - log(2*data_s[i0])
  
  sum(logLi*data_w)
}

flogLi <- function(sims,data,data_s,data_w=1)
{ 
  Ri         <- (sims - data) / data_s
  logLi      <- log(1-exp(-0.5*Ri^2)) - 2*log(abs(Ri)) - 0.5*log(2*pi) - log(data_s)
  i0         <- which(!is.finite(logLi)) # Simon this is more robust I think
  logLi[i0]  <- -0.5*log(2*pi) - log(2*data_s[i0])
  
  return(logLi*data_w)
}
