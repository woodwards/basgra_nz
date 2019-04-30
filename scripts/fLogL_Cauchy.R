flogL <- function(sims,data,data_s,data_w=1)
{ 
  Ri         <- (sims - data) / data_s
  logLi      <- - log(1+Ri^2) - log(pi*data_s)
  i0         <- which(!is.finite(logLi))
  logLi[i0]  <- 0
  
  sum(logLi*data_w)
}

flogLi <- function(sims,data,data_s,data_w=1)
{ 
  Ri         <- (sims - data) / data_s
  logLi      <- - log(1+Ri^2) - log(pi*data_s)
  i0         <- which(!is.finite(logLi))
  logLi[i0]  <- 0
  
  return(logLi*data_w)
}
