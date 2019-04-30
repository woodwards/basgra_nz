flogL <- function(sims,data,data_s,data_w=1,dof=7)
{ 
  Ri         <- (sims - data) / data_s
  logLi      <- - ((dof+1)/2)*log(1+(Ri^2)/dof) - 0.5*log(dof*pi*data_s^2) + log(gamma((dof+1)/2)/gamma(dof/2))
  i0         <- which(!is.finite(logLi))
  logLi[i0]  <- 0
  
  sum(logLi*data_w)
}

flogLi <- function(sims,data,data_s,data_w=1,dof=7)
{ 
  Ri         <- (sims - data) / data_s
  logLi      <- - ((dof+1)/2)*log(1+(Ri^2)/dof) - 0.5*log(dof*pi*data_s^2) + log(gamma((dof+1)/2)/gamma(dof/2))
  i0         <- which(!is.finite(logLi))
  logLi[i0]  <- 0
  
  return(logLi*data_w)
}
