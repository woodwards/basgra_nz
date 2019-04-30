flogL <- function(sims,data,data_s,data_w=1)
{ 
logLi <- (-1.*0.5)*((sims-data)/data_s)^2 - 0.5*log(2*pi) - log(data_s)

sum(logLi*data_w)
}

flogLi <- function(sims,data,data_s,data_w=1)
{ 
logLi <- (-1.*0.5)*((sims-data)/data_s)^2 - 0.5*log(2*pi) - log(data_s)

return(logLi*data_w)
}
