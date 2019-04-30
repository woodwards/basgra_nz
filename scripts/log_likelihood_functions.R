# log likelihood

Ri <-  1 # z value
data_s <- 1000 # stadnard error
dof <- 7 # 

# warning need non-standardised versions!

logLi <- (-1.*0.5)*(Ri)^2 - 0.5*log(2*pi*data_s^2) # Gauss
logLi <- log(1-exp(-0.5*Ri^2)) - 2*log(abs(Ri)) - 0.5*log(2*pi*data_s^2) # Sivia
logLi <- - log(1+Ri^2) - log(pi*data_s) # Cauchy
logLi <- - ((dof+1)/2)*log(1+(Ri^2)/dof) - 0.5*log(dof*pi*data_s^2) + log(gamma((dof+1)/2)/gamma(dof/2)) # Students
