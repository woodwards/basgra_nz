\name{run_model}
\alias{run_model}
\title{BASGRA model in Fortran}
\usage{
run_model(p=params, w=matrix_weather, h=days_harvest, n=NDAYS, v=NOUT, y=matrix(0,n,v))
}
\arguments{
  \item{p}{Vector of scaled parameters}
  \item{w}{Matrix of weather}
  \item{h}{Matrix of harvest information}
  \item{n}{Number of days}
  \item{v}{Number of outputs}
  \item{y}{Output matrix}
}
\description{
Runs Fortran BASGRA via R and C.
}
\examples{
None.
}
