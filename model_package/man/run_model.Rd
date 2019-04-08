\name{run_model}
\alias{run_model}
\title{BASGRA model in Fortran}
\usage{
run_model(p=params, w=matrix_weather, h=days_harvest, n=NDAYS, v=NOUT)
}
\arguments{
  \item{p}{Vector of scaled parameters}
  \item{w}{Matrix of weather}
  \item{h}{Matrix of harvest information}
  \item{n}{Number of days}
  \item{v}{Number of outputs}
}
\description{
Runs Fortran BASGRA via R and C.
}
\examples{
None.
}
