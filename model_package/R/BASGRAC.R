run_model <- function(p = params,
                      w = matrix_weather,
                      h = days_harvest,
                      n = NDAYS,
                      v = NOUT) {
  if (!is.double(p)) {storage.mode(p) <- 'double'}
  if (!is.double(w)) {storage.mode(w) <- 'double'}
  if (!is.integer(h)) {storage.mode(h) <- 'integer'}
  if (!is.integer(n)) {storage.mode(n) <- 'integer'}
  if (!is.integer(v)) {storage.mode(v) <- 'integer'}
  .Call(c_BASGRA, p, w, h, n, v)
}

# original
# run_model <- function(p = params,
#                       w = matrix_weather,
#                       h = days_harvest,
#                       n = NDAYS) {
#   .Fortran('BASGRA', p,w,h,n, NOUT,matrix(0,n,NOUT))[[6]]
# }

