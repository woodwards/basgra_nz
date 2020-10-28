run_model <- function(p = params,
                      w = matrix_weather,
                      h = days_harvest,
                      n = NDAYS,
                      v = NOUT,
                      y = matrix(0,n,v)) {
  if (!is.double(p)) {storage.mode(p) <- 'double'}
  if (!is.double(w)) {storage.mode(w) <- 'double'}
  if (!is.integer(h)) {storage.mode(h) <- 'integer'}
  if (!is.integer(n)) {storage.mode(n) <- 'integer'}
  if (!is.integer(v)) {storage.mode(v) <- 'integer'}
  if (!is.double(y)) {storage.mode(y) <- 'double'}
  if (nrow(y) !=n | ncol(y) != v){
    stop("Error: y matrix is wrong shape in run_model()")
  }
  .Call(c_BASGRA, p, w, h, n, v, y)
}

# http://r-pkgs.had.co.nz/src.html
.onUnload <- function (libpath){
  library.dynam.unload("mypackage", libpath)
}

# original
# run_model <- function(p = params,
#                       w = matrix_weather,
#                       h = days_harvest,
#                       n = NDAYS) {
#   .Fortran('BASGRA', p,w,h,n, NOUT,matrix(0,n,NOUT))[[6]]
# }

