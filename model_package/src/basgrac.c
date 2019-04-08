#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>

void F77_NAME(BASGRA)(double *PARAMS,double *MATRIX_WEATHER,int *DAYS_HARVEST,int NDAYS,int NOUT,double *y);

extern SEXP c_BASGRA(SEXP PARAMS, SEXP MATRIX_WEATHER, SEXP DAYS_HARVEST, SEXP NDAYS, SEXP NOUT){
  const int nd = INTEGER(NDAYS)[0];
  const int no = INTEGER(NOUT)[0];
  SEXP y;
  PROTECT(y = allocMatrix(REALSXP, nd, no));
  F77_CALL(BASGRA)(REAL(PARAMS), REAL(MATRIX_WEATHER), INTEGER(DAYS_HARVEST), nd, no, REAL(y));
  UNPROTECT(1);
  return(y);
}

static const R_CallMethodDef CallEntries[] = {
  {"c_BASGRA",   (DL_FUNC) &c_BASGRA,  5},
  {NULL,          NULL,                0}
};

void R_init_BASGRA (DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
