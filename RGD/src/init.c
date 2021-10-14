#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _RGD_cplexcoef(SEXP, SEXP, SEXP);
extern SEXP _RGD_GD(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RGD_picksamples(SEXP, SEXP, SEXP);
extern SEXP _RGD_RGD(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RGD_runif_in_pball(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_RGD_cplexcoef",      (DL_FUNC) &_RGD_cplexcoef,      3},
    {"_RGD_GD",             (DL_FUNC) &_RGD_GD,             6},
    {"_RGD_picksamples",    (DL_FUNC) &_RGD_picksamples,    3},
    {"_RGD_RGD",            (DL_FUNC) &_RGD_RGD,            6},
    {"_RGD_runif_in_pball", (DL_FUNC) &_RGD_runif_in_pball, 4},
    {NULL, NULL, 0}
};

void R_init_RGD(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
