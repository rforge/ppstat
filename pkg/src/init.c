#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void F77_NAME(sgram)(double *sg0, double *sg1, double *sg2,double *sg3, double *tb, int *nb);
		   
SEXP computeModelMatrix(SEXP t, SEXP B, SEXP delta, SEXP s);
		   
static const R_CallMethodDef CallEntries[] = {
    {"computeModelMatrix", (DL_FUNC) &computeModelMatrix, 4},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortEntries[] = {
    {"sgram", (DL_FUNC) &F77_SUB(sgram), 6},
    {NULL, NULL, 0}
};


void  R_init_ppstat(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, FortEntries, NULL);
}
