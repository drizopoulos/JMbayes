#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP JMbayes_dmvnorm2(SEXP, SEXP, SEXP, SEXP);
extern SEXP JMbayes_gradient_logPosterior(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP JMbayes_gradient_logPosterior_nogammas(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP JMbayes_lap_rwm_C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP JMbayes_lap_rwm_C_woRE(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP JMbayes_lap_rwm_C_woRE_nogammas(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP JMbayes_logPosterior(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP JMbayes_logPosterior_nogammas(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"JMbayes_dmvnorm2",                       (DL_FUNC) &JMbayes_dmvnorm2,                        4},
    {"JMbayes_gradient_logPosterior",          (DL_FUNC) &JMbayes_gradient_logPosterior,          22},
    {"JMbayes_gradient_logPosterior_nogammas", (DL_FUNC) &JMbayes_gradient_logPosterior_nogammas, 17},
    {"JMbayes_lap_rwm_C",                      (DL_FUNC) &JMbayes_lap_rwm_C,                       7},
    {"JMbayes_lap_rwm_C_woRE",                 (DL_FUNC) &JMbayes_lap_rwm_C_woRE,                  6},
    {"JMbayes_lap_rwm_C_woRE_nogammas",        (DL_FUNC) &JMbayes_lap_rwm_C_woRE_nogammas,         6},
    {"JMbayes_logPosterior",                   (DL_FUNC) &JMbayes_logPosterior,                   22},
    {"JMbayes_logPosterior_nogammas",          (DL_FUNC) &JMbayes_logPosterior_nogammas,          17},
    {NULL, NULL, 0}
};

void R_init_JMbayes(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}