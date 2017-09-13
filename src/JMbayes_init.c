#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _JMbayes_dmvnorm2(SEXP, SEXP, SEXP, SEXP);
extern SEXP _JMbayes_gradient_logPosterior(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _JMbayes_gradient_logPosterior_nogammas(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _JMbayes_lap_rwm_C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _JMbayes_lap_rwm_C_woRE(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _JMbayes_lap_rwm_C_woRE_nogammas(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _JMbayes_log_post_RE_svft(SEXP, SEXP);
extern SEXP _JMbayes_logPosterior(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _JMbayes_logPosterior_nogammas(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _JMbayes_survPred_svft_2(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_JMbayes_dmvnorm2",                       (DL_FUNC) &_JMbayes_dmvnorm2,                        4},
    {"_JMbayes_gradient_logPosterior",          (DL_FUNC) &_JMbayes_gradient_logPosterior,          22},
    {"_JMbayes_gradient_logPosterior_nogammas", (DL_FUNC) &_JMbayes_gradient_logPosterior_nogammas, 17},
    {"_JMbayes_lap_rwm_C",                      (DL_FUNC) &_JMbayes_lap_rwm_C,                       7},
    {"_JMbayes_lap_rwm_C_woRE",                 (DL_FUNC) &_JMbayes_lap_rwm_C_woRE,                  6},
    {"_JMbayes_lap_rwm_C_woRE_nogammas",        (DL_FUNC) &_JMbayes_lap_rwm_C_woRE_nogammas,         6},
    {"_JMbayes_log_post_RE_svft",               (DL_FUNC) &_JMbayes_log_post_RE_svft,                2},
    {"_JMbayes_logPosterior",                   (DL_FUNC) &_JMbayes_logPosterior,                   22},
    {"_JMbayes_logPosterior_nogammas",          (DL_FUNC) &_JMbayes_logPosterior_nogammas,          17},
    {"_JMbayes_survPred_svft_2",                (DL_FUNC) &_JMbayes_survPred_svft_2,                 2},
    {NULL, NULL, 0}
};

void R_init_JMbayes(DllInfo *dll) {
   R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
   R_useDynamicSymbols(dll, FALSE);
}



