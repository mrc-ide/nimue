#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void vaccine_initmod_desolve(void *);
extern void vaccine_output_dde(void *);
extern void vaccine_rhs_dde(void *);
extern void vaccine_rhs_desolve(void *);

/* .Call calls */
extern SEXP vaccine_contents(SEXP);
extern SEXP vaccine_create(SEXP);
extern SEXP vaccine_initial_conditions(SEXP, SEXP);
extern SEXP vaccine_metadata(SEXP);
extern SEXP vaccine_rhs_r(SEXP, SEXP, SEXP);
extern SEXP vaccine_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP vaccine_set_user(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"vaccine_initmod_desolve", (DL_FUNC) &vaccine_initmod_desolve, 1},
    {"vaccine_output_dde",      (DL_FUNC) &vaccine_output_dde,      1},
    {"vaccine_rhs_dde",         (DL_FUNC) &vaccine_rhs_dde,         1},
    {"vaccine_rhs_desolve",     (DL_FUNC) &vaccine_rhs_desolve,     1},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"vaccine_contents",           (DL_FUNC) &vaccine_contents,           1},
    {"vaccine_create",             (DL_FUNC) &vaccine_create,             1},
    {"vaccine_initial_conditions", (DL_FUNC) &vaccine_initial_conditions, 2},
    {"vaccine_metadata",           (DL_FUNC) &vaccine_metadata,           1},
    {"vaccine_rhs_r",              (DL_FUNC) &vaccine_rhs_r,              3},
    {"vaccine_set_initial",        (DL_FUNC) &vaccine_set_initial,        4},
    {"vaccine_set_user",           (DL_FUNC) &vaccine_set_user,           2},
    {NULL, NULL, 0}
};

void R_init_nimue(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
