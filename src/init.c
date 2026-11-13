#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP precmed_contrast_nr(SEXP, SEXP, SEXP, SEXP);
extern SEXP precmed_contrast_loss(SEXP, SEXP);
extern SEXP precmed_rmst(SEXP, SEXP, SEXP);
extern SEXP precmed_survival_contrast_nr(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP precmed_survival_contrast_loss(SEXP, SEXP);
extern SEXP precmed_ps_screen(SEXP, SEXP);
static const R_CallMethodDef methods[] = {
    {"precmed_contrast_nr", (DL_FUNC) &precmed_contrast_nr, 4},
    {"precmed_contrast_loss", (DL_FUNC) &precmed_contrast_loss, 2},
    {"precmed_rmst", (DL_FUNC) &precmed_rmst, 3},
    {"precmed_survival_contrast_nr", (DL_FUNC) &precmed_survival_contrast_nr, 5},
    {"precmed_survival_contrast_loss", (DL_FUNC) &precmed_survival_contrast_loss, 2},
    {"precmed_ps_screen", (DL_FUNC) &precmed_ps_screen, 2},
    {NULL, NULL, 0}
};
void attribute_visible R_init_precmed(DllInfo *dll) {
    R_registerRoutines(dll, NULL, methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
