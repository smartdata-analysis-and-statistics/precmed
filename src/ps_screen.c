#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <limits.h>
#include <float.h>

/* Copyright Smart Data Analysis and Statistics B.V.; Apache-2.0.
 * Screen interaction products without constructing the n-by-choose(p,2)
 * design. Pair ordering matches glm.ps's original i-then-j loops. */
SEXP precmed_ps_screen(SEXP xs, SEXP ts) {
    SEXP dim = Rf_getAttrib(xs, R_DimSymbol);
    if (TYPEOF(xs) != REALSXP || TYPEOF(dim) != INTSXP || LENGTH(dim) != 2 || TYPEOF(ts) != REALSXP)
        Rf_error("PS screening requires a double matrix and treatment vector");
    int n = INTEGER(dim)[0], p = INTEGER(dim)[1];
    if (XLENGTH(ts) != n) Rf_error("PS treatment length does not match the design");
    const double *x = REAL(xs), *t = REAL(ts);
    int count[2] = {0, 0};
    for (int k = 0; k < n; ++k) {
        if (t[k] != 0 && t[k] != 1) Rf_error("PS treatment must contain only 0 and 1");
        ++count[t[k] == 1];
    }
    if (p > 1 && (count[0] < 2 || count[1] < 2))
        Rf_error("Interaction screening requires at least two observations per treatment group");
    R_xlen_t npair = (R_xlen_t)p * (p - 1) / 2, pair = 0;
    unsigned char *keep = (unsigned char *) R_alloc(npair, sizeof(unsigned char));
    int selected = 0;
    for (int i = 0; i < p - 1; ++i) for (int j = i + 1; j < p; ++j, ++pair) {
        if (pair % 1024 == 0) R_CheckUserInterrupt();
        long double sum[2] = {0, 0}, correction[2] = {0, 0}, ss[2] = {0, 0}, delta[2] = {0, 0};
        double mean[2];
        for (int k = 0; k < n; ++k) {
            double z = x[k + (size_t)n*i] * x[k + (size_t)n*j];
            if (!R_FINITE(z)) Rf_error("PS interaction products must be finite");
            sum[t[k] == 1] += z;
        }
        for (int g = 0; g < 2; ++g) mean[g] = (double)(sum[g] / count[g]);
        for (int k = 0; k < n; ++k) {
            double z = x[k + (size_t)n*i] * x[k + (size_t)n*j];
            correction[t[k] == 1] += z - mean[t[k] == 1];
        }
        for (int g = 0; g < 2; ++g) mean[g] += (double)(correction[g] / count[g]);
        for (int k = 0; k < n; ++k) {
            double z = x[k + (size_t)n*i] * x[k + (size_t)n*j];
            int g = t[k] == 1;
            long double d = (long double)z - mean[g];
            ss[g] += d*d; delta[g] += d;
        }
        double v0 = (double)((ss[0] - delta[0]*delta[0]/count[0]) / (count[0]-1));
        double v1 = (double)((ss[1] - delta[1]*delta[1]/count[1]) / (count[1]-1));
        double smd = sqrt(2)*fabs(mean[1]-mean[0]) / sqrt(v1+v0+1e-6);
        /* Recheck borderline decisions with R's exact mean/var routines. */
        double uncertainty = 1e-12 + 16*DBL_EPSILON*(fabs(mean[1])+fabs(mean[0])) /
                             sqrt(v1+v0+1e-6);
        int near = R_FINITE(smd) && fabs(smd - .1) < uncertainty;
        keep[pair] = near ? 2 : (smd > .1 && v1+v0 > 0);
        if (keep[pair]) {
            if (selected == INT_MAX) Rf_error("Too many retained PS interactions");
            ++selected;
        }
    }
    SEXP pairs = PROTECT(Rf_allocMatrix(INTSXP, 2, selected));
    SEXP recheck = PROTECT(Rf_allocVector(LGLSXP, selected));
    pair = 0;
    int out = 0;
    for (int i = 0; i < p - 1; ++i) for (int j = i + 1; j < p; ++j, ++pair) {
        if (!keep[pair]) continue;
        INTEGER(pairs)[2*(R_xlen_t)out] = i + 1;
        INTEGER(pairs)[2*(R_xlen_t)out+1] = j + 1;
        LOGICAL(recheck)[out] = keep[pair] == 2;
        ++out;
    }
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, pairs); SET_VECTOR_ELT(result, 1, recheck);
    UNPROTECT(3);
    return result;
}
