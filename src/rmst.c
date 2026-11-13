#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <Rmath.h>

/* Copyright Smart Data Analysis and Statistics B.V.; Apache-2.0.
 * Integrate one Cox survival curve at a time. Auxiliary memory is O(1),
 * excluding the input event grid and the returned patient estimates. */
SEXP precmed_rmst(SEXP survival, SEXP hr, SEXP gap) {
    if (TYPEOF(survival) != REALSXP || TYPEOF(hr) != REALSXP || TYPEOF(gap) != REALSXP)
        Rf_error("RMST inputs must be double vectors");
    R_xlen_t k = XLENGTH(survival), n = XLENGTH(hr);
    if (XLENGTH(gap) != k + 1) Rf_error("RMST interval lengths do not match event times");
    const double *s = REAL(survival), *h = REAL(hr), *g = REAL(gap);
    for (R_xlen_t j = 0; j < k; ++j)
        if (!R_FINITE(s[j]) || s[j] < 0 || s[j] > 1)
            Rf_error("RMST baseline survival must be finite and in [0, 1]");
    for (R_xlen_t j = 0; j <= k; ++j)
        if (!R_FINITE(g[j]) || g[j] < 0)
            Rf_error("RMST interval lengths must be finite and non-negative");
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
    for (R_xlen_t i = 0; i < n; ++i) {
        if (i % 1024 == 0) R_CheckUserInterrupt();
        if (ISNAN(h[i]) || h[i] < 0) Rf_error("RMST hazard ratios must be non-negative and not missing");
        long double total = g[0];
        for (R_xlen_t j = 0; j < k; ++j) {
            /* Round products to double before accumulating, as in R's
             * elementwise multiplication followed by colSums(). */
            double product = R_pow(s[j], h[i]) * g[j + 1];
            total += product;
        }
        REAL(out)[i] = (double) total;
    }
    UNPROTECT(1);
    return out;
}
