#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <float.h>
#include <math.h>
#include <string.h>

/* Copyright Smart Data Analysis and Statistics B.V.; Apache-2.0. */
typedef struct {
    int n, p;
    const double *x, *y, *trt, *ps, *f1, *f0, *weights;
    double *eta, *error, *weighted, *score, *slope;
} Contrast;

static Contrast prepare(SEXP data, int survival) {
    if (TYPEOF(data) != VECSXP || XLENGTH(data) != (survival ? 7 : 6))
        Rf_error("Invalid contrast-regression data");
    SEXP x = VECTOR_ELT(data, 0), dim = Rf_getAttrib(x, R_DimSymbol);
    if (TYPEOF(x) != REALSXP || TYPEOF(dim) != INTSXP || LENGTH(dim) != 2)
        Rf_error("The contrast design must be a double matrix");
    Contrast c;
    c.n = INTEGER(dim)[0]; c.p = INTEGER(dim)[1];
    if (c.n < 1 || c.p < 1) Rf_error("Empty contrast design");
    c.x = REAL(x);
    for (int j = 1; j < (survival ? 7 : 6); ++j)
        if (TYPEOF(VECTOR_ELT(data, j)) != REALSXP || XLENGTH(VECTOR_ELT(data, j)) != c.n)
            Rf_error("Invalid contrast vector length or type");
    c.y = REAL(VECTOR_ELT(data, 1)); c.trt = REAL(VECTOR_ELT(data, 2));
    c.ps = REAL(VECTOR_ELT(data, 3)); c.f1 = REAL(VECTOR_ELT(data, 4));
    c.f0 = REAL(VECTOR_ELT(data, 5));
    c.weights = survival ? REAL(VECTOR_ELT(data, 6)) : NULL;
    c.eta = (double *) R_alloc(c.n, sizeof(double));
    c.error = (double *) R_alloc(c.n, sizeof(double));
    c.weighted = (double *) R_alloc((size_t)c.n * c.p, sizeof(double));
    c.score = (double *) R_alloc(c.p, sizeof(double));
    c.slope = (double *) R_alloc((size_t)c.p * c.p, sizeof(double));
    return c;
}

static double evaluate(Contrast *c, const double *beta, int slope, int stop_on_overflow) {
    double eta_max = 0;
    for (int i = 0; i < c->n; ++i) {
        double linear = 0;
        for (int j = 0; j < c->p; ++j) linear += c->x[i + (size_t)c->n*j] * beta[j];
        c->eta[i] = exp(linear);
        if (c->eta[i] > eta_max) eta_max = c->eta[i];
    }
    /* Survival NR stops before calculating a score when eta overflows. */
    if (stop_on_overflow && eta_max == R_PosInf) return eta_max;
    for (int i = 0; i < c->n; ++i) {
        double eta = c->eta[i], ps = c->ps[i], trt = c->trt[i];
        double den = eta * ps + (1 - ps);
        c->error[i] = (trt * (c->y[i] - eta*c->f0[i]/2 - c->f1[i]/2) * (1-ps)
                     - (1-trt) * (eta*c->y[i] - eta*c->f0[i]/2 - c->f1[i]/2) * ps) / den;
        if (slope) {
            double w = (c->y[i] + c->f0[i]*(trt/ps-1)/2
                       + c->f1[i]*((1-trt)/(1-ps)-1)/2) * eta*ps*(1-ps)/(den*den);
            for (int j = 0; j < c->p; ++j) {
                double value = c->x[i + (size_t)c->n*j];
                if (c->weights) value *= c->weights[i];
                c->weighted[i + (size_t)c->n*j] = value * w;
            }
        }
    }
    for (int j = 0; j < c->p; ++j) {
        long double total = 0;
        for (int i = 0; i < c->n; ++i) {
            double value = c->x[i + (size_t)c->n*j];
            if (c->weights) value *= c->weights[i];
            total += value * c->error[i];
        }
        c->score[j] = (double) total;
    }
    if (slope) {
        double one = 1, zero = 0;
        F77_CALL(dgemm)("T", "N", &c->p, &c->p, &c->n, &one,
                        c->weighted, &c->n, c->x, &c->n, &zero,
                        c->slope, &c->p FCONE FCONE);
    }
    return eta_max;
}

static SEXP contrast_loss(SEXP data, SEXP beta, int survival) {
    Contrast c = prepare(data, survival);
    if (TYPEOF(beta) != REALSXP || XLENGTH(beta) != c.p)
        Rf_error("Invalid contrast coefficients");
    evaluate(&c, REAL(beta), 0, 0);
    long double loss = 0;
    for (int j = 0; j < c.p; ++j) loss += c.score[j] * c.score[j];
    return Rf_ScalarReal((double) loss);
}

static SEXP contrast_nr(SEXP data, SEXP tolerance, SEXP maxiter, SEXP tuning,
                        int survival, double tau0) {
    Contrast c = prepare(data, survival);
    double tol = Rf_asReal(tolerance);
    int limit = Rf_asInteger(maxiter);
    if (!R_FINITE(tol) || tol <= 0 || limit < 0 || TYPEOF(tuning) != REALSXP || LENGTH(tuning) != 2)
        Rf_error("Invalid contrast solver controls");
    const double *tune = REAL(tuning);
    if (!R_FINITE(tune[0]) || !R_FINITE(tune[1]) || tune[0] <= 0 || tune[1] <= 0)
        Rf_error("Invalid contrast tuning values");
    SEXP beta = PROTECT(Rf_allocVector(REALSXP, c.p));
    memset(REAL(beta), 0, c.p*sizeof(double));
    size_t pp = (size_t)c.p*c.p;
    double *a = (double *) R_alloc(pp, sizeof(double));
    double *inverse = (double *) R_alloc(pp, sizeof(double));
    double *eigenvalues = (double *) R_alloc(c.p, sizeof(double));
    int lwork = 3*c.p > 1 ? 3*c.p : 1, info;
    double *work = (double *) R_alloc(4*(size_t)c.p + lwork, sizeof(double));
    int *pivots = (int *) R_alloc(c.p, sizeof(int));
    int *iwork = (int *) R_alloc(c.p, sizeof(int));
    double mae = R_PosInf, emin = R_PosInf, epsilon0 = 0, eta_max = 0;
    int iter = 0;
    while (mae > tol && iter <= limit &&
           (survival ? eta_max < 1000 * tau0 : emin > 0)) {
        R_CheckUserInterrupt();
        eta_max = evaluate(&c, REAL(beta), 1, survival);
        if (survival && eta_max == R_PosInf) break;
        for (size_t j = 0; j < pp; ++j)
            if (!R_FINITE(c.slope[j])) Rf_error("infinite or missing values in contrast slope");
        /* Count checks positivity every iteration; survival computes the
         * regularization once and uses its historical eta cap instead. */
        if (!survival || iter == 0) {
            memcpy(a, c.slope, pp*sizeof(double));
            F77_CALL(dsyev)("N", "L", &c.p, a, &c.p, eigenvalues, work, &lwork, &info FCONE FCONE);
            if (info != 0) Rf_error("Contrast slope eigendecomposition failed");
            emin = eigenvalues[0];
            if (iter == 0) epsilon0 = emin + emin*(emin < 0);
        }
        memcpy(a, c.slope, pp*sizeof(double));
        memset(inverse, 0, pp*sizeof(double));
        double norm = 0;
        for (int j = 0; j < c.p; ++j) {
            a[j + (size_t)c.p*j] += tune[1]*fabs(epsilon0);
            inverse[j + (size_t)c.p*j] = 1;
            double total = 0;
            for (int i = 0; i < c.p; ++i) total += fabs(a[i + (size_t)c.p*j]);
            if (total > norm) norm = total;
        }
        /* Preserve the R path's inverse-times-score update and singularity guard. */
        F77_CALL(dgesv)(&c.p, &c.p, a, &c.p, pivots, inverse, &c.p, &info);
        if (info != 0) Rf_error("Contrast slope is singular");
        double rcond;
        F77_CALL(dgecon)("1", &c.p, a, &c.p, &norm, &rcond, work, iwork, &info FCONE);
        if (info != 0 || rcond < DBL_EPSILON) Rf_error("Contrast slope is computationally singular");
        for (int j = 0; j < c.p; ++j) {
            double step = 0;
            for (int k = 0; k < c.p; ++k) step += inverse[j + (size_t)c.p*k] * c.score[k];
            REAL(beta)[j] += step*tune[0];
        }
        long double total = 0;
        for (int j = 0; j < c.p; ++j) total += fabs(c.score[j]);
        mae = (double) total;
        if (!R_FINITE(mae)) Rf_error("non-finite contrast score");
        ++iter;
    }
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(out, 0, beta);
    SET_VECTOR_ELT(out, 1, Rf_ScalarLogical(mae <= tol));
    UNPROTECT(2);
    return out;
}

/* Registered entry points keep the count ABI unchanged. */
SEXP precmed_contrast_nr(SEXP data, SEXP tolerance, SEXP maxiter, SEXP tuning) {
    return contrast_nr(data, tolerance, maxiter, tuning, 0, 0);
}
SEXP precmed_contrast_loss(SEXP data, SEXP beta) {
    return contrast_loss(data, beta, 0);
}
SEXP precmed_survival_contrast_nr(SEXP data, SEXP tolerance, SEXP maxiter,
                                 SEXP tuning, SEXP tau) {
    double tau0 = Rf_asReal(tau);
    if (!R_FINITE(tau0) || tau0 <= 0) Rf_error("tau0 must be finite and positive");
    return contrast_nr(data, tolerance, maxiter, tuning, 1, tau0);
}
SEXP precmed_survival_contrast_loss(SEXP data, SEXP beta) {
    return contrast_loss(data, beta, 1);
}
