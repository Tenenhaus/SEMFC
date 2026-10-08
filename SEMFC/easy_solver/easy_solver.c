#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifndef FCONE
#define FCONE
#endif

/* Ordre des parametres dans le vecteur ctrl (identique cote R) */
enum
{
    P_MAXIT,
    P_TOL,
    P_KKT_TOL,
    P_STEPSIZE,
    P_STEPSIZE_MIN,
    P_STEPSIZE_MAX,
    P_MAX_STEP,
    P_BACKTRACK,
    P_F_MIN,
    P_EQ_TOL,
    P_ACTIVE_TOL,
    P_MAX_EQ,
    P_TRACE,
    P_USE_BFGS,
    P_CURV_TOL,
    P_BFGS_SCALE,
    P_GRAD_EPS,
    P_COUNT
};

typedef struct
{
    SEXP f, h, rho;
    int n, m;
    const double *lb, *eqB;
    double tol, kkt_tol, stepsize_min, stepsize_max, max_step, backtrack,
        f_min, eq_tol, active_tol, curv_tol, grad_eps;
    int max_eq;
} Ctx;

#define ALLOC(k) ((double *)R_alloc((size_t)((k) > 0 ? (k) : 1), sizeof(double)))
#define IALLOC(k) ((int *)R_alloc((size_t)((k) > 0 ? (k) : 1), sizeof(int)))

/* ---------------------------------------------------------------- */
/* Utilitaires                                                      */
/* ---------------------------------------------------------------- */

static double maxabs(const double *v, int len)
{
    double mx = 0.0;
    for (int i = 0; i < len; i++)
    {
        double a = fabs(v[i]);
        if (ISNAN(a))
            return R_NaN;
        if (a > mx)
            mx = a;
    }
    return mx;
}

static double maxv(double a, double b)
{
    if (ISNAN(a) || ISNAN(b))
        return R_NaN;
    return a > b ? a : b;
}

/* y = M x, M est nr x nc (colonnes) */
static void matvec(const double *M, int nr, int nc, const double *x, double *y)
{
    for (int i = 0; i < nr; i++)
        y[i] = 0.0;
    for (int j = 0; j < nc; j++)
    {
        double xj = x[j];
        for (int i = 0; i < nr; i++)
            y[i] += M[i + (size_t)j * nr] * xj;
    }
}

/* y = M' x */
static void matTvec(const double *M, int nr, int nc, const double *x, double *y)
{
    for (int j = 0; j < nc; j++)
    {
        double s = 0.0;
        for (int i = 0; i < nr; i++)
            s += M[i + (size_t)j * nr] * x[i];
        y[j] = s;
    }
}

static void gemm(const char *ta, const char *tb, int M, int N, int K,
                 const double *A, int lda, const double *B, int ldb,
                 double *Cm, int ldc)
{
    double one = 1.0, zero = 0.0;
    F77_CALL(dgemm)(ta, tb, &M, &N, &K, &one, A, &lda, B, &ldb, &zero, Cm, &ldc FCONE FCONE);
}

/* Moindres carres de norme minimale par SVD tronquee (seuil 1e-10 * dmax)
   A : m x k (colonnes), b : m, x : k  -> x = pinv(A) b                   */
static void lstsq(const double *A, int m, int k, const double *b, double *x)
{
    for (int j = 0; j < k; j++)
        x[j] = 0.0;
    if (m <= 0 || k <= 0)
        return;

    for (size_t i = 0; i < (size_t)m * k; i++)
        if (!R_FINITE(A[i]))
            error("non-finite values in matrix passed to least squares");
    for (int i = 0; i < m; i++)
        if (!R_FINITE(b[i]))
            error("non-finite values in vector passed to least squares");

    int r = m < k ? m : k, info = 0, lwork = -1;
    double wq = 0.0;
    double *Ac = (double *)malloc(sizeof(double) * (size_t)m * k);
    double *S = (double *)malloc(sizeof(double) * (size_t)r);
    double *U = (double *)malloc(sizeof(double) * (size_t)m * r);
    double *VT = (double *)malloc(sizeof(double) * (size_t)r * k);
    int *iw = (int *)malloc(sizeof(int) * 8 * (size_t)r);
    double *work = NULL;
    int failed = 0;

    if (!Ac || !S || !U || !VT || !iw)
    {
        failed = 1;
    }
    else
    {
        memcpy(Ac, A, sizeof(double) * (size_t)m * k);
        F77_CALL(dgesdd)("S", &m, &k, Ac, &m, S, U, &m, VT, &r, &wq, &lwork, iw, &info FCONE);
        if (info != 0)
            failed = 1;
        else
        {
            lwork = (int)wq;
            work = (double *)malloc(sizeof(double) * (size_t)lwork);
            if (!work)
                failed = 1;
            else
            {
                F77_CALL(dgesdd)("S", &m, &k, Ac, &m, S, U, &m, VT, &r, work, &lwork, iw, &info FCONE);
                if (info != 0)
                    failed = 1;
                else if (S[0] > 0.0)
                {
                    double thr = S[0] * 1e-10;
                    for (int i = 0; i < r; i++)
                    {
                        if (S[i] > thr)
                        {
                            double dot = 0.0;
                            for (int p = 0; p < m; p++)
                                dot += U[p + (size_t)i * m] * b[p];
                            double cf = dot / S[i];
                            for (int j = 0; j < k; j++)
                                x[j] += cf * VT[i + (size_t)j * r];
                        }
                    }
                }
            }
        }
    }
    free(Ac);
    free(S);
    free(U);
    free(VT);
    free(iw);
    free(work);
    if (failed)
        error("SVD failed in least squares");
}

/* ---------------------------------------------------------------- */
/* Appels des closures R                                            */
/* ---------------------------------------------------------------- */

/* Retourne le resultat PROTEGE (l'appelant fait UNPROTECT(1)) */
static SEXP call1(SEXP fn, SEXP rho, const double *x, int n)
{
    SEXP xv = PROTECT(allocVector(REALSXP, n));
    memcpy(REAL(xv), x, sizeof(double) * (size_t)n);
    SEXP cl = PROTECT(lang2(fn, xv));
    SEXP res = eval(cl, rho);
    UNPROTECT(2);
    PROTECT(res);
    return res;
}

static void copy_real(SEXP r, double *out, int len, const char *what)
{
    if (xlength(r) != (R_xlen_t)len)
        error("%s has wrong length", what);
    SEXP rd = PROTECT(coerceVector(r, REALSXP));
    if (len > 0)
        memcpy(out, REAL(rd), sizeof(double) * (size_t)len);
    UNPROTECT(1);
}

static double eval_f(Ctx *c, const double *x)
{
    SEXP r = call1(c->f, c->rho, x, c->n);
    if (xlength(r) != 1)
        error("fun must return a single number");
    double v = asReal(r);
    UNPROTECT(1);
    return v;
}

static void eval_h(Ctx *c, const double *x, double *out)
{
    if (c->m == 0)
        return;
    SEXP r = call1(c->h, c->rho, x, c->n);
    copy_real(r, out, c->m, "eqfun result");
    UNPROTECT(1);
}

/* gradient (et jacobienne si m > 0) en un appel */
/* Differences finies centrees : g (si g != NULL) et J (si m > 0) */
static void fd_derivs(Ctx *c, const double *x, double *g, double *J)
{
    int n = c->n, m = c->m;
    const void *vmax = vmaxget();
    double *xp = ALLOC(n), *hp = ALLOC(m), *hm = ALLOC(m);
    memcpy(xp, x, sizeof(double) * (size_t)n);

    for (int j = 0; j < n; j++)
    {
        R_CheckUserInterrupt();
        double ax = fabs(x[j]);
        double dx = c->grad_eps * (ax > 1.0 ? ax : 1.0);
        double fp = 0.0, fm = 0.0;

        xp[j] = x[j] + dx;
        if (g)
            fp = eval_f(c, xp);
        if (m > 0)
            eval_h(c, xp, hp);

        xp[j] = x[j] - dx;
        if (g)
            fm = eval_f(c, xp);
        if (m > 0)
            eval_h(c, xp, hm);

        xp[j] = x[j];

        if (g)
        {
            if (!R_FINITE(fp) || !R_FINITE(fm))
                error("Non finite objective during gradient calculation");
            g[j] = (fp - fm) / (2.0 * dx);
        }
        for (int i = 0; i < m; i++)
            J[i + (size_t)j * m] = (hp[i] - hm[i]) / (2.0 * dx);
    }
    vmaxset(vmax);
}

static void eval_gj(Ctx *c, const double *x, double *g, double *J)
{
    fd_derivs(c, x, g, J);
}

static void eval_jac(Ctx *c, const double *x, double *J)
{
    fd_derivs(c, x, NULL, J);
}

/* ---------------------------------------------------------------- */
/* KKT                                                              */
/* ---------------------------------------------------------------- */

typedef struct
{
    double kkt, eq_res;
    int nact;
} Kkt;

static void kkt_measure(Ctx *c, const double *x, const double *g, const double *J,
                        const double *hx, double *r, double *lambda, Kkt *out)
{
    int n = c->n, m = c->m;
    const void *vmax = vmaxget();
    int *act = IALLOC(n), *fr = IALLOC(n);
    int na = 0, nf = 0;

    for (int j = 0; j < n; j++)
    {
        if (R_FINITE(c->lb[j]) && x[j] <= c->lb[j] + c->active_tol)
            act[na++] = j;
        else
            fr[nf++] = j;
    }

    for (int j = 0; j < n; j++)
        r[j] = g[j];
    for (int i = 0; i < m; i++)
        lambda[i] = 0.0;

    if (m > 0 && nf > 0)
    {
        double *A = ALLOC((size_t)nf * m), *b = ALLOC(nf);
        for (int a = 0; a < nf; a++)
        {
            b[a] = -g[fr[a]];
            for (int i = 0; i < m; i++)
                A[a + (size_t)i * nf] = J[i + (size_t)fr[a] * m];
        }
        lstsq(A, nf, m, b, lambda);
        for (int j = 0; j < n; j++)
        {
            double s = 0.0;
            for (int i = 0; i < m; i++)
                s += J[i + (size_t)j * m] * lambda[i];
            r[j] += s;
        }
    }

    double free_res = 0.0, act_viol = 0.0, eq_res = 0.0;
    for (int a = 0; a < nf; a++)
        free_res = maxv(free_res, fabs(r[fr[a]]));
    for (int a = 0; a < na; a++)
        act_viol = maxv(act_viol, (-r[act[a]] > 0.0) ? -r[act[a]] : 0.0);
    if (m > 0)
    {
        double *t = ALLOC(m);
        for (int i = 0; i < m; i++)
            t[i] = hx[i] - c->eqB[i];
        eq_res = maxabs(t, m);
    }

    out->kkt = maxv(maxv(free_res, act_viol), eq_res);
    out->eq_res = eq_res;
    out->nact = na;
    vmaxset(vmax);
}

/* ---------------------------------------------------------------- */
/* Direction admissible (equalites + bornes inf, active set)        */
/* H == NULL : projection euclidienne ; sinon metrique H (BFGS)     */
/* ---------------------------------------------------------------- */

static void feasible_direction(Ctx *c, const double *x, const double *g,
                               const double *J, const double *H, double *d)
{
    int n = c->n, m = c->m;
    const void *vmax = vmaxget();
    int *active = IALLOC(n), *fixed = IALLOC(n);
    int na = 0, nfx = 0;

    for (int j = 0; j < n; j++)
        if (R_FINITE(c->lb[j]) && x[j] <= c->lb[j] + c->active_tol)
            active[na++] = j;

    int kmax = m + n;
    double *A = ALLOC((size_t)kmax * n), *v = ALLOC(n), *w = ALLOC(n);
    double *t = ALLOC(n), *t2 = ALLOC(n);
    double *Av = ALLOC(kmax), *lam = ALLOC(kmax);
    double *AAt = ALLOC((size_t)kmax * kmax), *AH = ALLOC((size_t)kmax * n);

    for (int it = 0; it <= n; it++)
    {
        int k = m + nfx;

        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < m; i++)
                A[i + (size_t)j * k] = J[i + (size_t)j * m];
            for (int q = 0; q < nfx; q++)
                A[(m + q) + (size_t)j * k] = 0.0;
        }
        for (int q = 0; q < nfx; q++)
            A[(m + q) + (size_t)fixed[q] * k] = 1.0;

        for (int j = 0; j < n; j++)
            v[j] = -g[j];

        if (H)
        {
            matvec(H, n, n, v, w);
            if (k == 0)
            {
                memcpy(d, w, sizeof(double) * (size_t)n);
            }
            else
            {
                gemm("N", "N", k, n, n, A, k, H, n, AH, k);
                gemm("N", "T", k, k, n, AH, k, A, k, AAt, k);
                matvec(AH, k, n, v, Av);
                lstsq(AAt, k, k, Av, lam);
                matTvec(A, k, n, lam, t);
                matvec(H, n, n, t, t2);
                for (int j = 0; j < n; j++)
                    d[j] = w[j] - t2[j];
            }
        }
        else
        {
            if (k == 0)
            {
                memcpy(d, v, sizeof(double) * (size_t)n);
            }
            else
            {
                gemm("N", "T", k, k, n, A, k, A, k, AAt, k);
                matvec(A, k, n, v, Av);
                lstsq(AAt, k, k, Av, lam);
                matTvec(A, k, n, lam, t);
                for (int j = 0; j < n; j++)
                    d[j] = v[j] - t[j];
            }
        }

        if (na == 0)
            break;

        int newbad = 0;
        for (int a = 0; a < na; a++)
        {
            int j = active[a];
            if (d[j] < -1e-14)
            {
                int already = 0;
                for (int q = 0; q < nfx; q++)
                    if (fixed[q] == j)
                    {
                        already = 1;
                        break;
                    }
                if (!already)
                {
                    fixed[nfx++] = j;
                    newbad++;
                }
            }
        }
        if (newbad == 0)
            break;
    }

    double nd = 0.0;
    for (int j = 0; j < n; j++)
        nd += d[j] * d[j];
    nd = sqrt(nd);

    if (!R_FINITE(nd) || nd == 0.0)
    {
        for (int j = 0; j < n; j++)
            d[j] = 0.0;
    }
    else if (nd > c->max_step)
    {
        double sc = c->max_step / nd;
        for (int j = 0; j < n; j++)
            d[j] *= sc;
    }
    vmaxset(vmax);
}

/* ---------------------------------------------------------------- */
/* Correction d'egalite (Gauss-Newton + Broyden, jacobienne exacte  */
/* recalculee au point de depart)                                   */
/* ---------------------------------------------------------------- */

static int eq_correction(Ctx *c, const double *x_in, double *x_out)
{
    int n = c->n, m = c->m;
    memcpy(x_out, x_in, sizeof(double) * (size_t)n);
    if (m == 0)
        return 1;

    const void *vmax = vmaxget();
    double *J = ALLOC((size_t)m * n), *hcur = ALLOC(m), *ht = ALLOC(m);
    double *r = ALLOC(m), *mr = ALLOC(m), *tmp = ALLOC(m), *Js = ALLOC(m);
    double *dx = ALLOC(n), *xt = ALLOC(n);
    int have_J = 0, J_fresh = 0, ok = 0;

    eval_h(c, x_out, hcur);
    for (int i = 0; i < m; i++)
        r[i] = hcur[i] - c->eqB[i];

    for (int it = 0; it < c->max_eq; it++)
    {
        double r_norm = maxabs(r, m);
        if (r_norm <= c->eq_tol)
        {
            ok = 1;
            goto done;
        }

        if (!have_J)
        {
            eval_jac(c, x_out, J);
            have_J = 1;
            J_fresh = 1;
        }

        for (int i = 0; i < m; i++)
            mr[i] = -r[i];
        lstsq(J, m, n, mr, dx);

        double ndx = 0.0;
        for (int j = 0; j < n; j++)
            ndx += dx[j] * dx[j];
        ndx = sqrt(ndx);
        if (!R_FINITE(ndx))
        {
            ok = 0;
            goto done;
        }
        if (ndx > c->max_step)
        {
            double sc = c->max_step / ndx;
            for (int j = 0; j < n; j++)
                dx[j] *= sc;
        }

        int accepted = 0;
        double alpha = 1.0;
        while (alpha >= 1e-6)
        {
            for (int j = 0; j < n; j++)
            {
                xt[j] = x_out[j] + alpha * dx[j];
                if (xt[j] < c->lb[j])
                    xt[j] = c->lb[j];
            }
            eval_h(c, xt, ht);
            for (int i = 0; i < m; i++)
                tmp[i] = ht[i] - c->eqB[i];
            if (maxabs(tmp, m) < r_norm)
            {
                accepted = 1;
                break;
            }
            alpha *= c->backtrack;
        }

        if (accepted)
        {
            double ss = 0.0;
            for (int j = 0; j < n; j++)
            {
                dx[j] = xt[j] - x_out[j];
                ss += dx[j] * dx[j];
            }
            for (int i = 0; i < m; i++)
                tmp[i] = ht[i] - hcur[i];
            matvec(J, m, n, dx, Js);
            if (ss > 0.0)
            {
                for (int j = 0; j < n; j++)
                    for (int i = 0; i < m; i++)
                        J[i + (size_t)j * m] += (tmp[i] - Js[i]) * dx[j] / ss;
            }
            memcpy(x_out, xt, sizeof(double) * (size_t)n);
            for (int i = 0; i < m; i++)
            {
                hcur[i] = ht[i];
                r[i] = ht[i] - c->eqB[i];
            }
            J_fresh = 0;
        }
        else if (J_fresh)
        {
            ok = 0;
            goto done;
        }
        else
        {
            have_J = 0;
        }
    }
    ok = (maxabs(r, m) <= c->eq_tol);

done:
    vmaxset(vmax);
    return ok;
}

/* ---------------------------------------------------------------- */
/* Historique et sorties                                            */
/* ---------------------------------------------------------------- */

static void add_hist(double *hist, int maxit, int *nh, int iter, double f,
                     double kkt, double eq, double gn, double step, int acc)
{
    int r = *nh;
    hist[r] = iter;
    hist[maxit + r] = f;
    hist[2 * maxit + r] = kkt;
    hist[3 * maxit + r] = eq;
    hist[4 * maxit + r] = gn;
    hist[5 * maxit + r] = step;
    hist[6 * maxit + r] = acc;
    (*nh)++;
}

static SEXP mkvec(const double *p, int len)
{
    SEXP v = allocVector(REALSXP, len);
    if (len > 0)
        memcpy(REAL(v), p, sizeof(double) * (size_t)len);
    return v;
}

/* ---------------------------------------------------------------- */
/* Solveur principal                                                */
/* ---------------------------------------------------------------- */

SEXP easy_solver_c(SEXP s_x, SEXP s_lb, SEXP s_eqB, SEXP s_ctrl,
                   SEXP s_f, SEXP s_h, SEXP s_rho)
{
    Ctx C;
    int n = LENGTH(s_x), m = LENGTH(s_eqB);
    const double *p = REAL(s_ctrl);

    C.f = s_f;
    C.h = s_h;
    C.rho = s_rho;
    C.grad_eps = p[P_GRAD_EPS];
    C.n = n;
    C.m = m;
    C.lb = REAL(s_lb);
    C.eqB = REAL(s_eqB);
    C.tol = p[P_TOL];
    C.kkt_tol = p[P_KKT_TOL];
    C.stepsize_min = p[P_STEPSIZE_MIN];
    C.stepsize_max = p[P_STEPSIZE_MAX];
    C.max_step = p[P_MAX_STEP];
    C.backtrack = p[P_BACKTRACK];
    C.f_min = p[P_F_MIN];
    C.eq_tol = p[P_EQ_TOL];
    C.active_tol = p[P_ACTIVE_TOL];
    C.curv_tol = p[P_CURV_TOL];
    C.max_eq = (int)p[P_MAX_EQ];

    int maxit = (int)p[P_MAXIT];
    int trace = (int)p[P_TRACE];
    int use_bfgs = p[P_USE_BFGS] != 0.0;
    int bfgs_scale = p[P_BFGS_SCALE] != 0.0;
    double stepsize = p[P_STEPSIZE];

    double *x = ALLOC(n), *g = ALLOC(n), *x_old = ALLOC(n), *g_old = ALLOC(n);
    double *J = ALLOC((size_t)m * n), *J_old = ALLOC((size_t)m * n);
    double *H = ALLOC((size_t)n * n), *d = ALLOC(n);
    double *xtr = ALLOC(n), *xtr2 = ALLOC(n);
    double *r = ALLOC(n), *lambda = ALLOC(m), *hx = ALLOC(m);
    double *s = ALLOC(n), *y = ALLOC(n), *Hy = ALLOC(n);
    double *hist = ALLOC((size_t)7 * (maxit > 0 ? maxit : 1));
    int nh = 0;

    memcpy(x, REAL(s_x), sizeof(double) * (size_t)n);
    for (int j = 0; j < n; j++)
        if (x[j] < C.lb[j])
            x[j] = C.lb[j];

    double f_cur = eval_f(&C, x);
    if (!R_FINITE(f_cur))
        error("Initial objective is not finite");
    if (f_cur < C.f_min)
        error("Initial objective is below f_min");

    if (m > 0)
    {
        if (eq_correction(&C, x, xtr))
        {
            double fc = eval_f(&C, xtr);
            if (R_FINITE(fc) && fc >= C.f_min)
            {
                memcpy(x, xtr, sizeof(double) * (size_t)n);
                f_cur = fc;
            }
        }
    }

    for (int i = 0; i < n * n; i++)
        H[i] = 0.0;
    for (int j = 0; j < n; j++)
        H[j + (size_t)j * n] = 1.0;
    int H_is_id = 1, have_old = 0, have_gj = 0, converged = 0, iter_done = 0;
    const char *reason = "maxit reached";
    Kkt km;

    for (int iter = 1; iter <= maxit; iter++)
    {
        R_CheckUserInterrupt();
        iter_done = iter;

        if (!have_gj)
        {
            eval_gj(&C, x, g, J);
            have_gj = 1;
        }
        eval_h(&C, x, hx);
        kkt_measure(&C, x, g, J, hx, r, lambda, &km);

        double gn = 0.0;
        for (int j = 0; j < n; j++)
            gn += g[j] * g[j];
        gn = sqrt(gn);

        if (trace > 0 && (iter == 1 || iter % trace == 0))
        {
            Rprintf("iter %4d | f = % .6e | kkt = %.3e | eq = %.3e | grad = %.3e | step = %.3e\n",
                    iter, f_cur, km.kkt, km.eq_res, gn, stepsize);
        }

        if (km.kkt <= C.kkt_tol)
        {
            converged = 1;
            reason = "KKT tolerance reached";
            break;
        }

        /* ---- mise a jour BFGS inverse (gradient du Lagrangien) ---- */
        if (use_bfgs && have_old)
        {
            double sy = 0.0, ss = 0.0, yy = 0.0;
            for (int j = 0; j < n; j++)
            {
                s[j] = x[j] - x_old[j];
                y[j] = g[j] - g_old[j];
            }
            if (m > 0)
            {
                for (int j = 0; j < n; j++)
                {
                    double acc = 0.0;
                    for (int i = 0; i < m; i++)
                        acc += (J[i + (size_t)j * m] - J_old[i + (size_t)j * m]) * lambda[i];
                    y[j] += acc;
                }
            }
            for (int j = 0; j < n; j++)
            {
                sy += s[j] * y[j];
                ss += s[j] * s[j];
                yy += y[j] * y[j];
            }

            if (R_FINITE(sy) && sy > C.curv_tol * sqrt(ss) * sqrt(yy))
            {
                if (bfgs_scale && H_is_id)
                {
                    for (int i = 0; i < n * n; i++)
                        H[i] = 0.0;
                    for (int j = 0; j < n; j++)
                        H[j + (size_t)j * n] = sy / yy;
                }
                matvec(H, n, n, y, Hy);
                double yHy = 0.0;
                for (int j = 0; j < n; j++)
                    yHy += y[j] * Hy[j];
                double c1 = (sy + yHy) / (sy * sy);
                for (int j = 0; j < n; j++)
                    for (int i = 0; i < n; i++)
                        H[i + (size_t)j * n] += c1 * s[i] * s[j] - (Hy[i] * s[j] + s[i] * Hy[j]) / sy;
                H_is_id = 0;
            }
        }

        /* ---- direction ---- */
        feasible_direction(&C, x, g, J, use_bfgs ? H : NULL, d);
        double nd = 0.0;
        for (int j = 0; j < n; j++)
            nd += d[j] * d[j];
        nd = sqrt(nd);

        if (!R_FINITE(nd) || nd <= C.tol)
        {
            reason = "No usable feasible descent direction";
            break;
        }

        /* ---- backtracking ---- */
        double alpha = use_bfgs ? C.stepsize_max
                                : (stepsize < C.stepsize_max ? stepsize : C.stepsize_max);
        int accepted = 0;
        double f_trial = f_cur;

        while (alpha >= C.stepsize_min)
        {
            for (int j = 0; j < n; j++)
            {
                xtr[j] = x[j] + alpha * d[j];
                if (xtr[j] < C.lb[j])
                    xtr[j] = C.lb[j];
            }
            f_trial = eval_f(&C, xtr);

            if (!R_FINITE(f_trial) || f_trial < C.f_min)
            {
                alpha *= C.backtrack;
                continue;
            }
            if (f_trial >= f_cur)
            {
                alpha *= C.backtrack;
                continue;
            }

            if (m > 0)
            {
                if (!eq_correction(&C, xtr, xtr2))
                {
                    alpha *= C.backtrack;
                    continue;
                }
                double f2 = eval_f(&C, xtr2);
                if (!R_FINITE(f2) || f2 < C.f_min)
                {
                    alpha *= C.backtrack;
                    continue;
                }
                if (f2 > f_cur)
                {
                    alpha *= C.backtrack;
                    continue;
                }
                memcpy(xtr, xtr2, sizeof(double) * (size_t)n);
                f_trial = f2;
            }
            accepted = 1;
            break;
        }

        if (!accepted)
        {
            add_hist(hist, maxit, &nh, iter, f_cur, km.kkt, km.eq_res, gn, 0.0, 0);
            if (use_bfgs)
            {
                if (H_is_id)
                {
                    reason = "Minimum stepsize reached";
                    break;
                }
                for (int i = 0; i < n * n; i++)
                    H[i] = 0.0;
                for (int j = 0; j < n; j++)
                    H[j + (size_t)j * n] = 1.0;
                H_is_id = 1;
                have_old = 0;
                continue;
            }
            stepsize *= C.backtrack;
            if (stepsize < C.stepsize_min)
            {
                reason = "Minimum stepsize reached";
                break;
            }
            continue;
        }

        /* ---- pas accepte ---- */
        if (use_bfgs)
        {
            memcpy(x_old, x, sizeof(double) * (size_t)n);
            memcpy(g_old, g, sizeof(double) * (size_t)n);
            if (m > 0)
                memcpy(J_old, J, sizeof(double) * (size_t)m * n);
            have_old = 1;
        }
        memcpy(x, xtr, sizeof(double) * (size_t)n);
        f_cur = f_trial;
        have_gj = 0;

        add_hist(hist, maxit, &nh, iter, f_cur, km.kkt, km.eq_res, gn, alpha, 1);

        if (alpha >= 0.9 * stepsize)
        {
            stepsize *= 1.1;
            if (stepsize > C.stepsize_max)
                stepsize = C.stepsize_max;
        }
        if (alpha < 0.9 * stepsize)
        {
            stepsize *= 0.7;
            if (stepsize < C.stepsize_min)
                stepsize = C.stepsize_min;
        }
    }

    /* ---- diagnostics finaux ---- */
    if (!have_gj)
    {
        eval_gj(&C, x, g, J);
        have_gj = 1;
    }
    eval_h(&C, x, hx);
    kkt_measure(&C, x, g, J, hx, r, lambda, &km);

    const char *names[] = {"pars", "value", "convergence", "message", "gradient",
                           "kkt", "kkt_residual", "lambda", "eq", "eq_residual",
                           "active", "iterations", "stepsize", "history", ""};
    SEXP res = PROTECT(mkNamed(VECSXP, names));
    SET_VECTOR_ELT(res, 0, mkvec(x, n));
    SET_VECTOR_ELT(res, 1, ScalarReal(f_cur));
    SET_VECTOR_ELT(res, 2, ScalarInteger(converged ? 0 : 1));
    SET_VECTOR_ELT(res, 3, mkString(reason));
    SET_VECTOR_ELT(res, 4, mkvec(g, n));
    SET_VECTOR_ELT(res, 5, ScalarReal(km.kkt));
    SET_VECTOR_ELT(res, 6, mkvec(r, n));
    SET_VECTOR_ELT(res, 7, mkvec(lambda, m));
    SET_VECTOR_ELT(res, 8, mkvec(hx, m));
    SET_VECTOR_ELT(res, 9, ScalarReal(km.eq_res));

    SEXP act = PROTECT(allocVector(INTSXP, km.nact));
    int q = 0;
    for (int j = 0; j < n; j++)
        if (R_FINITE(C.lb[j]) && x[j] <= C.lb[j] + C.active_tol)
            INTEGER(act)
    [q++] = j + 1;
    SET_VECTOR_ELT(res, 10, act);
    UNPROTECT(1);

    SET_VECTOR_ELT(res, 11, ScalarInteger(iter_done));
    SET_VECTOR_ELT(res, 12, ScalarReal(stepsize));

    SEXP hm = PROTECT(allocMatrix(REALSXP, nh, 7));
    for (int col = 0; col < 7; col++)
        for (int row = 0; row < nh; row++)
            REAL(hm)
    [row + (size_t)col * nh] = hist[row + (size_t)col * maxit];
    SET_VECTOR_ELT(res, 13, hm);
    UNPROTECT(2);
    return res;
}