/*
 * gPlace placement engine (libgplace_core.so).
 *
 * Exact re-implementation of the numerical core of pplacer v1.1.alpha19 (commit 807f6f3, the
 * build GTDB-Tk uses) for its Gcat (FastTree CAT) and Gmix (discrete gamma) models, protein and
 * DNA: every floating-point result keeps pplacer's operands, operation order and rounding, so
 * placements are bit-identical. Built twice: NS=20 (protein, libgplace_core.so) and NS=4 (DNA,
 * libgplace_core_nt.so).
 *
 *   - dediagonalize / multi_exp            (pplacer_src/linear_c.c, diagd.ml)
 *   - evolve_into = per-site gsl dgemv      (gcat_model.ml, gslcblas loop order)
 *   - pairwise / statd-pairwise products    (linear_c.c)
 *   - perhaps_pull_exponent                 (gcat_model.ml)
 *   - mat_log_like3 / mat_masked_logdot     (linear_c.c)
 *   - Minimization.brent + bisection start  (minimization.ml) driving GSL 1.16's Brent
 *   - Three_tax.optimize                    (three_tax.ml; the tree keeps the LAST evaluated
 *                                            branch lengths; an assertion inside a GSL callback
 *                                            becomes NaN -> GSL error -> retry, as ocamlgsl does)
 *
 * SIMD (AVX2) lanes only ever process DIFFERENT results side by side; no reduction is
 * reordered and mul/add stay separate (build with -ffp-contract=off). The scalar originals are
 * kept as *_ref for self-tests. Link with GSL 1.16, glibc-2.19 hypot and log1p, and CORE-MATH
 * exp/log/pow/sin (see gplace/engine/build.py).
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 */
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <limits.h>

/* The number of states is a compile-time parameter (-DNS=4 for DNA). The AVX2 kernels are
 * written for 20 states (5 registers of 4 lanes); other state counts use the scalar forms below,
 * which perform the same operations in the same order (lanes never reorder a reduction).
 * -DXP_FORCE_SCALAR builds a 20-state engine on the scalar forms (to check them). */
#ifndef NS
#define NS 20
#endif
/* AVX2 kernels only when the compiler targets AVX2 (the x86-64 "avx2" engine, -mavx2 -mfma); NEON
 * kernels on ARM64 (the "neon" engine; Advanced SIMD is part of every ARMv8-A CPU); otherwise (the
 * "portable" engine, -DXP_FORCE_SCALAR on ARM64) the scalar forms. All give the same results: lanes
 * hold different results, each computed with the same operations in the same order (a multiply
 * and an add, never a fused multiply-add), and no reduction is reordered. */
#if defined(__AVX2__) && !defined(XP_FORCE_SCALAR)
#include <immintrin.h>
#define XP_AVX2 1
#else
#define XP_AVX2 0
#endif
#if defined(__aarch64__) && defined(__ARM_NEON) && !defined(XP_FORCE_SCALAR)
#include <arm_neon.h>
#define XP_NEON 1
#else
#define XP_NEON 0
#endif
#if NS == 20 && (XP_AVX2 || XP_NEON)
#define XP_SIMD20 1
#else
#define XP_SIMD20 0
#endif

/* XP_CR: use correctly-rounded exp/log (CORE-MATH). Stock pplacer is statically linked with
 * glibc 2.19, whose IBM-derived exp/log are correctly rounded; modern glibc's are not always. */
#ifdef XP_CR
double cr_exp(double); double cr_log(double);
#define exp cr_exp
#define log cr_log
#endif
double xp_exp(double v) { return exp(v); }
double xp_log(double v) { return log(v); }
/* pow and sin are only used by GSL (gamma functions); the link maps GSL's exp/log/pow/sin to the
 * correctly-rounded CORE-MATH functions (glibc 2.19's are correctly rounded): see build.py. */
double xp_pow(double a, double b) { return pow(a, b); }
double xp_sin(double v) { return sin(v); }
double xp_hypot(double a, double b) { return hypot(a, b); }
/* log1p is only used by GSL (gsl_cdf_gamma_Pinv, 21+ gamma categories): glibc 2.19's s_log1p.c is
 * linked in (build.py); exported for the self-tests. */
double xp_log1p(double v) { return log1p(v); }
#define XP_STR2(x) #x
#define XP_STR(x) XP_STR2(x)
#if XP_AVX2
#define XP_KERNELS "AVX2"
#elif XP_NEON
#define XP_KERNELS "NEON"
#else
#define XP_KERNELS "portable, scalar"
#endif
#if defined(__x86_64__)
#define XP_ARCH "x86-64"
#elif defined(__aarch64__)
#define XP_ARCH "arm64"
#else
#define XP_ARCH "other"
#endif
const char *xp_engine_version(void) { return "gplace-core 4 (" XP_STR(NS) " states): pplacer v1.1.alpha19 exact (GSL 1.16, glibc-2.19 hypot/log1p, CORE-MATH exp/log, " XP_ARCH " " XP_KERNELS "; C)"; }
int xp_simd(void) { return XP_AVX2 ? 1 : XP_NEON ? 2 : 0; }   /* 0 scalar, 1 AVX2, 2 NEON */
int xp_n_states(void) { return NS; }

/* ------------------------------------------------------------------ model -- */

#include <stdio.h>

/* pplacer turns every GSL error into an OCaml exception (Gsl_error.init () in ppatteries.ml), and its
 * discrete gamma catches none: the first GSL error raised while the rates are computed ends the run
 * (e.g. Gsl_exn(6, "inverse failed to converge") for alpha 1e6). The engine records that first error
 * so that gPlace stops with the same message; outside the rates GSL errors are ignored (the optimiser
 * uses GSL's return codes, as pplacer's handlers of Gsl_exn do). The rates are computed before any
 * worker starts, so nothing else calls GSL meanwhile. */
static int xp_gsl_recording = 0;
static int xp_gsl_errno = 0;
static char xp_gsl_reason[256];

static void xp_gsl_error_handler(const char *reason, const char *file, int line, int gsl_errno)
{
  (void)file; (void)line;
  if (xp_gsl_recording && !xp_gsl_errno) {
    xp_gsl_errno = gsl_errno;
    snprintf(xp_gsl_reason, sizeof xp_gsl_reason, "%s", reason ? reason : "");
  }
}

/* GSL must report errors as return codes (the optimiser relies on them), not abort. The handler
 * ignores them, like gsl_set_error_handler_off, except while the discrete gamma rates are computed:
 * there it records the first one. */
int xp_init(void) { gsl_set_error_handler(xp_gsl_error_handler); return 0; }

/* The first GSL error of the last xp_discrete_gamma call: its gsl_errno (0 if none), with its message
 * copied into buf (NUL-terminated, truncated to len bytes). */
int xp_gamma_gsl_error(char *buf, size_t len)
{
  if (!xp_gsl_errno) return 0;
  if (len > 0) {
    size_t n = strlen(xp_gsl_reason);
    if (n > len - 1) n = len - 1;
    memcpy(buf, xp_gsl_reason, n);
    buf[n] = 0;
  }
  return xp_gsl_errno;
}

/* Targets of the libm functions GSL references but gPlace never calls (build.py LIBM_TRAPS): abort
 * instead of silently using an inexact system implementation. */
#define XP_TRAP(f) void xp_libm_trap_##f(void) { fputs("gPlace engine: unexpected call to libm " #f \
  "; results would not be bit-exact. Please report this.\n", stderr); abort(); }
XP_TRAP(cos) XP_TRAP(cosh) XP_TRAP(sinh) XP_TRAP(tan) XP_TRAP(tanh) XP_TRAP(log1p) XP_TRAP(sincos)
XP_TRAP(expm1) XP_TRAP(atan) XP_TRAP(atan2) XP_TRAP(asin) XP_TRAP(acos) XP_TRAP(log10) XP_TRAP(log2)
XP_TRAP(exp2) XP_TRAP(cbrt) XP_TRAP(erf) XP_TRAP(erfc) XP_TRAP(lgamma) XP_TRAP(tgamma) XP_TRAP(fma)
XP_TRAP(acosh) XP_TRAP(atanh) XP_TRAP(asinh)

/* Gsl_eigen.symmv: unsorted eigen-decomposition of a symmetric matrix (lower triangle used). */
int xp_symmv(const double *m_in, double *eval, double *evec, int n)
{
  gsl_matrix *m = gsl_matrix_alloc(n, n);
  gsl_vector *ev = gsl_vector_alloc(n);
  gsl_matrix *vec = gsl_matrix_alloc(n, n);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
  memcpy(m->data, m_in, sizeof(double) * n * n);
  int st = gsl_eigen_symmv(m, ev, vec, w);
  memcpy(eval, ev->data, sizeof(double) * n);
  memcpy(evec, vec->data, sizeof(double) * n * n);
  gsl_eigen_symmv_free(w); gsl_matrix_free(vec); gsl_vector_free(ev); gsl_matrix_free(m);
  return st;
}

/* Reference (scalar) kernels from xp.c: row-major P[i][j]. Kept for verification only. */
/* dst.{i,j} = sum_k lambda.{k} * u.{i,k} * uit.{j,k}   (linear_c.c dediagonalize) */
void xp_dediag_ref(double *dst, const double *u, const double *lambda, const double *uit, int n)
{
  int i, j, k;
  for (i = 0; i < n; i++) {
    const double *uit_p = uit;
    for (j = 0; j < n; j++) {
      double d = 0;
      for (k = 0; k < n; k++) d += lambda[k] * u[i * n + k] * uit_p[k];
      dst[i * n + j] = d;
      uit_p += n;
    }
  }
}

/* Diagd.multi_exp: tensor[r] = X exp(rates[r]*bl*Lambda) X^-1 for rates with mask[r]. */
void xp_multi_exp_ref(double *tensor, const double *x, const double *l, const double *xit,
                  const double *rates, int nrates, const uint8_t *mask, double bl)
{
  double util[NS];
  memset(tensor, 0, sizeof(double) * nrates * NS * NS);
  for (int r = 0; r < nrates; r++) {
    if (!mask[r]) continue;
    for (int i = 0; i < NS; i++) util[i] = exp(rates[r] * bl * l[i]);
    xp_dediag_ref(tensor + (size_t)r * NS * NS, x, util, xit, NS);
  }
}

/* Gcat evolve_into: per site, dst = P_{cat(site)} * src via gslcblas dgemv (NoTrans):
 *   temp += X[j] * A[i,j];  Y[i] = 0 + alpha*temp.  Exponents are copied by the caller. */
void xp_evolve_ref(double *dst, const double *src, int nsites, const int32_t *cat, const double *tensor)
{
  for (int s = 0; s < nsites; s++) {
    const double *P = tensor + (size_t)cat[s] * NS * NS;
    const double *xv = src + (size_t)s * NS;
    double *yv = dst + (size_t)s * NS;
    for (int i = 0; i < NS; i++) {
      double temp = 0.0;
      for (int j = 0; j < NS; j++) temp += xv[j] * P[i * NS + j];
      yv[i] = 0.0 + 1.0 * temp;
    }
  }
}


/* dediagonalize, TRANSPOSED output: dstT[j*n+i] = P[i][j] = sum_k (lambda[k]*u[i*n+k]) * uit[j*n+k],
 * k ascending from +0.0 -- exactly the scalar expression d += lambda[k] * u[i*n+k] * uit_p[k].
 * The product lambda[k]*u[i][k] does not depend on j, so it is computed once (same value).
 * Lanes run over i (4 different P[i][j] at once). */
void xp_dediag(double *dstT, const double *u, const double *lambda, const double *uit, int n)
{
  if (n != NS || !XP_SIMD20) {                     /* generic scalar path, transposed store */
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
        double d = 0;
        for (int k = 0; k < n; k++) d += lambda[k] * u[i * n + k] * uit[j * n + k];
        dstT[j * n + i] = d;
      }
    return;
  }
#if XP_SIMD20
  double A[NS * NS];                               /* A[k][i] = lambda[k] * u[i][k] */
  for (int k = 0; k < NS; k++)
    for (int i = 0; i < NS; i++) A[k * NS + i] = lambda[k] * u[i * NS + k];
  for (int j = 0; j < NS; j++) {
    const double *uj = uit + j * NS;
#if XP_AVX2
    __m256d a0 = _mm256_setzero_pd(), a1 = a0, a2 = a0, a3 = a0, a4 = a0;
    for (int k = 0; k < NS; k++) {
      const __m256d w = _mm256_set1_pd(uj[k]);
      const double *Ak = A + k * NS;
      a0 = _mm256_add_pd(a0, _mm256_mul_pd(_mm256_loadu_pd(Ak + 0), w));
      a1 = _mm256_add_pd(a1, _mm256_mul_pd(_mm256_loadu_pd(Ak + 4), w));
      a2 = _mm256_add_pd(a2, _mm256_mul_pd(_mm256_loadu_pd(Ak + 8), w));
      a3 = _mm256_add_pd(a3, _mm256_mul_pd(_mm256_loadu_pd(Ak + 12), w));
      a4 = _mm256_add_pd(a4, _mm256_mul_pd(_mm256_loadu_pd(Ak + 16), w));
    }
    double *dj = dstT + j * NS;
    _mm256_storeu_pd(dj + 0, a0); _mm256_storeu_pd(dj + 4, a1); _mm256_storeu_pd(dj + 8, a2);
    _mm256_storeu_pd(dj + 12, a3); _mm256_storeu_pd(dj + 16, a4);
#else                                              /* NEON: lanes run over i, 2 per register */
    float64x2_t a0 = vdupq_n_f64(0.0), a1 = a0, a2 = a0, a3 = a0, a4 = a0, a5 = a0, a6 = a0, a7 = a0,
                a8 = a0, a9 = a0;
    for (int k = 0; k < NS; k++) {
      const float64x2_t w = vdupq_n_f64(uj[k]);
      const double *Ak = A + k * NS;
#define XP_NDD(r) a##r = vaddq_f64(a##r, vmulq_f64(vld1q_f64(Ak + 2 * r), w))
      XP_NDD(0); XP_NDD(1); XP_NDD(2); XP_NDD(3); XP_NDD(4); XP_NDD(5); XP_NDD(6); XP_NDD(7); XP_NDD(8); XP_NDD(9);
#undef XP_NDD
    }
    double *dj = dstT + j * NS;
    vst1q_f64(dj + 0, a0); vst1q_f64(dj + 2, a1); vst1q_f64(dj + 4, a2); vst1q_f64(dj + 6, a3);
    vst1q_f64(dj + 8, a4); vst1q_f64(dj + 10, a5); vst1q_f64(dj + 12, a6); vst1q_f64(dj + 14, a7);
    vst1q_f64(dj + 16, a8); vst1q_f64(dj + 18, a9);
#endif
  }
#endif
}

/* Diagd.multi_exp: tensorT[r] = (X exp(rates[r]*bl*Lambda) X^-1)^T for rates with mask[r]. */
void xp_multi_exp(double *tensor, const double *x, const double *l, const double *xit,
                  const double *rates, int nrates, const uint8_t *mask, double bl)
{
  double util[NS];
  memset(tensor, 0, sizeof(double) * nrates * NS * NS);
  for (int r = 0; r < nrates; r++) {
    if (!mask[r]) continue;
    for (int i = 0; i < NS; i++) util[i] = exp(rates[r] * bl * l[i]);
    xp_dediag(tensor + (size_t)r * NS * NS, x, util, xit, NS);
  }
}

/* ------------------------------------------------------------------- glvs -- */

/* Gcat evolve_into with the TRANSPOSED tensor from xp_multi_exp: per site and output state i,
 *   temp = +0.0; for j ascending: temp = temp + x[j]*P[i][j];  y[i] = 0.0 + temp
 * (gslcblas dgemv NoTrans with alpha=1, beta=0; 1.0*temp == temp exactly).
 * Lanes run over i: 4 different y[i] per instruction. Exponents are copied by the caller. */
static inline void evolve_site(double *yv, const double *xv, const double *PT)
{
#if !XP_SIMD20
  for (int i = 0; i < NS; i++) {                   /* same as the lanes: temp = +0.0 + x[j]*P[i][j]... */
    double temp = 0.0;
    for (int j = 0; j < NS; j++) temp = temp + xv[j] * PT[j * NS + i];
    yv[i] = 0.0 + temp;
  }
#elif XP_AVX2
  const __m256d z = _mm256_setzero_pd();
  __m256d a0 = z, a1 = z, a2 = z, a3 = z, a4 = z;
  for (int j = 0; j < NS; j++) {
    const __m256d xj = _mm256_set1_pd(xv[j]);
    const double *Pj = PT + j * NS;
    a0 = _mm256_add_pd(a0, _mm256_mul_pd(xj, _mm256_loadu_pd(Pj + 0)));
    a1 = _mm256_add_pd(a1, _mm256_mul_pd(xj, _mm256_loadu_pd(Pj + 4)));
    a2 = _mm256_add_pd(a2, _mm256_mul_pd(xj, _mm256_loadu_pd(Pj + 8)));
    a3 = _mm256_add_pd(a3, _mm256_mul_pd(xj, _mm256_loadu_pd(Pj + 12)));
    a4 = _mm256_add_pd(a4, _mm256_mul_pd(xj, _mm256_loadu_pd(Pj + 16)));
  }
  _mm256_storeu_pd(yv + 0, _mm256_add_pd(z, a0)); _mm256_storeu_pd(yv + 4, _mm256_add_pd(z, a1));
  _mm256_storeu_pd(yv + 8, _mm256_add_pd(z, a2)); _mm256_storeu_pd(yv + 12, _mm256_add_pd(z, a3));
  _mm256_storeu_pd(yv + 16, _mm256_add_pd(z, a4));
#else                                              /* NEON: 2 different y[i] per register */
  const float64x2_t z = vdupq_n_f64(0.0);
  float64x2_t a0 = z, a1 = z, a2 = z, a3 = z, a4 = z, a5 = z, a6 = z, a7 = z, a8 = z, a9 = z;
  for (int j = 0; j < NS; j++) {
    const float64x2_t xj = vdupq_n_f64(xv[j]);
    const double *Pj = PT + j * NS;
#define XP_NEV(r) a##r = vaddq_f64(a##r, vmulq_f64(xj, vld1q_f64(Pj + 2 * r)))
    XP_NEV(0); XP_NEV(1); XP_NEV(2); XP_NEV(3); XP_NEV(4); XP_NEV(5); XP_NEV(6); XP_NEV(7); XP_NEV(8); XP_NEV(9);
#undef XP_NEV
  }
  vst1q_f64(yv + 0, vaddq_f64(z, a0)); vst1q_f64(yv + 2, vaddq_f64(z, a1)); vst1q_f64(yv + 4, vaddq_f64(z, a2));
  vst1q_f64(yv + 6, vaddq_f64(z, a3)); vst1q_f64(yv + 8, vaddq_f64(z, a4)); vst1q_f64(yv + 10, vaddq_f64(z, a5));
  vst1q_f64(yv + 12, vaddq_f64(z, a6)); vst1q_f64(yv + 14, vaddq_f64(z, a7)); vst1q_f64(yv + 16, vaddq_f64(z, a8));
  vst1q_f64(yv + 18, vaddq_f64(z, a9));
#endif
}

void xp_evolve(double *dst, const double *src, int nsites, const int32_t *cat, const double *tensor)
{
  for (int s = 0; s < nsites; s++)
    evolve_site(dst + (size_t)s * NS, src + (size_t)s * NS, tensor + (size_t)cat[s] * NS * NS);
}

/* dst = x .* y (elementwise), exponents summed */
void xp_pairwise(double *dst, int64_t *de, const double *x, const int64_t *xe,
                 const double *y, const int64_t *ye, int nsites)
{
  size_t n = (size_t)nsites * NS;
  for (size_t i = 0; i < n; i++) dst[i] = x[i] * y[i];
  for (int s = 0; s < nsites; s++) de[s] = xe[s] + ye[s];
}

/* dst = statd .* a .* b   (mat_statd_pairwise_prod), exponents summed */
void xp_statd_pairwise(double *dst, int64_t *de, const double *statd, const double *a,
                       const int64_t *ae, const double *b, const int64_t *be, int nsites)
{
  for (int s = 0; s < nsites; s++) {
    for (int st = 0; st < NS; st++) {
      *dst = statd[st] * (*a) * (*b);
      dst++; a++; b++;
    }
    de[s] = ae[s] + be[s];
  }
}

/* Glv.perhaps_pull_exponent: per site, if frexp-exponent of the max entry < min_twoexp,
 * scale the site by 2^-twoexp and add twoexp to the stored exponent. */
void xp_pull_exponent(double *a, int64_t *e, int nsites, int min_twoexp)
{
  for (int s = 0; s < nsites; s++) {
    double *v = a + (size_t)s * NS;
    double mx = v[0];
    for (int i = 1; i < NS; i++) if (v[i] > mx) mx = v[i];
    int twoexp;
    frexp(mx, &twoexp);
    if (twoexp < min_twoexp) {
      double f = ldexp(0.5, -twoexp + 1);
      for (int i = 0; i < NS; i++) v[i] *= f;
      e[s] += twoexp;
    }
  }
}

double xp_total_twoexp(const int64_t *e, int nsites)
{
  double tot = 0.;
  for (int s = 0; s < nsites; s++) tot += (double)e[s];
  return tot;
}

/* mat_log_like3: sum_s log(sum_state statd*x*y*z) */
double xp_log_like3(const double *statd, const double *x, const double *y, const double *z, int nsites)
{
  double ll_tot = 0;
  for (int s = 0; s < nsites; s++) {
    double util = 0;
    for (int st = 0; st < NS; st++) util += statd[st] * x[st] * y[st] * z[st];
    ll_tot += log(util);
    x += NS; y += NS; z += NS;
  }
  return ll_tot;
}

/* Exact edge-major prefilter for many queries at once.
 * pplacer scores query q at edge v as  sum_{s in mask_q} log( sum_st qev_q[s][st] * snode[s][st] )
 *                                      + log2*(0 + sum_{s in mask_q} snode_e[s]).
 * qev_q[s] depends only on (site category, query character code), so we first build
 *   W[s][c] = log( sum_st qevtab[cat[s]][c][st] * snode[s][st] )   (identical operands/order)
 * and then each query is a sum of W lookups in site order (bit-identical to pplacer).
 * codes[q*nsites+s] = code index, or 255 for non-informative ('-','?').
 * out[q*stride] receives the score. */
void xp_edge_scores(double *out, long stride, const double *snode, const int64_t *se,
                    const double *qevtab, const int32_t *cat, int nsites, int ncode,
                    const uint8_t *codes, int nq, double *W)
{
  for (int s = 0; s < nsites; s++) {
    const double *sn = snode + (size_t)s * NS;
    const double *qt = qevtab + (size_t)cat[s] * ncode * NS;
    for (int c = 0; c < ncode; c++) {
      const double *qv = qt + (size_t)c * NS;
      double util = 0;
      for (int st = 0; st < NS; st++) util += qv[st] * sn[st];
      W[(size_t)s * ncode + c] = log(util);
    }
  }
  for (int q = 0; q < nq; q++) {
    const uint8_t *cq = codes + (size_t)q * nsites;
    double ll = 0; long long et = 0;
    for (int s = 0; s < nsites; s++) {
      uint8_t c = cq[s];
      if (c != 255) { ll += W[(size_t)s * ncode + c]; et += se[s]; }
    }
    out[(size_t)q * stride] = ll + (M_LN2 * (0. + (double)et));
  }
}

/* Exact prefilter for 8 edges at once (lanes = edges). For each lane (edge) and query q the
 * score is built exactly as xp_edge_scores does it:
 *   W[s][c] = log( sum_st qevtab[cat[s]][c][st] * snode[s][st] )     (same operands/order)
 *   ll = +0.0; for s ascending with code != gap: ll = ll + W[s][code];  et = sum se[s] (integer)
 *   score = ll + log2*(0 + (double)et)
 * Gap sites add +0.0 instead of being skipped: ll can never be -0.0 (log never returns -0.0
 * and x + (-x) = +0.0 under round-to-nearest), so ll + 0.0 == ll bit for bit.
 * codes[q*nsites+s] in 0..ncode-1, or ncode for a gap ('-' or '?').
 * out[q*ldo + col0 + lane] receives the score for lane < ne.
 * Site blocks keep the W block in L2; 4 queries are interleaved to overlap the add latency.
 * Exponents are summed in int32 lanes within one 128-site block, then added into int64 totals:
 * on very large trees (the full GTDB bacterial tree) per-site exponents reach ~2^18, which the
 * int32 totals of earlier versions could not hold. Integer addition is exact in any order, so
 * the result equals the scalar int64 sum.
 * Returns 0, or -1 on bad arguments. */
#define XP_EB 8
#if XP_AVX2
/* add 8 int32 lane sums into 8 int64 totals (exact) */
#define ADD64(E, e) do { \
    __m256i lo_ = _mm256_cvtepi32_epi64(_mm256_castsi256_si128(e)); \
    __m256i hi_ = _mm256_cvtepi32_epi64(_mm256_extracti128_si256(e, 1)); \
    _mm256_storeu_si256((__m256i *)(E), _mm256_add_epi64(_mm256_loadu_si256((const __m256i *)(E)), lo_)); \
    _mm256_storeu_si256((__m256i *)((E) + 4), _mm256_add_epi64(_mm256_loadu_si256((const __m256i *)((E) + 4)), hi_)); \
  } while (0)
#endif
#if XP_NEON
/* add 8 int32 lane sums (lo: lanes 0-3, hi: 4-7) into 8 int64 totals (exact) */
#define NADD64(E, lo, hi) do { \
    vst1q_s64((E) + 0, vaddq_s64(vld1q_s64((E) + 0), vmovl_s32(vget_low_s32(lo)))); \
    vst1q_s64((E) + 2, vaddq_s64(vld1q_s64((E) + 2), vmovl_high_s32(lo))); \
    vst1q_s64((E) + 4, vaddq_s64(vld1q_s64((E) + 4), vmovl_s32(vget_low_s32(hi)))); \
    vst1q_s64((E) + 6, vaddq_s64(vld1q_s64((E) + 6), vmovl_high_s32(hi))); \
  } while (0)
#endif
static void es8_scalar(double *out, long ldo, long col0, int ne, const double *W,
                       const int64_t *se8, int nsites, int ncode, const uint8_t *codes, int nq)
{
  const int NC = ncode + 1;
  for (int lane = 0; lane < ne; lane++)
    for (int q = 0; q < nq; q++) {
      const uint8_t *cq = codes + (size_t)q * nsites;
      double ll = 0; long long et = 0;
      for (int s = 0; s < nsites; s++) {
        unsigned c = cq[s];
        if (c != (unsigned)ncode) { ll += W[((size_t)s * NC + c) * XP_EB + lane]; et += se8[(size_t)lane * nsites + s]; }
      }
      out[(size_t)q * ldo + col0 + lane] = ll + (M_LN2 * (0. + (double)et));
    }
}

int xp_edge_scores8(double *out, long ldo, long col0, int ne,
                    const double *sn8, const int64_t *se8,
                    const double *qevtab, const int32_t *cat, int nsites, int ncode,
                    const uint8_t *codes, int nq)
{
  if (ne < 1 || ne > XP_EB || ncode < 1 || ncode > 254) return -1;
  const int NC = ncode + 1;
  double *W = (double *)aligned_alloc(64, sizeof(double) * (size_t)nsites * NC * XP_EB);
  int32_t *SE = (int32_t *)aligned_alloc(64, sizeof(int32_t) * (size_t)nsites * 2 * XP_EB);
  double *accD = (double *)aligned_alloc(64, sizeof(double) * ((size_t)nq * XP_EB + 8));
  int64_t *accE = (int64_t *)aligned_alloc(64, sizeof(int64_t) * ((size_t)nq * XP_EB + 8));
  if (!W || !SE || !accD || !accE) { free(W); free(SE); free(accD); free(accE); return -1; }
  int small = 1;       /* int32 sums over one 128-site block cannot overflow (|se| <= 2^23) */
  for (int s = 0; s < nsites; s++) {
    const double *qt = qevtab + (size_t)cat[s] * ncode * NS;
    double *Ws = W + (size_t)s * NC * XP_EB;
    for (int c = 0; c < ncode; c++) {
      const double *qv = qt + (size_t)c * NS;
      for (int lane = 0; lane < XP_EB; lane++) {
        if (lane < ne) {
          const double *sn = sn8 + ((size_t)lane * nsites + s) * NS;
          double util = 0;
          for (int st = 0; st < NS; st++) util += qv[st] * sn[st];
          Ws[c * XP_EB + lane] = log(util);
        } else Ws[c * XP_EB + lane] = 0.0;
      }
    }
    for (int lane = 0; lane < XP_EB; lane++) Ws[ncode * XP_EB + lane] = 0.0;
    int32_t *SEs = SE + (size_t)s * 2 * XP_EB;
    for (int lane = 0; lane < XP_EB; lane++) {
      int64_t v = lane < ne ? se8[(size_t)lane * nsites + s] : 0;
      if (v > (1 << 23) || v < -(1 << 23)) small = 0;
      SEs[lane] = 0; SEs[XP_EB + lane] = (int32_t)v;
    }
  }
  if (!small || !(XP_AVX2 || XP_NEON)) {   /* exact scalar form: portable engine, or |se| > 2^23 */
    es8_scalar(out, ldo, col0, ne, W, se8, nsites, ncode, codes, nq);
    free(W); free(SE); free(accD); free(accE); return 0;
  }
#if XP_AVX2
  memset(accD, 0, sizeof(double) * (size_t)nq * XP_EB);      /* +0.0 */
  memset(accE, 0, sizeof(int64_t) * (size_t)nq * XP_EB);
  const unsigned GAP = (unsigned)ncode;
  const int BS = 128;
  for (int s0 = 0; s0 < nsites; s0 += BS) {
    const int s1 = s0 + BS < nsites ? s0 + BS : nsites;
    int q = 0;
    for (; q + 4 <= nq; q += 4) {
      double *A0 = accD + (size_t)q * XP_EB, *A1 = A0 + XP_EB, *A2 = A1 + XP_EB, *A3 = A2 + XP_EB;
      int64_t *E0 = accE + (size_t)q * XP_EB, *E1 = E0 + XP_EB, *E2 = E1 + XP_EB, *E3 = E2 + XP_EB;
      __m256d a0l = _mm256_loadu_pd(A0), a0h = _mm256_loadu_pd(A0 + 4);
      __m256d a1l = _mm256_loadu_pd(A1), a1h = _mm256_loadu_pd(A1 + 4);
      __m256d a2l = _mm256_loadu_pd(A2), a2h = _mm256_loadu_pd(A2 + 4);
      __m256d a3l = _mm256_loadu_pd(A3), a3h = _mm256_loadu_pd(A3 + 4);
      __m256i e0 = _mm256_setzero_si256(), e1 = e0, e2 = e0, e3 = e0;   /* block partial sums */
      const uint8_t *c0 = codes + (size_t)q * nsites, *c1 = c0 + nsites, *c2 = c1 + nsites, *c3 = c2 + nsites;
      for (int s = s0; s < s1; s++) {
        const double *Ws = W + (size_t)s * NC * XP_EB;
        const int32_t *SEs = SE + (size_t)s * 2 * XP_EB;
        const unsigned k0 = c0[s], k1 = c1[s], k2 = c2[s], k3 = c3[s];
        const double *w0 = Ws + k0 * XP_EB, *w1 = Ws + k1 * XP_EB, *w2 = Ws + k2 * XP_EB, *w3 = Ws + k3 * XP_EB;
        a0l = _mm256_add_pd(a0l, _mm256_loadu_pd(w0)); a0h = _mm256_add_pd(a0h, _mm256_loadu_pd(w0 + 4));
        a1l = _mm256_add_pd(a1l, _mm256_loadu_pd(w1)); a1h = _mm256_add_pd(a1h, _mm256_loadu_pd(w1 + 4));
        a2l = _mm256_add_pd(a2l, _mm256_loadu_pd(w2)); a2h = _mm256_add_pd(a2h, _mm256_loadu_pd(w2 + 4));
        a3l = _mm256_add_pd(a3l, _mm256_loadu_pd(w3)); a3h = _mm256_add_pd(a3h, _mm256_loadu_pd(w3 + 4));
        e0 = _mm256_add_epi32(e0, _mm256_loadu_si256((const __m256i *)(SEs + (k0 != GAP) * XP_EB)));
        e1 = _mm256_add_epi32(e1, _mm256_loadu_si256((const __m256i *)(SEs + (k1 != GAP) * XP_EB)));
        e2 = _mm256_add_epi32(e2, _mm256_loadu_si256((const __m256i *)(SEs + (k2 != GAP) * XP_EB)));
        e3 = _mm256_add_epi32(e3, _mm256_loadu_si256((const __m256i *)(SEs + (k3 != GAP) * XP_EB)));
      }
      _mm256_storeu_pd(A0, a0l); _mm256_storeu_pd(A0 + 4, a0h); _mm256_storeu_pd(A1, a1l); _mm256_storeu_pd(A1 + 4, a1h);
      _mm256_storeu_pd(A2, a2l); _mm256_storeu_pd(A2 + 4, a2h); _mm256_storeu_pd(A3, a3l); _mm256_storeu_pd(A3 + 4, a3h);
      ADD64(E0, e0); ADD64(E1, e1); ADD64(E2, e2); ADD64(E3, e3);
    }
    for (; q < nq; q++) {
      double *A0 = accD + (size_t)q * XP_EB; int64_t *E0 = accE + (size_t)q * XP_EB;
      __m256d a0l = _mm256_loadu_pd(A0), a0h = _mm256_loadu_pd(A0 + 4);
      __m256i e0 = _mm256_setzero_si256();
      const uint8_t *c0 = codes + (size_t)q * nsites;
      for (int s = s0; s < s1; s++) {
        const double *w0 = W + ((size_t)s * NC + c0[s]) * XP_EB;
        const int32_t *SEs = SE + (size_t)s * 2 * XP_EB;
        a0l = _mm256_add_pd(a0l, _mm256_loadu_pd(w0)); a0h = _mm256_add_pd(a0h, _mm256_loadu_pd(w0 + 4));
        e0 = _mm256_add_epi32(e0, _mm256_loadu_si256((const __m256i *)(SEs + (c0[s] != GAP) * XP_EB)));
      }
      _mm256_storeu_pd(A0, a0l); _mm256_storeu_pd(A0 + 4, a0h); ADD64(E0, e0);
    }
  }
  for (int q = 0; q < nq; q++)
    for (int lane = 0; lane < ne; lane++)
      out[(size_t)q * ldo + col0 + lane] =
        accD[(size_t)q * XP_EB + lane] + (M_LN2 * (0. + (double)accE[(size_t)q * XP_EB + lane]));
#elif XP_NEON
  /* as the AVX2 loop: per (query, lane) the W values are added in site order from +0.0 and the
   * exponents summed exactly; 8 lanes = 4 registers per query, 4 queries interleaved */
  memset(accD, 0, sizeof(double) * (size_t)nq * XP_EB);      /* +0.0 */
  memset(accE, 0, sizeof(int64_t) * (size_t)nq * XP_EB);
  const unsigned GAP = (unsigned)ncode;
  const int BS = 128;
  for (int s0 = 0; s0 < nsites; s0 += BS) {
    const int s1 = s0 + BS < nsites ? s0 + BS : nsites;
    int q = 0;
    for (; q + 4 <= nq; q += 4) {
      double *A = accD + (size_t)q * XP_EB;                     /* 4 x 8 accumulators */
      int64_t *E = accE + (size_t)q * XP_EB;
#define XP_NLD(Q) float64x2_t a##Q##0 = vld1q_f64(A + 8 * Q), a##Q##1 = vld1q_f64(A + 8 * Q + 2), \
                              a##Q##2 = vld1q_f64(A + 8 * Q + 4), a##Q##3 = vld1q_f64(A + 8 * Q + 6); \
                  int32x4_t e##Q##l = vdupq_n_s32(0), e##Q##h = vdupq_n_s32(0)
      XP_NLD(0); XP_NLD(1); XP_NLD(2); XP_NLD(3);
#undef XP_NLD
      const uint8_t *c0 = codes + (size_t)q * nsites, *c1 = c0 + nsites, *c2 = c1 + nsites, *c3 = c2 + nsites;
      for (int s = s0; s < s1; s++) {
        const double *Ws = W + (size_t)s * NC * XP_EB;
        const int32_t *SEs = SE + (size_t)s * 2 * XP_EB;
        const unsigned k0 = c0[s], k1 = c1[s], k2 = c2[s], k3 = c3[s];
#define XP_NAC(Q) do { const double *w_ = Ws + k##Q * XP_EB; \
          const int32_t *x_ = SEs + (k##Q != GAP) * XP_EB; \
          a##Q##0 = vaddq_f64(a##Q##0, vld1q_f64(w_)); a##Q##1 = vaddq_f64(a##Q##1, vld1q_f64(w_ + 2)); \
          a##Q##2 = vaddq_f64(a##Q##2, vld1q_f64(w_ + 4)); a##Q##3 = vaddq_f64(a##Q##3, vld1q_f64(w_ + 6)); \
          e##Q##l = vaddq_s32(e##Q##l, vld1q_s32(x_)); e##Q##h = vaddq_s32(e##Q##h, vld1q_s32(x_ + 4)); } while (0)
        XP_NAC(0); XP_NAC(1); XP_NAC(2); XP_NAC(3);
#undef XP_NAC
      }
#define XP_NST(Q) do { vst1q_f64(A + 8 * Q, a##Q##0); vst1q_f64(A + 8 * Q + 2, a##Q##1); \
          vst1q_f64(A + 8 * Q + 4, a##Q##2); vst1q_f64(A + 8 * Q + 6, a##Q##3); \
          NADD64(E + 8 * Q, e##Q##l, e##Q##h); } while (0)
      XP_NST(0); XP_NST(1); XP_NST(2); XP_NST(3);
#undef XP_NST
    }
    for (; q < nq; q++) {
      double *A = accD + (size_t)q * XP_EB; int64_t *E = accE + (size_t)q * XP_EB;
      float64x2_t a0 = vld1q_f64(A), a1 = vld1q_f64(A + 2), a2 = vld1q_f64(A + 4), a3 = vld1q_f64(A + 6);
      int32x4_t el = vdupq_n_s32(0), eh = vdupq_n_s32(0);
      const uint8_t *c0 = codes + (size_t)q * nsites;
      for (int s = s0; s < s1; s++) {
        const double *w0 = W + ((size_t)s * NC + c0[s]) * XP_EB;
        const int32_t *x0 = SE + (size_t)s * 2 * XP_EB + (c0[s] != GAP) * XP_EB;
        a0 = vaddq_f64(a0, vld1q_f64(w0)); a1 = vaddq_f64(a1, vld1q_f64(w0 + 2));
        a2 = vaddq_f64(a2, vld1q_f64(w0 + 4)); a3 = vaddq_f64(a3, vld1q_f64(w0 + 6));
        el = vaddq_s32(el, vld1q_s32(x0)); eh = vaddq_s32(eh, vld1q_s32(x0 + 4));
      }
      vst1q_f64(A, a0); vst1q_f64(A + 2, a1); vst1q_f64(A + 4, a2); vst1q_f64(A + 6, a3);
      NADD64(E, el, eh);
    }
  }
  for (int q = 0; q < nq; q++)
    for (int lane = 0; lane < ne; lane++)
      out[(size_t)q * ldo + col0 + lane] =
        accD[(size_t)q * XP_EB + lane] + (M_LN2 * (0. + (double)accE[(size_t)q * XP_EB + lane]));
#endif
  free(W); free(SE); free(accD); free(accE);
  return 0;
}

/* ---------------------------------------------------------- three-taxon tree -- */

/* Memory-lean three-taxon tree (xp3). Same values as the xp.c TT, computed with less traffic:
 *  - the pendant (query) edge: a query site vector is one of nu <= 64 distinct LUT vectors, and
 *    evolve() is per-site, so pend_v[s] == QP[cat[s]][qc[s]] where QP[r][u] = evolve(U[u]) under
 *    P_r(pend_bl) -- identical arithmetic, computed for nu vectors instead of m sites;
 *  - the distal/proximal edges: PD[s][st] = (statd[st]*prox_v[s][st])*dist_v[s][st], which is the
 *    left part of mat_log_like3's statd*x*y*z (x=prox, y=dist, z=pend) with the same association;
 *    dist_v/prox_v themselves are never stored;
 *  - log-likelihood: ll = sum_s log(sum_st PD[s][st]*z[st]) in site order, + log2*etot. */
typedef struct {
  int m;                 /* number of (query-informative) masked sites */
  const int32_t *cat;    /* site category per masked site */
  const double *statd;
  const double *x, *l, *xit, *rates; int nrates;
  uint8_t rmask[64];     /* rates occupied by these masked sites */
  double *tensorD, *tensorP, *tensorQ;   /* nrates x NS x NS (transposed) for dist/prox/pend bl */
  double *dist_o, *prox_o;               /* m x NS originals at the masked sites */
  int64_t *dist_e, *prox_e;
  double etot;                           /* (tot(prox_e) + tot(dist_e)) + tot(pend_e = 0) */
  int nu; double *U; uint8_t *qc;        /* distinct query vectors, per-site index */
  double *QP;                            /* nrates x nu x NS evolved query vectors at pend_bl */
  double *PD;                            /* m x NS */
  double ll_core;                        /* mat_log_like3 value at the current state */
  double prox_bl, dist_bl, pend_bl;
  long n_calls;
  int assert_fail;
  int gmix, R;           /* gamma mixture (every site under all R rates), see xp_gtt_new */
} TT;

/* Gamma mixture (Gmix_model) three-taxon passes. Glv layout here is site-major [m][R][NS]
 * (pplacer's is rate-major [R][m][NS]); only the per-site accumulation order matters and it is
 * pplacer's: util[s] += ((statd*x)*y)*z over rates ascending, then states ascending
 * (ten_log_like3_c), then ll = sum_s log(util[s]) - (float)m * log((float)R). */
static void gtt_ll_pass(TT *t)
{
  double ll_tot = 0;
  const int nu = t->nu, R = t->R;
  for (int s = 0; s < t->m; s++) {
    double util = 0.0;
    for (int r = 0; r < R; r++) {
      const double *pd = t->PD + ((size_t)s * R + r) * NS;
      const double *zz = t->QP + ((size_t)r * nu + t->qc[s]) * NS;
      for (int st = 0; st < NS; st++) util += pd[st] * zz[st];
    }
    ll_tot += log(util);
  }
  ll_tot -= ((float) t->m) * log((float) R);
  t->ll_core = ll_tot;
}

static void gtt_dist_pass(TT *t)
{
  const int R = t->R;
  for (int s = 0; s < t->m; s++)
    for (int r = 0; r < R; r++) {
      double dv[NS], pv[NS];
      const size_t o = ((size_t)s * R + r) * NS, off = (size_t)r * NS * NS;
      evolve_site(dv, t->dist_o + o, t->tensorD + off);
      evolve_site(pv, t->prox_o + o, t->tensorP + off);
      double *pd = t->PD + o;
      for (int st = 0; st < NS; st++) pd[st] = (t->statd[st] * pv[st]) * dv[st];
    }
  gtt_ll_pass(t);
}

void xp_tt_free(TT *t);

static void tt_build_QP(TT *t, double v)
{
  t->pend_bl = v;
  xp_multi_exp(t->tensorQ, t->x, t->l, t->xit, t->rates, t->nrates, t->rmask, v);
  for (int r = 0; r < t->nrates; r++) {
    if (!t->rmask[r]) continue;
    const double *PT = t->tensorQ + (size_t)r * NS * NS;
    double *q = t->QP + (size_t)r * t->nu * NS;
    for (int u = 0; u < t->nu; u++) evolve_site(q + (size_t)u * NS, t->U + (size_t)u * NS, PT);
  }
}

static void tt_ll_pass(TT *t)                   /* ll_core from PD and QP */
{
  if (t->gmix) { gtt_ll_pass(t); return; }
  double ll_tot = 0;
  const int nu = t->nu;
  for (int s = 0; s < t->m; s++) {
    const double *pd = t->PD + (size_t)s * NS;
    const double *zz = t->QP + ((size_t)t->cat[s] * nu + t->qc[s]) * NS;
    double util = 0;
    for (int st = 0; st < NS; st++) util += pd[st] * zz[st];
    ll_tot += log(util);
  }
  t->ll_core = ll_tot;
}

static void tt_dist_pass(TT *t)                 /* PD and ll_core from dist_o/prox_o, QP */
{
  if (t->gmix) { gtt_dist_pass(t); return; }
  double ll_tot = 0;
  const int nu = t->nu;
#if !XP_SIMD20
  for (int s = 0; s < t->m; s++) {
    double dv[NS], pv[NS];
    const size_t off = (size_t)t->cat[s] * NS * NS;
    evolve_site(dv, t->dist_o + (size_t)s * NS, t->tensorD + off);
    evolve_site(pv, t->prox_o + (size_t)s * NS, t->tensorP + off);
    double *pd = t->PD + (size_t)s * NS;
    for (int st = 0; st < NS; st++) pd[st] = (t->statd[st] * pv[st]) * dv[st];
    const double *zz = t->QP + ((size_t)t->cat[s] * nu + t->qc[s]) * NS;
    double util = 0;
    for (int st = 0; st < NS; st++) util += pd[st] * zz[st];
    ll_tot += log(util);
  }
  t->ll_core = ll_tot;
#elif XP_AVX2
  const __m256d s0 = _mm256_loadu_pd(t->statd + 0), s1 = _mm256_loadu_pd(t->statd + 4),
                s2 = _mm256_loadu_pd(t->statd + 8), s3 = _mm256_loadu_pd(t->statd + 12),
                s4 = _mm256_loadu_pd(t->statd + 16);
  for (int s = 0; s < t->m; s++) {
    double dv[NS], pv[NS];
    const size_t off = (size_t)t->cat[s] * NS * NS;
    evolve_site(dv, t->dist_o + (size_t)s * NS, t->tensorD + off);
    evolve_site(pv, t->prox_o + (size_t)s * NS, t->tensorP + off);
    double *pd = t->PD + (size_t)s * NS;
    _mm256_storeu_pd(pd + 0, _mm256_mul_pd(_mm256_mul_pd(s0, _mm256_loadu_pd(pv + 0)), _mm256_loadu_pd(dv + 0)));
    _mm256_storeu_pd(pd + 4, _mm256_mul_pd(_mm256_mul_pd(s1, _mm256_loadu_pd(pv + 4)), _mm256_loadu_pd(dv + 4)));
    _mm256_storeu_pd(pd + 8, _mm256_mul_pd(_mm256_mul_pd(s2, _mm256_loadu_pd(pv + 8)), _mm256_loadu_pd(dv + 8)));
    _mm256_storeu_pd(pd + 12, _mm256_mul_pd(_mm256_mul_pd(s3, _mm256_loadu_pd(pv + 12)), _mm256_loadu_pd(dv + 12)));
    _mm256_storeu_pd(pd + 16, _mm256_mul_pd(_mm256_mul_pd(s4, _mm256_loadu_pd(pv + 16)), _mm256_loadu_pd(dv + 16)));
    const double *zz = t->QP + ((size_t)t->cat[s] * nu + t->qc[s]) * NS;
    double util = 0;
    for (int st = 0; st < NS; st++) util += pd[st] * zz[st];
    ll_tot += log(util);
  }
  t->ll_core = ll_tot;
#else                                              /* NEON: PD = (statd * prox) * dist, 2 states per op */
  const double *sd = t->statd;
  const float64x2_t s0 = vld1q_f64(sd + 0), s1 = vld1q_f64(sd + 2), s2 = vld1q_f64(sd + 4), s3 = vld1q_f64(sd + 6),
                    s4 = vld1q_f64(sd + 8), s5 = vld1q_f64(sd + 10), s6 = vld1q_f64(sd + 12), s7 = vld1q_f64(sd + 14),
                    s8 = vld1q_f64(sd + 16), s9 = vld1q_f64(sd + 18);
  for (int s = 0; s < t->m; s++) {
    double dv[NS], pv[NS];
    const size_t off = (size_t)t->cat[s] * NS * NS;
    evolve_site(dv, t->dist_o + (size_t)s * NS, t->tensorD + off);
    evolve_site(pv, t->prox_o + (size_t)s * NS, t->tensorP + off);
    double *pd = t->PD + (size_t)s * NS;
#define XP_NPD(r) vst1q_f64(pd + 2 * r, vmulq_f64(vmulq_f64(s##r, vld1q_f64(pv + 2 * r)), vld1q_f64(dv + 2 * r)))
    XP_NPD(0); XP_NPD(1); XP_NPD(2); XP_NPD(3); XP_NPD(4); XP_NPD(5); XP_NPD(6); XP_NPD(7); XP_NPD(8); XP_NPD(9);
#undef XP_NPD
    const double *zz = t->QP + ((size_t)t->cat[s] * nu + t->qc[s]) * NS;
    double util = 0;
    for (int st = 0; st < NS; st++) util += pd[st] * zz[st];
    ll_tot += log(util);
  }
  t->ll_core = ll_tot;
#endif
}

static double tt_log_like(TT *t) { return t->ll_core + (M_LN2 * t->etot); }

static void tt_set_pend(TT *t, double v) { tt_build_QP(t, v); tt_ll_pass(t); }

static void tt_set_dist(TT *t, double v)
{
  double cut = t->dist_bl + t->prox_bl;
  double prox = cut - v;
  if (!(prox >= 0.)) t->assert_fail = 1;           /* alpha19: assert(prox_bl >= 0.) */
  t->dist_bl = v;
  xp_multi_exp(t->tensorD, t->x, t->l, t->xit, t->rates, t->nrates, t->rmask, v);
  t->prox_bl = prox;
  xp_multi_exp(t->tensorP, t->x, t->l, t->xit, t->rates, t->nrates, t->rmask, prox);
  tt_dist_pass(t);
}

typedef struct { TT *t; int which; int in_gsl; } OptCtx;  /* which: 0 = pendant, 1 = distal */

#include <stdio.h>
static FILE *xp_tr = NULL;                 /* optional bit-level trace (debugging only) */
void xp_trace_open(const char *path) { xp_tr = fopen(path, "w"); }
void xp_trace_close(void) { if (xp_tr) fclose(xp_tr); xp_tr = NULL; }

static double opt_fun(double v, void *p)
{
  OptCtx *c = (OptCtx *)p;
  c->t->n_calls++;
  if (xp_tr) fprintf(xp_tr, "%c %a\n", c->which == 0 ? 'P' : 'D', v);
  if (c->which == 1 && !((c->t->dist_bl + c->t->prox_bl) - v >= 0.)) {
    /* alpha19 set_dist_bl: assert(prox_bl >= 0.) raises BEFORE any branch length is changed.
     * Inside a GSL callback, ocamlgsl's gslfun_callback uses caml_callback_exn and returns NaN
     * on exception; GSL's SAFE_FUNC_CALL then fails with GSL_EBADFUNC -> Gsl_exn in OCaml, which
     * safe_ml_optimize_location catches (retry from default branch lengths). Outside GSL (the
     * OCaml bisection start finder) the Assert_failure is uncaught and would abort pplacer. */
    if (c->in_gsl) return NAN;
    c->t->assert_fail = 1; return NAN;
  }
  if (c->which == 0) tt_set_pend(c->t, v); else tt_set_dist(c->t, v);
  if (xp_tr) fprintf(xp_tr, "L %a\n", c->t->ll_core);
  return -tt_log_like(c->t);
}

/* status codes */
#define XP_OK 0
#define XP_EXCEEDED 1
#define XP_GSL 2
#define XP_FINDSTART 3
#define XP_INVALIDSTART 4
#define XP_ASSERT 5

/* Minimization.brent with bisection_start_finder (max_iters=100), on GSL's Brent. */
static int xp_brent(OptCtx *c, double raw_start, double left, double right, double tol, double *res)
{
  const int max_iters = 100;
  if (left >= raw_start || raw_start >= right) return XP_INVALIDSTART;
  /* OCaml evaluates the tuple (f left, f right) right-to-left */
  double righty = opt_fun(right, c);
  double lefty = opt_fun(left, c);
  double smaller = (lefty < righty) ? left : right;
  double miny = (lefty <= righty) ? lefty : righty;
  double prev = raw_start, start;
  int it = 1;
  for (;;) {
    double pv = opt_fun(prev, c);
    if (pv < miny) { start = prev; break; }
    else if (fabs(prev - smaller) < tol) { *res = (pv < miny) ? prev : smaller; return XP_OK; }
    else if (it > max_iters) return XP_FINDSTART;
    prev = (prev + smaller) / 2.; it++;
  }
  gsl_function F; F.function = &opt_fun; F.params = c;
  gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
  c->in_gsl = 1;
  int st = gsl_min_fminimizer_set(s, &F, start, left, right);
  if (st) { c->in_gsl = 0; gsl_min_fminimizer_free(s); return XP_GSL; }
  for (int step = 1; ; step++) {
    if (step > max_iters) { c->in_gsl = 0; gsl_min_fminimizer_free(s); return XP_EXCEEDED; }
    st = gsl_min_fminimizer_iterate(s);
    if (st) { c->in_gsl = 0; gsl_min_fminimizer_free(s); return XP_GSL; }
    double lo = gsl_min_fminimizer_x_lower(s), hi = gsl_min_fminimizer_x_upper(s);
    if (hi - lo < tol) { c->in_gsl = 0; *res = gsl_min_fminimizer_x_minimum(s); gsl_min_fminimizer_free(s); return XP_OK; }
  }
}

/* Three_tax.optimize (alternating pendant then distal; stop when both move <= tol). */
int xp_tt_optimize(TT *t, double tol, double max_pend, int max_iter)
{
  OptCtx cp = { t, 0, 0 }, cd = { t, 1, 0 };
  double prev_pend = t->pend_bl, prev_dist = t->dist_bl;
  t->n_calls = 0;
  if (xp_tr) fprintf(xp_tr, "O\n");
  for (int step = 1; ; step++) {
    if (step > max_iter) return XP_EXCEEDED;
    double cur_pend, cur_dist; int st;
    /* optimize_pend_bl: if min = max then start */
    if (1e-8 == max_pend) cur_pend = t->pend_bl;
    else { st = xp_brent(&cp, t->pend_bl, 1e-8, max_pend, tol, &cur_pend); if (st) return st; }
    if (t->assert_fail) return XP_ASSERT;
    double cut = t->dist_bl + t->prox_bl;
    if (0. == cut) cur_dist = t->dist_bl;
    else { st = xp_brent(&cd, t->dist_bl, 0., cut, tol, &cur_dist); if (st) return st; }
    if (t->assert_fail) return XP_ASSERT;
    if (fabs(prev_pend - cur_pend) > tol || fabs(prev_dist - cur_dist) > tol) {
      prev_pend = cur_pend; prev_dist = cur_dist;
    } else return XP_OK;
  }
}

/* Allocate a TT for one query (masked sites m). Buffers owned by TT. Returns NULL if the query
 * has more than 64 distinct site vectors (impossible for protein LUT input). */
TT *xp_tt_new(int m, const int32_t *cat, const double *statd, const double *x, const double *l,
              const double *xit, const double *rates, int nrates, const double *query_masked)
{
  TT *t = (TT *)calloc(1, sizeof(TT));
  t->m = m; t->cat = cat; t->statd = statd; t->x = x; t->l = l; t->xit = xit;
  t->rates = rates; t->nrates = nrates;
  for (int s = 0; s < m; s++) t->rmask[cat[s]] = 1;
  size_t n = (size_t)m * NS;
  t->tensorD = (double *)malloc(sizeof(double) * nrates * NS * NS);
  t->tensorP = (double *)malloc(sizeof(double) * nrates * NS * NS);
  t->tensorQ = (double *)malloc(sizeof(double) * nrates * NS * NS);
  t->dist_o = malloc(sizeof(double) * n); t->prox_o = malloc(sizeof(double) * n);
  t->PD = malloc(sizeof(double) * n);
  t->dist_e = calloc(m, sizeof(int64_t)); t->prox_e = calloc(m, sizeof(int64_t));
  t->U = malloc(sizeof(double) * 64 * NS); t->qc = malloc(m > 0 ? m : 1);
  t->nu = 0;
  for (int s = 0; s < m; s++) {
    const double *q = query_masked + (size_t)s * NS;
    int u = 0;
    while (u < t->nu && memcmp(t->U + (size_t)u * NS, q, sizeof(double) * NS) != 0) u++;
    if (u == t->nu) {
      if (t->nu == 64) { xp_tt_free(t); return NULL; }
      memcpy(t->U + (size_t)u * NS, q, sizeof(double) * NS); t->nu++;
    }
    t->qc[s] = (uint8_t)u;
  }
  t->QP = malloc(sizeof(double) * (size_t)nrates * (t->nu > 0 ? t->nu : 1) * NS);
  return t;
}

void xp_tt_free(TT *t)
{
  free(t->tensorD); free(t->tensorP); free(t->tensorQ); free(t->dist_o); free(t->prox_o);
  free(t->PD); free(t->dist_e); free(t->prox_e); free(t->U); free(t->qc); free(t->QP); free(t);
}

/* set_tt_edges: copy masked sites of darr[loc]/parr[loc] (full-site glvs) into the edges,
 * then set pendant, distal, proximal branch lengths (in pplacer's order). */
void xp_tt_set(TT *t, const double *darr_a, const int64_t *darr_e, const double *parr_a,
               const int64_t *parr_e, const int32_t *reind, double cut_bl,
               double pendant, double distal)
{
  for (int i = 0; i < t->m; i++) {
    memcpy(t->dist_o + (size_t)i * NS, darr_a + (size_t)reind[i] * NS, sizeof(double) * NS);
    memcpy(t->prox_o + (size_t)i * NS, parr_a + (size_t)reind[i] * NS, sizeof(double) * NS);
    t->dist_e[i] = darr_e[reind[i]];
    t->prox_e[i] = parr_e[reind[i]];
  }
  /* tt_log_like in xp.c: (tot(prox_e) + tot(dist_e)) + tot(pend_e), pend_e == 0 */
  t->etot = (xp_total_twoexp(t->prox_e, t->m) + xp_total_twoexp(t->dist_e, t->m)) + 0.0;
  /* pendant first (Glv_edge.set_bl pend_edge), then dist edge, then prox edge */
  tt_build_QP(t, pendant);
  t->dist_bl = distal;
  xp_multi_exp(t->tensorD, t->x, t->l, t->xit, t->rates, t->nrates, t->rmask, distal);
  t->prox_bl = cut_bl - distal;
  xp_multi_exp(t->tensorP, t->x, t->l, t->xit, t->rates, t->nrates, t->rmask, cut_bl - distal);
  tt_dist_pass(t);
  t->assert_fail = 0;
}

void xp_tt_results(TT *t, double *out)   /* (log_like, pend_bl, dist_bl, n_calls) */
{
  out[0] = tt_log_like(t); out[1] = t->pend_bl; out[2] = t->dist_bl; out[3] = (double)t->n_calls;
}

/* ======================================================================================
 * Gamma mixture model (pplacer gmix_model.ml, gamma.ml, linear_c.c "tensor" functions).
 * Glvs are stored site-major: a[s][r][st] (pplacer: rate-major a[r][s][st]); one exponent per
 * site shared by all rates. Every value below uses pplacer's operands and per-value order.
 * ====================================================================================== */

/* Gamma.discrete_gamma: inverse of the gamma CDF, below a tiny epsilon -> 0. */
static double inverse_gamma_cdf(double alpha, double beta, double y)
{
  const double epsilon = 1e-7;
  beta = 1. / beta;
  if (y < gsl_cdf_gamma_P(epsilon, alpha, beta)) return 0.;
  return gsl_cdf_gamma_Pinv(y, alpha, beta);
}

/* Gamma.discrete_gamma n_cat alpha (Yang 1994, mean rates), then make_mean_one. Returns 0, -1 for
 * n_cat <= 0, or 2 if GSL raised an error (where pplacer stops; xp_gamma_gsl_error tells which). */
int xp_discrete_gamma(int n_cat, double alpha, double *out)
{
  if (n_cat <= 0) return -1;
  double *freq = (double *)malloc(sizeof(double) * (n_cat + 1));
  xp_gsl_errno = 0;
  xp_gsl_reason[0] = 0;
  xp_gsl_recording = 1;
  for (int i = 0; i <= n_cat; i++) {
    if (i == 0) freq[i] = 0.;
    else if (i == n_cat) freq[i] = 1.;
    else {
      double cutpoint = inverse_gamma_cdf(alpha, 1., (double)i / (double)n_cat);
      freq[i] = gsl_sf_gamma_inc_P(alpha + 1., cutpoint);
    }
  }
  xp_gsl_recording = 0;
  if (xp_gsl_errno) {
    free(freq);
    return 2;
  }
  double tot = 0.;
  for (int i = 0; i < n_cat; i++) { out[i] = freq[i + 1] - freq[i]; tot = tot + out[i]; }
  const double f_n = (double)n_cat;
  for (int i = 0; i < n_cat; i++) out[i] = f_n * out[i] / tot;
  free(freq);
  return 0;
}

/* ten_pairwise_prod + iba1_pairwise_sum */
void xp_g_pairwise(double *dst, int64_t *de, const double *x, const int64_t *xe,
                   const double *y, const int64_t *ye, int nsites, int R)
{
  size_t n = (size_t)nsites * R * NS;
  for (size_t i = 0; i < n; i++) dst[i] = x[i] * y[i];
  for (int s = 0; s < nsites; s++) de[s] = xe[s] + ye[s];
}

/* ten_statd_pairwise_prod: statd[st] * a * b */
void xp_g_statd_pairwise(double *dst, int64_t *de, const double *statd, const double *a,
                         const int64_t *ae, const double *b, const int64_t *be, int nsites, int R)
{
  size_t n = (size_t)nsites * R;
  for (size_t k = 0; k < n; k++)
    for (int st = 0; st < NS; st++) {
      *dst = statd[st] * (*a) * (*b);
      dst++; a++; b++;
    }
  for (int s = 0; s < nsites; s++) de[s] = ae[s] + be[s];
}

/* Gmix Glv.perhaps_pull_exponent: per site, the largest frexp exponent over the rates' maxima;
 * if below min_twoexp, every rate's slice is scaled by of_twoexp(-max) = ldexp(0.5, -max+1). */
/* Optional fix (NOT pplacer behaviour; gplace --fix_gamma_underflow): frexp(0.0) gives exponent 0,
 * so in pplacer a rate whose values have all underflowed to 0.0 blocks rescaling of the whole site;
 * on large trees the other rates then underflow too (site likelihood 0 -> -inf). With xp_gamma_fix
 * set, all-zero rates are skipped when the exponent is chosen. The setting is per process (forked
 * workers inherit it). */
static int xp_gamma_fix = 0;
void xp_set_gamma_fix(int on) { xp_gamma_fix = on; }

void xp_g_pull_exponent(double *a, int64_t *e, int nsites, int R, int min_twoexp)
{
  for (int s = 0; s < nsites; s++) {
    double *site = a + (size_t)s * R * NS;
    int max_twoexp = -INT_MAX;
    for (int r = 0; r < R; r++) {
      const double *v = site + (size_t)r * NS;
      double mx = v[0];                          /* gsl_vector_max */
      for (int i = 1; i < NS; i++) if (v[i] > mx) mx = v[i];
      if (xp_gamma_fix && mx == 0.0) continue;
      int twoexp; frexp(mx, &twoexp);
      if (twoexp > max_twoexp) max_twoexp = twoexp;
    }
    if (max_twoexp == -INT_MAX) continue;        /* only with the fix: the whole site is zero */
    if (max_twoexp < min_twoexp) {
      double f = ldexp(0.5, -max_twoexp + 1);
      for (int i = 0; i < R * NS; i++) site[i] *= f;
      e[s] += max_twoexp;
    }
  }
}

/* Gmix prefilter (Glv.masked_logdot of the evolved query and the supernode), 8 edges at once:
 *   W[s][c] = log( util ), util = sum over r ascending, st ascending of qev[r][c][st]*sn[s][r][st]
 *   ll = +0.0 + sum_{s informative} W[s][code]           (ten_masked_logdot_c ...)
 *   ll = ll - (float)nused[q] * log((float)R)            (... its final subtraction)
 *   score = ll + log2*(0 + sum_{s informative} se[s])    (masked_total_twoexp, query e = 0)
 * sn8: 8 x nsites x R x NS; se8: 8 x nsites; qevtab: R x ncode x NS; codes as xp_edge_scores8. */
int xp_g_edge_scores8(double *out, long ldo, long col0, int ne, const double *sn8, const int64_t *se8,
                      const double *qevtab, int R, int nsites, int ncode, const uint8_t *codes,
                      const int32_t *nused, int nq)
{
  if (ne < 1 || ne > XP_EB || ncode < 1 || ncode > 254) return -1;
  const int NC = ncode + 1;
  double *W = (double *)malloc(sizeof(double) * (size_t)nsites * NC * XP_EB);
  if (!W) return -1;
  for (int s = 0; s < nsites; s++) {
    double *Ws = W + (size_t)s * NC * XP_EB;
    for (int c = 0; c < ncode; c++)
      for (int lane = 0; lane < XP_EB; lane++) {
        if (lane < ne) {
          const double *sn = sn8 + ((size_t)lane * nsites + s) * R * NS;
          double util = 0.0;
          for (int r = 0; r < R; r++) {
            const double *qv = qevtab + ((size_t)r * ncode + c) * NS;
            const double *sv = sn + (size_t)r * NS;
            for (int st = 0; st < NS; st++) util += qv[st] * sv[st];
          }
          Ws[c * XP_EB + lane] = log(util);
        } else Ws[c * XP_EB + lane] = 0.0;
      }
  }
  const double logR = log((float) R);
  for (int q = 0; q < nq; q++) {
    const uint8_t *cq = codes + (size_t)q * nsites;
    for (int lane = 0; lane < ne; lane++) {
      double ll = 0; long long et = 0;
      for (int s = 0; s < nsites; s++) {
        unsigned c = cq[s];
        if (c != (unsigned)ncode) { ll += W[((size_t)s * NC + c) * XP_EB + lane]; et += se8[(size_t)lane * nsites + s]; }
      }
      ll -= ((float) nused[q]) * logR;
      out[(size_t)q * ldo + col0 + lane] = ll + (M_LN2 * (0. + (double)et));
    }
  }
  free(W);
  return 0;
}

/* Gmix three-taxon tree: m masked sites, all R rates at every site (rmask all set). */
TT *xp_gtt_new(int m, int R, const double *statd, const double *x, const double *l,
               const double *xit, const double *rates, const double *query_masked)
{
  if (R > 64) return NULL;
  static int32_t zero_cat[1] = {0};
  TT *t = xp_tt_new(0, zero_cat, statd, x, l, xit, rates, R, query_masked);   /* tensors */
  free(t->dist_o); free(t->prox_o); free(t->PD); free(t->dist_e); free(t->prox_e);
  free(t->U); free(t->qc); free(t->QP);
  t->gmix = 1; t->R = R; t->m = m; t->cat = NULL;
  for (int r = 0; r < R; r++) t->rmask[r] = 1;
  size_t n = (size_t)m * R * NS;
  t->dist_o = malloc(sizeof(double) * (n ? n : 1)); t->prox_o = malloc(sizeof(double) * (n ? n : 1));
  t->PD = malloc(sizeof(double) * (n ? n : 1));
  t->dist_e = calloc(m ? m : 1, sizeof(int64_t)); t->prox_e = calloc(m ? m : 1, sizeof(int64_t));
  t->U = malloc(sizeof(double) * 64 * NS); t->qc = malloc(m > 0 ? m : 1);
  t->nu = 0;
  for (int s = 0; s < m; s++) {            /* the query glv is constant across rates */
    const double *q = query_masked + (size_t)s * NS;
    int u = 0;
    while (u < t->nu && memcmp(t->U + (size_t)u * NS, q, sizeof(double) * NS) != 0) u++;
    if (u == t->nu) {
      if (t->nu == 64) { xp_tt_free(t); return NULL; }
      memcpy(t->U + (size_t)u * NS, q, sizeof(double) * NS); t->nu++;
    }
    t->qc[s] = (uint8_t)u;
  }
  t->QP = malloc(sizeof(double) * (size_t)R * (t->nu > 0 ? t->nu : 1) * NS);
  return t;
}

/* xp_tt_set for gmix: full-site glvs darr/parr are nsites x R x NS */
void xp_gtt_set(TT *t, const double *darr_a, const int64_t *darr_e, const double *parr_a,
                const int64_t *parr_e, const int32_t *reind, double cut_bl, double pendant, double distal)
{
  const size_t w = (size_t)t->R * NS;
  for (int i = 0; i < t->m; i++) {
    memcpy(t->dist_o + (size_t)i * w, darr_a + (size_t)reind[i] * w, sizeof(double) * w);
    memcpy(t->prox_o + (size_t)i * w, parr_a + (size_t)reind[i] * w, sizeof(double) * w);
    t->dist_e[i] = darr_e[reind[i]];
    t->prox_e[i] = parr_e[reind[i]];
  }
  t->etot = (xp_total_twoexp(t->prox_e, t->m) + xp_total_twoexp(t->dist_e, t->m)) + 0.0;
  tt_build_QP(t, pendant);
  t->dist_bl = distal;
  xp_multi_exp(t->tensorD, t->x, t->l, t->xit, t->rates, t->nrates, t->rmask, distal);
  t->prox_bl = cut_bl - distal;
  xp_multi_exp(t->tensorP, t->x, t->l, t->xit, t->rates, t->nrates, t->rmask, cut_bl - distal);
  tt_dist_pass(t);
  t->assert_fail = 0;
}
