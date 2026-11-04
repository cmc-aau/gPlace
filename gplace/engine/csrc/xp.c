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
#include <pthread.h>
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

/* Opt-in exact speedups, read from the environment once (default off). */
static int xp_switches_read = 0, xp_vlog = 0, xp_vdot = 0;

/* Site-block lengths of the prefilter kernels (xp_set_pf_bs; GPLACE_PF_BS): the sites whose table
 * rows are built before the queries are run over them. Scheduling only: every query's per-edge
 * values are added (f64) or multiplied (f32) in site order whatever the block length, so the scores
 * do not depend on it. Within a block the f64 kernel's int32 exponent sums are flushed into the
 * int64 totals every XP_EBLK sites (|se| <= 2^23, so they cannot overflow) and the f32 products are
 * renormalised every XP_RBLK sites (64 mantissas in [1, 2) stay below 2^64; renormalising is an
 * exact power-of-two scaling). NEON's f64 loop flushes once per block: its blocks stay <= 255. */
#define XP_EBLK 128
#define XP_RBLK 64
#define XP_PF_BS_MAX 4096
static int xp_pf_bs64 = 128, xp_pf_bs32 = 64;
void xp_set_pf_bs(int bs64, int bs32)
{
  if (bs64 >= 1 && bs64 <= XP_PF_BS_MAX) xp_pf_bs64 = bs64;
  if (bs32 >= 1 && bs32 <= XP_PF_BS_MAX) xp_pf_bs32 = bs32;
}
static void xp_read_switches(void)
{
  if (xp_switches_read) return;
  const char *v = getenv("GPLACE_VLOG"), *d = getenv("GPLACE_VDOT");
  xp_vlog = (v && *v == '1') && (XP_AVX2 || XP_NEON);
  xp_vdot = (d && *d == '1') && NS == 20;
  xp_switches_read = 1;
}

#if XP_AVX2 || XP_NEON
#include "xp_vlog_tab.h"
double cr_log(double);
#endif
#if XP_AVX2
/* Four-lane correctly rounded log (GPLACE_VLOG=1): CORE-MATH's cr_log fast path with
 * its proven absolute bound 0x1.b6p-69 and rounding test, as transcribed operation for operation
 * by crmvec's log_fast (MIT, https://github.com/anun333/crmvec, crmvec.c at df66af3), gather form.
 * A lane the test cannot decide, and any x that is not a positive normal finite number other than
 * 1, goes to cr_log itself -- the function gPlace's log is (libm_bind.c) -- so every lane returns
 * cr_log's correctly rounded result. */
static inline __m256d xp_log4(__m256d x)
{
  const __m256i MANT = _mm256_set1_epi64x(0xfffffffffffffLL);
  __m256i u = _mm256_castpd_si256(x);
  __m256i ok = _mm256_and_si256(_mm256_cmpgt_epi64(u, _mm256_set1_epi64x(0x000fffffffffffffLL)),
                                _mm256_cmpgt_epi64(_mm256_set1_epi64x(0x7ff0000000000000LL), u));
  u = _mm256_blendv_epi8(_mm256_castpd_si256(_mm256_set1_pd(1.0)), u, ok);
  __m256i m = _mm256_or_si256(_mm256_and_si256(u, MANT), _mm256_set1_epi64x(1LL << 52));
  __m256i c = _mm256_cmpgt_epi64(m, _mm256_set1_epi64x(0x16a09e667f3bcdLL - 1));
  __m256i idx = _mm256_sub_epi64(_mm256_srlv_epi64(m, _mm256_sub_epi64(_mm256_set1_epi64x(43), c)),
                                 _mm256_set1_epi64x(362));
  __m256d vf = _mm256_castsi256_pd(_mm256_or_si256(_mm256_and_si256(u, MANT), _mm256_set1_epi64x(0x3ff0000000000000LL)));
  __m256d y = _mm256_blendv_pd(vf, _mm256_mul_pd(vf, _mm256_set1_pd(0.5)), _mm256_castsi256_pd(c));
  __m256i e = _mm256_sub_epi64(_mm256_sub_epi64(_mm256_srli_epi64(u, 52), _mm256_set1_epi64x(0x3ff)), c);
  __m256d ee = _mm256_sub_pd(_mm256_castsi256_pd(_mm256_add_epi64(e, _mm256_castpd_si256(_mm256_set1_pd(0x1.8p52)))),
                             _mm256_set1_pd(0x1.8p52));
  __m256d r = _mm256_i64gather_pd(LOG_INVERSE, idx, 8);
  __m256i i2 = _mm256_slli_epi64(idx, 1);
  __m256d l1 = _mm256_i64gather_pd(&LOG_INV[0][0], i2, 8), l2 = _mm256_i64gather_pd(&LOG_INV[0][0] + 1, i2, 8);
  __m256d z = _mm256_fmadd_pd(r, y, _mm256_set1_pd(-1.0));
  __m256d z2 = _mm256_mul_pd(z, z);
  __m256d p45 = _mm256_fmadd_pd(_mm256_set1_pd(-0x1.55362255e0f63p-3), z, _mm256_set1_pd(0x1.999a14758b084p-3));
  __m256d p23 = _mm256_fmadd_pd(_mm256_set1_pd(-0x1.0000000537df6p-2), z, _mm256_set1_pd(0x1.555555554f4d8p-2));
  __m256d ph = _mm256_fmadd_pd(p45, z2, p23);
  ph = _mm256_fmadd_pd(ph, z, _mm256_set1_pd(-0x1.ffffffffffffap-2));
  ph = _mm256_mul_pd(ph, z2);
  __m256d a = _mm256_fmadd_pd(ee, _mm256_set1_pd(0x1.62e42fefa38p-1), l1);
  __m256d h = _mm256_add_pd(a, z);
  __m256d l = _mm256_sub_pd(z, _mm256_sub_pd(h, a));
  l = _mm256_add_pd(ph, _mm256_add_pd(l, l2));
  l = _mm256_fmadd_pd(ee, _mm256_set1_pd(0x1.ef35793c7673p-45), l);
  const __m256d ERR = _mm256_set1_pd(0x1.b6p-69);
  __m256d left = _mm256_add_pd(h, _mm256_sub_pd(l, ERR)), right = _mm256_add_pd(h, _mm256_add_pd(l, ERR));
  __m256d bad = _mm256_or_pd(_mm256_cmp_pd(left, right, _CMP_NEQ_UQ),
                             _mm256_castsi256_pd(_mm256_xor_si256(ok, _mm256_set1_epi64x(-1))));
  int mk = _mm256_movemask_pd(bad);
  if (!mk) return left;
  double xs[4], ys[4]; _mm256_storeu_pd(xs, x); _mm256_storeu_pd(ys, left);
  for (int i = 0; i < 4; i++) if (mk >> i & 1) ys[i] = cr_log(xs[i]);
  return _mm256_loadu_pd(ys);
}

/* Four sites' utils at once (GPLACE_VDOT=1, with GPLACE_VLOG=1), bit for bit the scalar loop
 * `util = 0; for st: util += pd[st] * zz[st]`: each site's products are formed along its row
 * (the same products), five 4x4 transposes put one site in each lane, and each lane then adds
 * its products from 0 in state order -- the scalar sequence exactly. NS == 20. */
static inline __m256d xp_dot4(const double *const p[4], const double *const z[4])
{
  __m256d acc = _mm256_setzero_pd();
  for (int b = 0; b < 20; b += 4) {
    __m256d r0 = _mm256_mul_pd(_mm256_loadu_pd(p[0] + b), _mm256_loadu_pd(z[0] + b));
    __m256d r1 = _mm256_mul_pd(_mm256_loadu_pd(p[1] + b), _mm256_loadu_pd(z[1] + b));
    __m256d r2 = _mm256_mul_pd(_mm256_loadu_pd(p[2] + b), _mm256_loadu_pd(z[2] + b));
    __m256d r3 = _mm256_mul_pd(_mm256_loadu_pd(p[3] + b), _mm256_loadu_pd(z[3] + b));
    __m256d t0 = _mm256_unpacklo_pd(r0, r1), t1 = _mm256_unpackhi_pd(r0, r1);
    __m256d t2 = _mm256_unpacklo_pd(r2, r3), t3 = _mm256_unpackhi_pd(r2, r3);
    acc = _mm256_add_pd(acc, _mm256_permute2f128_pd(t0, t2, 0x20));   /* state b */
    acc = _mm256_add_pd(acc, _mm256_permute2f128_pd(t1, t3, 0x20));   /* b + 1 */
    acc = _mm256_add_pd(acc, _mm256_permute2f128_pd(t0, t2, 0x31));   /* b + 2 */
    acc = _mm256_add_pd(acc, _mm256_permute2f128_pd(t1, t3, 0x31));   /* b + 3 */
  }
  return acc;
}
#endif
#if XP_NEON
/* xp_log4 for NEON, two lanes: the same operations in the same order (vfmaq_f64 is the same fused
 * multiply-add as _mm256_fmadd_pd), the table entries read lane by lane. */
static inline float64x2_t xp_log2(float64x2_t x)
{
  const int64x2_t MANT = vdupq_n_s64(0xfffffffffffffLL);
  int64x2_t u = vreinterpretq_s64_f64(x);
  const uint64x2_t ok = vandq_u64(vcgtq_s64(u, vdupq_n_s64(0x000fffffffffffffLL)),
                                  vcgtq_s64(vdupq_n_s64(0x7ff0000000000000LL), u));
  u = vbslq_s64(ok, u, vreinterpretq_s64_f64(vdupq_n_f64(1.0)));
  const int64x2_t m = vorrq_s64(vandq_s64(u, MANT), vdupq_n_s64(1LL << 52));
  const uint64x2_t c = vcgtq_s64(m, vdupq_n_s64(0x16a09e667f3bcdLL - 1));
  const int64x2_t ci = vreinterpretq_s64_u64(c);                       /* -1 or 0, as AVX2's */
  const int64x2_t idx = vsubq_s64(vreinterpretq_s64_u64(vshlq_u64(vreinterpretq_u64_s64(m),
                                                                  vnegq_s64(vsubq_s64(vdupq_n_s64(43), ci)))),
                                  vdupq_n_s64(362));
  const float64x2_t vf = vreinterpretq_f64_s64(vorrq_s64(vandq_s64(u, MANT), vdupq_n_s64(0x3ff0000000000000LL)));
  const float64x2_t y = vbslq_f64(c, vmulq_f64(vf, vdupq_n_f64(0.5)), vf);
  const int64x2_t e = vsubq_s64(vsubq_s64(vreinterpretq_s64_u64(vshrq_n_u64(vreinterpretq_u64_s64(u), 52)),
                                          vdupq_n_s64(0x3ff)), ci);
  const float64x2_t ee = vsubq_f64(vreinterpretq_f64_s64(vaddq_s64(e, vreinterpretq_s64_f64(vdupq_n_f64(0x1.8p52)))),
                                   vdupq_n_f64(0x1.8p52));
  const int64_t i0 = vgetq_lane_s64(idx, 0), i1 = vgetq_lane_s64(idx, 1);
  const double rr[2] = { LOG_INVERSE[i0], LOG_INVERSE[i1] }, a1[2] = { LOG_INV[i0][0], LOG_INV[i1][0] },
               a2[2] = { LOG_INV[i0][1], LOG_INV[i1][1] };
  const float64x2_t r = vld1q_f64(rr), l1 = vld1q_f64(a1), l2 = vld1q_f64(a2);
  const float64x2_t z = vfmaq_f64(vdupq_n_f64(-1.0), r, y);
  const float64x2_t z2 = vmulq_f64(z, z);
  const float64x2_t p45 = vfmaq_f64(vdupq_n_f64(0x1.999a14758b084p-3), vdupq_n_f64(-0x1.55362255e0f63p-3), z);
  const float64x2_t p23 = vfmaq_f64(vdupq_n_f64(0x1.555555554f4d8p-2), vdupq_n_f64(-0x1.0000000537df6p-2), z);
  float64x2_t ph = vfmaq_f64(p23, p45, z2);
  ph = vfmaq_f64(vdupq_n_f64(-0x1.ffffffffffffap-2), ph, z);
  ph = vmulq_f64(ph, z2);
  const float64x2_t a = vfmaq_f64(l1, ee, vdupq_n_f64(0x1.62e42fefa38p-1));
  const float64x2_t h = vaddq_f64(a, z);
  float64x2_t l = vsubq_f64(z, vsubq_f64(h, a));
  l = vaddq_f64(ph, vaddq_f64(l, l2));
  l = vfmaq_f64(l, ee, vdupq_n_f64(0x1.ef35793c7673p-45));
  const float64x2_t ERR = vdupq_n_f64(0x1.b6p-69);
  const float64x2_t left = vaddq_f64(h, vsubq_f64(l, ERR)), right = vaddq_f64(h, vaddq_f64(l, ERR));
  const uint64x2_t good = vandq_u64(vceqq_f64(left, right), ok);      /* not (unordered or unequal), and ok */
  if (vgetq_lane_u64(good, 0) && vgetq_lane_u64(good, 1)) return left;
  double xs[2], ys[2]; vst1q_f64(xs, x); vst1q_f64(ys, left);
  if (!vgetq_lane_u64(good, 0)) ys[0] = cr_log(xs[0]);
  if (!vgetq_lane_u64(good, 1)) ys[1] = cr_log(xs[1]);
  return vld1q_f64(ys);
}

/* xp_dot4 for NEON, two sites: products along each row, a 2x2 transpose, and each lane adds its
 * products from 0 in state order -- the scalar sequence exactly. NS == 20. */
static inline float64x2_t xp_dot2(const double *p0, const double *z0, const double *p1, const double *z1)
{
  float64x2_t acc = vdupq_n_f64(0.0);
  for (int b = 0; b < 20; b += 2) {
    const float64x2_t r0 = vmulq_f64(vld1q_f64(p0 + b), vld1q_f64(z0 + b));
    const float64x2_t r1 = vmulq_f64(vld1q_f64(p1 + b), vld1q_f64(z1 + b));
    acc = vaddq_f64(acc, vzip1q_f64(r0, r1));                      /* state b */
    acc = vaddq_f64(acc, vzip2q_f64(r0, r1));                      /* b + 1 */
  }
  return acc;
}
#endif
/* Four logs / four sites' dots into memory, on either SIMD engine (GPLACE_VLOG, GPLACE_VDOT). */
#if XP_AVX2
static inline void xp_log4v(const double *u, double *l) { _mm256_storeu_pd(l, xp_log4(_mm256_loadu_pd(u))); }
static inline void xp_dot4v(const double *const p[4], const double *const z[4], double *u)
{ _mm256_storeu_pd(u, xp_dot4(p, z)); }
#elif XP_NEON
static inline void xp_log4v(const double *u, double *l)
{ vst1q_f64(l, xp_log2(vld1q_f64(u))); vst1q_f64(l + 2, xp_log2(vld1q_f64(u + 2))); }
static inline void xp_dot4v(const double *const p[4], const double *const z[4], double *u)
{ vst1q_f64(u, xp_dot2(p[0], z[0], p[1], z[1])); vst1q_f64(u + 2, xp_dot2(p[2], z[2], p[3], z[3])); }
#else                                              /* portable: xp_vlog is 0, the scalar forms */
static inline void xp_log4v(const double *u, double *l) { for (int k = 0; k < 4; k++) l[k] = log(u[k]); }
static inline void xp_dot4v(const double *const p[4], const double *const z[4], double *u)
{ for (int k = 0; k < 4; k++) { double v = 0; for (int st = 0; st < NS; st++) v += p[k][st] * z[k][st]; u[k] = v; } }
#endif
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
/* --full_model (NOT pplacer behaviour): a category of rate exactly 0 (invariable sites, +I) gets the exact
 * identity, as IQ-TREE and RAxML-NG compute it. X exp(0) X^-1 is the identity only up to rounding (~1e-16
 * off the diagonal), and at a variable site, where the identity's likelihood is exactly 0, that noise is
 * added to the other categories' likelihoods, which can be far smaller. pplacer models never set it (their
 * gamma rates can underflow to 0.0, and pplacer then evolves with X X^-1). Per process, as xp_gamma_fix. */
static int xp_zero_rate_identity = 0;
void xp_set_zero_rate_identity(int on) { xp_zero_rate_identity = on; }

static int zero_rate_identity(double *tensor_r, double rate)
{
  if (!xp_zero_rate_identity || rate != 0.0) return 0;
  for (int i = 0; i < NS; i++) tensor_r[i * NS + i] = 1.0;     /* off the diagonal: the caller's zeros */
  return 1;
}

void xp_multi_exp_ref(double *tensor, const double *x, const double *l, const double *xit,
                  const double *rates, int nrates, const uint8_t *mask, double bl)
{
  double util[NS];
  memset(tensor, 0, sizeof(double) * nrates * NS * NS);
  for (int r = 0; r < nrates; r++) {
    if (!mask[r]) continue;
    if (zero_rate_identity(tensor + (size_t)r * NS * NS, rates[r])) continue;
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
    if (zero_rate_identity(tensor + (size_t)r * NS * NS, rates[r])) continue;
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
  xp_read_switches();
  for (int s = 0; s < nsites; s++) {
    const double *qt = qevtab + (size_t)cat[s] * ncode * NS;
    double *Ws = W + (size_t)s * NC * XP_EB;
    for (int c = 0; c < ncode; c++) {
      const double *qv = qt + (size_t)c * NS;
#if XP_AVX2 || XP_NEON
      if (xp_vlog && ne == XP_EB) {
        /* Eight lanes' utils, then two four-lane logs: the same values in the same places. */
        double u8[XP_EB];
        for (int lane = 0; lane < XP_EB; lane++) {
          const double *sn = sn8 + ((size_t)lane * nsites + s) * NS;
          double util = 0;
          for (int st = 0; st < NS; st++) util += qv[st] * sn[st];
          u8[lane] = util;
        }
        xp_log4v(u8, Ws + c * XP_EB);
        xp_log4v(u8 + 4, Ws + c * XP_EB + 4);
        continue;
      }
#endif
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
  const int BS = xp_pf_bs64;
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
      const uint8_t *c0 = codes + (size_t)q * nsites, *c1 = c0 + nsites, *c2 = c1 + nsites, *c3 = c2 + nsites;
      for (int sb = s0; sb < s1; sb += XP_EBLK) {
      const int sb1 = sb + XP_EBLK < s1 ? sb + XP_EBLK : s1;
      __m256i e0 = _mm256_setzero_si256(), e1 = e0, e2 = e0, e3 = e0;   /* partial sums, <= XP_EBLK sites */
      for (int s = sb; s < sb1; s++) {
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
      ADD64(E0, e0); ADD64(E1, e1); ADD64(E2, e2); ADD64(E3, e3);
      }
      _mm256_storeu_pd(A0, a0l); _mm256_storeu_pd(A0 + 4, a0h); _mm256_storeu_pd(A1, a1l); _mm256_storeu_pd(A1 + 4, a1h);
      _mm256_storeu_pd(A2, a2l); _mm256_storeu_pd(A2 + 4, a2h); _mm256_storeu_pd(A3, a3l); _mm256_storeu_pd(A3 + 4, a3h);
    }
    for (; q < nq; q++) {
      double *A0 = accD + (size_t)q * XP_EB; int64_t *E0 = accE + (size_t)q * XP_EB;
      __m256d a0l = _mm256_loadu_pd(A0), a0h = _mm256_loadu_pd(A0 + 4);
      const uint8_t *c0 = codes + (size_t)q * nsites;
      for (int sb = s0; sb < s1; sb += XP_EBLK) {
        const int sb1 = sb + XP_EBLK < s1 ? sb + XP_EBLK : s1;
        __m256i e0 = _mm256_setzero_si256();
        for (int s = sb; s < sb1; s++) {
          const double *w0 = W + ((size_t)s * NC + c0[s]) * XP_EB;
          const int32_t *SEs = SE + (size_t)s * 2 * XP_EB;
          a0l = _mm256_add_pd(a0l, _mm256_loadu_pd(w0)); a0h = _mm256_add_pd(a0h, _mm256_loadu_pd(w0 + 4));
          e0 = _mm256_add_epi32(e0, _mm256_loadu_si256((const __m256i *)(SEs + (c0[s] != GAP) * XP_EB)));
        }
        ADD64(E0, e0);
      }
      _mm256_storeu_pd(A0, a0l); _mm256_storeu_pd(A0 + 4, a0h);
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
  const int BS = xp_pf_bs64 < 255 ? xp_pf_bs64 : 255;   /* one int32 flush per block */
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
  int mix;               /* mixture model (--full_model): rate category r has its own eigen-system and */
  const double *mx, *ml, *mxit, *mstatd;   /* frequencies: R x NS x NS, R x NS, R x NS x NS, R x NS */
  /* PD cache (GPLACE_EXPERIMENTAL_PD_CACHE=1; never part of a set of switches: it was slower than
   * no cache at 64 workers on an AMD EPYC 7763 node): PD depends only on the location's edge vectors and the
   * (distal, proximal) lengths, not on the query's pendant state, and the alternating Brent
   * search re-brackets at the same lengths every round. npd buffers (GPLACE_PD_SLOTS, 1 to
   * PD_SLOTS, default PD_SLOTS), keyed by the exact (dist, prox) pair, cleared by xp_tt_set;
   * t->PD always points at one of them. Fewer slots: less memory per worker, fewer hits. */
  int pdc, npd;
  double *pd_buf[4]; double pd_kd[4], pd_kp[4]; int pd_ok[4]; long pd_age[4], pd_clock;
  /* Distal grouping (GPLACE_DGRP=1): sites whose (distal row, category) are bitwise equal
   * evolve to the same dv, so each distinct pair is evolved once per distal pass and the sites
   * read its result. Classes are rebuilt by xp_tt_set from the new edge vectors. */
  int dgrp, nd;
  int32_t *dcls, *drep, *dhash;          /* class per site, first site of each class, table */
  double *DV;                            /* nd x NS evolved distal rows */
} TT;
#define PD_SLOTS 4

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
      const double *sd = t->mix ? t->mstatd + (size_t)r * NS : t->statd;
      evolve_site(dv, t->dist_o + o, t->tensorD + off);
      evolve_site(pv, t->prox_o + o, t->tensorP + off);
      double *pd = t->PD + o;
      for (int st = 0; st < NS; st++) pd[st] = (sd[st] * pv[st]) * dv[st];
    }
  gtt_ll_pass(t);
}

void xp_tt_free(TT *t);

/* The transition tensor of every occupied rate category at branch length v: from the one
 * eigen-system (pplacer), or, for a mixture, category by category from its own. */
static void tt_tensor(TT *t, double *tensor, double v)
{
  if (!t->mix) { xp_multi_exp(tensor, t->x, t->l, t->xit, t->rates, t->nrates, t->rmask, v); return; }
  for (int r = 0; r < t->nrates; r++)
    xp_multi_exp(tensor + (size_t)r * NS * NS, t->mx + (size_t)r * NS * NS, t->ml + (size_t)r * NS,
                 t->mxit + (size_t)r * NS * NS, t->rates + r, 1, t->rmask + r, v);
}

static void tt_build_QP(TT *t, double v)
{
  t->pend_bl = v;
  tt_tensor(t, t->tensorQ, v);
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
  int s0 = 0;
#if XP_AVX2 || XP_NEON
  if (xp_vlog) {
    /* Four sites' utils, one four-lane log, the four logs added in site order: the same sum. */
    for (; s0 + 4 <= t->m; s0 += 4) {
      double u4[4];
      if (xp_vdot) {
        const double *pp[4], *zp[4];
        for (int k = 0; k < 4; k++) {
          const int s = s0 + k;
          pp[k] = t->PD + (size_t)s * NS;
          zp[k] = t->QP + ((size_t)t->cat[s] * nu + t->qc[s]) * NS;
        }
        double l4[4]; xp_dot4v(pp, zp, u4); xp_log4v(u4, l4);
        ll_tot += l4[0]; ll_tot += l4[1]; ll_tot += l4[2]; ll_tot += l4[3];
        continue;
      }
      for (int k = 0; k < 4; k++) {
        const int s = s0 + k;
        const double *pd = t->PD + (size_t)s * NS;
        const double *zz = t->QP + ((size_t)t->cat[s] * nu + t->qc[s]) * NS;
        double util = 0;
        for (int st = 0; st < NS; st++) util += pd[st] * zz[st];
        u4[k] = util;
      }
      double l4[4]; xp_log4v(u4, l4);
      ll_tot += l4[0]; ll_tot += l4[1]; ll_tot += l4[2]; ll_tot += l4[3];
    }
  }
#endif
  for (int s = s0; s < t->m; s++) {
    const double *pd = t->PD + (size_t)s * NS;
    const double *zz = t->QP + ((size_t)t->cat[s] * nu + t->qc[s]) * NS;
    double util = 0;
    for (int st = 0; st < NS; st++) util += pd[st] * zz[st];
    ll_tot += log(util);
  }
  t->ll_core = ll_tot;
}

static void tt_group_distal(TT *t)
{
  const int m = t->m, hb = t->dgrp;              /* t->dgrp holds the table size (a power of 2) */
  for (int h = 0; h < hb; h++) t->dhash[h] = -1;
  t->nd = 0;
  for (int s = 0; s < m; s++) {
    const uint64_t *w = (const uint64_t *)(t->dist_o + (size_t)s * NS);
    uint64_t h = 1469598103934665603ULL ^ (uint64_t)t->cat[s];
    for (int k = 0; k < NS; k++) { h ^= w[k]; h *= 1099511628211ULL; h ^= h >> 29; }
    int b = (int)(h & (uint64_t)(hb - 1));
    for (;;) {
      const int c = t->dhash[b];
      if (c < 0) { t->dhash[b] = t->nd; t->drep[t->nd] = s; t->dcls[s] = t->nd++; break; }
      const int r = t->drep[c];
      if (t->cat[r] == t->cat[s] && memcmp(t->dist_o + (size_t)r * NS, w, sizeof(double) * NS) == 0) { t->dcls[s] = c; break; }
      b = (b + 1) & (hb - 1);
    }
  }
}

static void tt_dist_pass(TT *t)                 /* PD and ll_core from dist_o/prox_o, QP */
{
  if (t->gmix) { gtt_dist_pass(t); return; }
  double ll_tot = 0;
  const int nu = t->nu;
  /* PD = (statd * prox) * dist per state, on every engine */
#if XP_SIMD20 && XP_AVX2
  const __m256d s0 = _mm256_loadu_pd(t->statd + 0), s1 = _mm256_loadu_pd(t->statd + 4),
                s2 = _mm256_loadu_pd(t->statd + 8), s3 = _mm256_loadu_pd(t->statd + 12),
                s4 = _mm256_loadu_pd(t->statd + 16);
#define XP_PD_SITE(pd, pv, dv) do { \
    _mm256_storeu_pd(pd + 0, _mm256_mul_pd(_mm256_mul_pd(s0, _mm256_loadu_pd(pv + 0)), _mm256_loadu_pd(dv + 0))); \
    _mm256_storeu_pd(pd + 4, _mm256_mul_pd(_mm256_mul_pd(s1, _mm256_loadu_pd(pv + 4)), _mm256_loadu_pd(dv + 4))); \
    _mm256_storeu_pd(pd + 8, _mm256_mul_pd(_mm256_mul_pd(s2, _mm256_loadu_pd(pv + 8)), _mm256_loadu_pd(dv + 8))); \
    _mm256_storeu_pd(pd + 12, _mm256_mul_pd(_mm256_mul_pd(s3, _mm256_loadu_pd(pv + 12)), _mm256_loadu_pd(dv + 12))); \
    _mm256_storeu_pd(pd + 16, _mm256_mul_pd(_mm256_mul_pd(s4, _mm256_loadu_pd(pv + 16)), _mm256_loadu_pd(dv + 16))); \
  } while (0)
#elif XP_SIMD20 && XP_NEON                        /* 2 states per op */
  const double *sd = t->statd;
  const float64x2_t s0 = vld1q_f64(sd + 0), s1 = vld1q_f64(sd + 2), s2 = vld1q_f64(sd + 4), s3 = vld1q_f64(sd + 6),
                    s4 = vld1q_f64(sd + 8), s5 = vld1q_f64(sd + 10), s6 = vld1q_f64(sd + 12), s7 = vld1q_f64(sd + 14),
                    s8 = vld1q_f64(sd + 16), s9 = vld1q_f64(sd + 18);
#define XP_NPD(pd, pv, dv, r) vst1q_f64(pd + 2 * r, vmulq_f64(vmulq_f64(s##r, vld1q_f64(pv + 2 * r)), vld1q_f64(dv + 2 * r)))
#define XP_PD_SITE(pd, pv, dv) do { \
    XP_NPD(pd, pv, dv, 0); XP_NPD(pd, pv, dv, 1); XP_NPD(pd, pv, dv, 2); XP_NPD(pd, pv, dv, 3); XP_NPD(pd, pv, dv, 4); \
    XP_NPD(pd, pv, dv, 5); XP_NPD(pd, pv, dv, 6); XP_NPD(pd, pv, dv, 7); XP_NPD(pd, pv, dv, 8); XP_NPD(pd, pv, dv, 9); \
  } while (0)
#else
#define XP_PD_SITE(pd, pv, dv) do { for (int st = 0; st < NS; st++) pd[st] = (t->statd[st] * pv[st]) * dv[st]; } while (0)
#endif
  const int vl = xp_vlog, vd = xp_vlog && xp_vdot;  /* both off on the portable engine */
  double u4[4]; int nb = 0;
  const double *pp4[4], *zp4[4];
  if (t->dgrp)                                     /* each distinct distal row evolved once */
    for (int c = 0; c < t->nd; c++) {
      const int r = t->drep[c];
      evolve_site(t->DV + (size_t)c * NS, t->dist_o + (size_t)r * NS, t->tensorD + (size_t)t->cat[r] * NS * NS);
    }
  for (int s = 0; s < t->m; s++) {
    double dvb[NS], pv[NS];
    const double *dv = dvb;
    const size_t off = (size_t)t->cat[s] * NS * NS;
    if (t->dgrp) dv = t->DV + (size_t)t->dcls[s] * NS;
    else evolve_site(dvb, t->dist_o + (size_t)s * NS, t->tensorD + off);
    evolve_site(pv, t->prox_o + (size_t)s * NS, t->tensorP + off);
    double *pd = t->PD + (size_t)s * NS;
    XP_PD_SITE(pd, pv, dv);
    const double *zz = t->QP + ((size_t)t->cat[s] * nu + t->qc[s]) * NS;
    if (vd) {                                      /* four sites' dots and logs at once */
      pp4[nb] = pd; zp4[nb] = zz; nb++;
      if (nb == 4) {
        double l4[4]; xp_dot4v(pp4, zp4, u4); xp_log4v(u4, l4);
        ll_tot += l4[0]; ll_tot += l4[1]; ll_tot += l4[2]; ll_tot += l4[3];
        nb = 0;
      }
      continue;
    }
    double util = 0;
    for (int st = 0; st < NS; st++) util += pd[st] * zz[st];
    if (vl) {                                      /* four sites' logs at once, added in order */
      u4[nb++] = util;
      if (nb == 4) {
        double l4[4]; xp_log4v(u4, l4);
        ll_tot += l4[0]; ll_tot += l4[1]; ll_tot += l4[2]; ll_tot += l4[3];
        nb = 0;
      }
      continue;
    }
    ll_tot += log(util);
  }
#undef XP_PD_SITE
#undef XP_NPD
  if (vd)                                          /* the last few sites, still in order */
    for (int k = 0; k < nb; k++) { double u = 0; for (int st = 0; st < NS; st++) u += pp4[k][st] * zp4[k][st]; ll_tot += log(u); }
  else
    for (int k = 0; k < nb; k++) ll_tot += log(u4[k]);
  t->ll_core = ll_tot;
}

static double tt_log_like(TT *t) { return t->ll_core + (M_LN2 * t->etot); }

static void tt_set_pend(TT *t, double v) { tt_build_QP(t, v); tt_ll_pass(t); }

/* Point t->PD at a slot for (dist, prox): returns 1 on a hit (PD already holds those values),
 * 0 after claiming a free or the least recently used slot for the caller to fill. */
static int pd_claim(TT *t, double d, double p)
{
  int victim = -1;
  for (int k = 0; k < t->npd; k++) {
    if (t->pd_ok[k] && t->pd_kd[k] == d && t->pd_kp[k] == p) {
      t->PD = t->pd_buf[k]; t->pd_age[k] = ++t->pd_clock; return 1;
    }
    if (!t->pd_ok[k]) { if (victim < 0 || t->pd_ok[victim]) victim = k; }
    else if (victim < 0 || (t->pd_ok[victim] && t->pd_age[k] < t->pd_age[victim])) victim = k;
  }
  t->pd_ok[victim] = 1; t->pd_kd[victim] = d; t->pd_kp[victim] = p; t->pd_age[victim] = ++t->pd_clock;
  t->PD = t->pd_buf[victim];
  return 0;
}

static void tt_set_dist(TT *t, double v)
{
  double cut = t->dist_bl + t->prox_bl;
  double prox = cut - v;
  if (!(prox >= 0.)) t->assert_fail = 1;           /* alpha19: assert(prox_bl >= 0.) */
  t->dist_bl = v;
  t->prox_bl = prox;
  if (t->pdc && !t->gmix && pd_claim(t, v, prox)) {
    /* A hit: PD is what tt_dist_pass would recompute from the same inputs; tt_ll_pass is the
     * same per-site dot and log, in the same order, that tt_dist_pass ends with. */
    tt_ll_pass(t);
    return;
  }
  tt_tensor(t, t->tensorD, v);
  tt_tensor(t, t->tensorP, prox);
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
  xp_read_switches();
  {
    const char *e = getenv("GPLACE_EXPERIMENTAL_PD_CACHE");
    t->pdc = (e && *e == '1' && m > 0);            /* m = 0 is the gamma path's placeholder */
  }
  if (t->pdc) {
    const char *e = getenv("GPLACE_PD_SLOTS");
    const int want = e ? atoi(e) : PD_SLOTS;
    t->npd = want < 1 ? 1 : want > PD_SLOTS ? PD_SLOTS : want;
    for (int k = 0; k < t->npd; k++) t->pd_buf[k] = malloc(sizeof(double) * n);
    t->PD = t->pd_buf[0];
  } else t->PD = malloc(sizeof(double) * n);
  t->dist_e = calloc(m, sizeof(int64_t)); t->prox_e = calloc(m, sizeof(int64_t));
  {
    const char *g = getenv("GPLACE_DGRP");
    if (g && *g == '1' && m > 0 && XP_AVX2 && XP_SIMD20) {
      int hb = 1; while (hb < 2 * m) hb <<= 1;
      t->dgrp = hb;
      t->dcls = malloc(sizeof(int32_t) * m); t->drep = malloc(sizeof(int32_t) * m);
      t->dhash = malloc(sizeof(int32_t) * hb); t->DV = malloc(sizeof(double) * n);
    }
  }
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
  if (t->pdc) for (int k = 0; k < t->npd; k++) free(t->pd_buf[k]);
  else free(t->PD);
  if (t->dgrp) { free(t->dcls); free(t->drep); free(t->dhash); free(t->DV); }
  free(t->dist_e); free(t->prox_e); free(t->U); free(t->qc); free(t->QP); free(t);
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
  if (t->dgrp) tt_group_distal(t);
  /* pendant first (Glv_edge.set_bl pend_edge), then dist edge, then prox edge */
  tt_build_QP(t, pendant);
  t->dist_bl = distal;
  tt_tensor(t, t->tensorD, distal);
  t->prox_bl = cut_bl - distal;
  tt_tensor(t, t->tensorP, cut_bl - distal);
  if (t->pdc && !t->gmix) {                        /* new edge vectors: nothing cached applies */
    for (int k = 0; k < t->npd; k++) t->pd_ok[k] = 0;
    pd_claim(t, t->dist_bl, t->prox_bl);
  }
  tt_dist_pass(t);
  t->assert_fail = 0;
}

void xp_tt_results(TT *t, double *out)   /* (log_like, pend_bl, dist_bl, n_calls) */
{
  out[0] = tt_log_like(t); out[1] = t->pend_bl; out[2] = t->dist_bl; out[3] = (double)t->n_calls;
}

/* Batch evaluation (GPLACE_CEVAL=1): the tiled placement's per-query work of
 * Placer._query_evaluator, for every query of one rebuilt tile or pack of pieces, on pthreads.
 * The same calls in the same order as the Python evaluator (get_d / get_p with the leaf proximal
 * rebuilt as Placer._prox_into does, xp_tt_set, xp_tt_optimize, the retry from the default
 * lengths after a GSL error), so the results are identical; it only replaces the forked worker
 * pool and the Python around those calls (GIL, copy-on-write, pool start-up). */
typedef struct {
  int S, nrates, max_iter;
  double start_pend, max_pend;
  const double *darr, *parr; const int64_t *dexp, *pexp;  /* tile slots: [slot][S][NS], [slot][S] */
  const int32_t *slot;                                     /* node -> slot, -1 if not rebuilt */
  const int32_t *parent, *ch_off, *ch, *leaf_row;          /* children of v: ch[ch_off[v]..ch_off[v+1]) */
  const double *bl, *lut;                                  /* lut: [256][NS] */
  const uint8_t *refA, *QM;                                /* [leaf rows][S], [queries][S] */
  const int32_t *cats; const double *statd, *x, *l, *xit, *rates; const uint8_t *occupied;
} EvCtx;

typedef struct {
  const EvCtx *c; int nq; const int32_t *qis; const int64_t *joff; const int32_t *loc;
  const double *pend, *dist, *tol; double *out; int32_t *status;
  int next, err;
} EvBatch;

typedef struct { double *D, *Pp, *b4, *b5, *b6, *tensor; int64_t *De, *Pe, *zero; } EvBuf;

static void ev_leaf_lv(const EvCtx *c, int v, double *dst)
{
  const uint8_t *r = c->refA + (size_t)c->leaf_row[v] * c->S;
  for (int s = 0; s < c->S; s++) memcpy(dst + (size_t)s * NS, c->lut + (size_t)r[s] * NS, sizeof(double) * NS);
}

static int ev_is_leaf(const EvCtx *c, int v) { return c->ch_off[v + 1] == c->ch_off[v]; }

/* Placer.d_pre: the distal vector before the final pull (a leaf's in buf, exponents zero) */
static int ev_dpre(const EvCtx *c, const EvBuf *b, int v, double *buf, const double **a, const int64_t **e)
{
  if (!ev_is_leaf(c, v)) {
    int k = c->slot[v]; if (k < 0) return -1;
    *a = c->darr + (size_t)k * c->S * NS; *e = c->dexp + (size_t)k * c->S;
  } else { ev_leaf_lv(c, v, buf); *a = buf; *e = b->zero; }
  return 0;
}

static void ev_evolve(const EvCtx *c, const EvBuf *b, double *dst, const double *src, double bl)
{
  xp_multi_exp(b->tensor, c->x, c->l, c->xit, c->rates, c->nrates, c->occupied, bl);
  xp_evolve(dst, src, c->S, c->cats, b->tensor);
}

static int ev_get_d(const EvCtx *c, const EvBuf *b, int v)       /* Placer.get_d */
{
  const int S = c->S;
  if (!ev_is_leaf(c, v)) {
    int k = c->slot[v]; if (k < 0) return -1;
    memcpy(b->D, c->darr + (size_t)k * S * NS, sizeof(double) * S * NS);
    memcpy(b->De, c->dexp + (size_t)k * S, sizeof(int64_t) * S);
  } else { ev_leaf_lv(c, v, b->D); memset(b->De, 0, sizeof(int64_t) * S); }
  xp_pull_exponent(b->D, b->De, S, -10);
  return 0;
}

static int ev_get_p(const EvCtx *c, const EvBuf *b, int v)       /* Placer.get_p, Placer._prox_into */
{
  const int S = c->S;
  if (!ev_is_leaf(c, v)) {
    int k = c->slot[v]; if (k < 0) return -1;
    memcpy(b->Pp, c->parr + (size_t)k * S * NS, sizeof(double) * S * NS);
    memcpy(b->Pe, c->pexp + (size_t)k * S, sizeof(int64_t) * S);
  } else {
    int pv = c->parent[v]; if (pv < 0) return -1;
    int kp = c->slot[pv]; if (kp < 0) return -1;
    const int64_t *evpe = c->pexp + (size_t)kp * S;
    ev_evolve(c, b, b->b4, c->parr + (size_t)kp * S * NS, c->bl[pv]);
    int first = 1;
    for (int j = c->ch_off[pv]; j < c->ch_off[pv + 1]; j++) {
      int r = c->ch[j]; if (r == v) continue;
      const double *a; const int64_t *e;
      if (ev_dpre(c, b, r, b->b6, &a, &e)) return -1;
      ev_evolve(c, b, b->b5, a, c->bl[r]);
      if (first) xp_pairwise(b->Pp, b->Pe, b->b4, evpe, b->b5, e, S);
      else xp_pairwise(b->Pp, b->Pe, b->Pp, b->Pe, b->b5, e, S);
      first = 0;
    }
    if (first) { memcpy(b->Pp, b->b4, sizeof(double) * S * NS); memcpy(b->Pe, evpe, sizeof(int64_t) * S); }
    xp_pull_exponent(b->Pp, b->Pe, S, -50);
  }
  xp_pull_exponent(b->Pp, b->Pe, S, -10);
  return 0;
}

static int ev_query(const EvBatch *B, const EvBuf *b, int i)
{
  const EvCtx *c = B->c;
  const int S = c->S, q = B->qis[i];
  const uint8_t *row = c->QM + (size_t)q * S;
  int32_t *reind = malloc(sizeof(int32_t) * S), *catm = malloc(sizeof(int32_t) * S);
  int m = 0;
  for (int s = 0; s < S; s++) if (row[s] != '-' && row[s] != '?') reind[m++] = s;  /* core.informative */
  double *qm = malloc(sizeof(double) * (size_t)(m ? m : 1) * NS);
  for (int k = 0; k < m; k++) {
    memcpy(qm + (size_t)k * NS, c->lut + (size_t)row[reind[k]] * NS, sizeof(double) * NS);
    catm[k] = c->cats[reind[k]];
  }
  int rc = 0;
  TT *t = xp_tt_new(m, catm, c->statd, c->x, c->l, c->xit, c->rates, c->nrates, qm);
  if (!t) rc = -2;
  int cur = -1; double res[4];
  for (int64_t j = B->joff[i]; !rc && j < B->joff[i + 1]; j++) {
    int loc = B->loc[j];
    double pend = B->pend[j], dist = B->dist[j];
    for (int attempt = 0; ; attempt++) {
      if (cur != loc) {
        if (ev_get_d(c, b, loc) || ev_get_p(c, b, loc)) { rc = -1; break; }
        cur = loc;
      }
      xp_tt_set(t, b->D, b->De, b->Pp, b->Pe, reind, c->bl[loc], pend, dist);
      int st = xp_tt_optimize(t, B->tol[j], c->max_pend, c->max_iter);
      if (st == XP_GSL && attempt == 0) {          /* safe_ml_optimize_location: retry once */
        pend = c->start_pend; dist = c->bl[loc] / 2.; continue;
      }
      B->status[j] = st;
      if (st == XP_OK || st == XP_EXCEEDED) {
        xp_tt_results(t, res);
        B->out[3 * j] = res[0]; B->out[3 * j + 1] = res[1]; B->out[3 * j + 2] = res[2];
      }
      break;
    }
  }
  if (t) xp_tt_free(t);
  free(qm); free(catm); free(reind);
  return rc;
}

static void *ev_worker(void *arg)
{
  EvBatch *B = (EvBatch *)arg;
  const EvCtx *c = B->c;
  const size_t V = (size_t)c->S * NS;
  EvBuf b;
  b.D = malloc(sizeof(double) * V); b.Pp = malloc(sizeof(double) * V);
  b.b4 = malloc(sizeof(double) * V); b.b5 = malloc(sizeof(double) * V); b.b6 = malloc(sizeof(double) * V);
  b.tensor = malloc(sizeof(double) * (size_t)c->nrates * NS * NS);
  b.De = malloc(sizeof(int64_t) * c->S); b.Pe = malloc(sizeof(int64_t) * c->S);
  b.zero = calloc(c->S, sizeof(int64_t));
  for (;;) {
    if (__atomic_load_n(&B->err, __ATOMIC_RELAXED)) break;
    int i = __atomic_fetch_add(&B->next, 1, __ATOMIC_RELAXED);
    if (i >= B->nq) break;
    int rc = ev_query(B, &b, i);
    if (rc) { __atomic_store_n(&B->err, rc, __ATOMIC_RELAXED); break; }
  }
  free(b.D); free(b.Pp); free(b.b4); free(b.b5); free(b.b6); free(b.tensor);
  free(b.De); free(b.Pe); free(b.zero);
  return NULL;
}

/* Queries qis[i] with jobs joff[i]..joff[i+1] (loc, pend, dist, tol); writes out[3j] = (like,
 * pend, dist) and status[j]. Returns 0, -1 (a vector the rebuild did not provide) or -2 (TT
 * allocation failed). */
int xp_eval_batch(const EvCtx *c, int nq, const int32_t *qis, const int64_t *joff, const int32_t *loc,
                  const double *pend, const double *dist, const double *tol, double *out,
                  int32_t *status, int nthreads)
{
  EvBatch B = { c, nq, qis, joff, loc, pend, dist, tol, out, status, 0, 0 };
  if (nthreads < 1) nthreads = 1;
  if (nthreads > nq) nthreads = nq > 0 ? nq : 1;
  pthread_t *th = malloc(sizeof(pthread_t) * nthreads);
  int started = 0;
  for (int k = 1; k < nthreads; k++)
    if (pthread_create(&th[k], NULL, ev_worker, &B) == 0) started = k; else break;
  ev_worker(&B);
  for (int k = 1; k <= started; k++) pthread_join(th[k], NULL);
  free(th);
  return B.err;
}

/* Candidate merge (GPLACE_CMERGE=1): placer.merge_candidates in C, rows on pthreads. Row r keeps
 * the first K of its K candidates (cs, ce) and the E new scores sc[r] of `edges`, in ranking
 * order: score descending, ties by ascending edge number, NaN after everything. The K kept are
 * written in no particular order (every caller orders the rows by (score, edge) at the end), so
 * the rows hold the same sets as the numpy version's. */
static inline int mg_better(double as, int64_t ae, double bs, int64_t be)
{
  int an = as != as, bn = bs != bs;
  if (an | bn) return an != bn ? bn : ae < be;
  return as != bs ? as > bs : ae < be;
}

/* Heap of K entries with the worst at the root: replace the root and sift it down. */
static void mg_sift(double *s, int64_t *e, int K, int i)
{
  for (;;) {
    int w = i, l = 2 * i + 1, r = l + 1;
    if (l < K && mg_better(s[w], e[w], s[l], e[l])) w = l;
    if (r < K && mg_better(s[w], e[w], s[r], e[r])) w = r;
    if (w == i) return;
    double ts = s[i]; s[i] = s[w]; s[w] = ts;
    int64_t te = e[i]; e[i] = e[w]; e[w] = te;
    i = w;
  }
}

typedef struct {
  double *cs; int64_t *ce; int K; int64_t nrows; const double *sc; int64_t ld;
  const int64_t *edges; int64_t E; int64_t next;
} MgJob;

static void *mg_worker(void *arg)
{
  MgJob *J = (MgJob *)arg;
  const int K = J->K;
  for (;;) {
    int64_t r0 = __atomic_fetch_add(&J->next, 64, __ATOMIC_RELAXED);
    if (r0 >= J->nrows) break;
    int64_t r1 = r0 + 64 < J->nrows ? r0 + 64 : J->nrows;
    for (int64_t r = r0; r < r1; r++) {
      double *s = J->cs + r * K; int64_t *e = J->ce + r * K;
      for (int i = K / 2 - 1; i >= 0; i--) mg_sift(s, e, K, i);
      const double *x = J->sc + r * J->ld;
      const int64_t *ed = J->edges, E = J->E;
      int64_t j = 0;
      /* Most scores are below the worst kept one (the root, which only rises): skip 8 at a time
       * when every one is strictly below it (ordered compare; NaN and ties take the exact test). */
#if XP_AVX2
      for (; j + 8 <= E; j += 8) {
        __m256d th = _mm256_set1_pd(s[0]);
        __m256d c = _mm256_or_pd(_mm256_cmp_pd(_mm256_loadu_pd(x + j), th, _CMP_NLT_UQ),
                                 _mm256_cmp_pd(_mm256_loadu_pd(x + j + 4), th, _CMP_NLT_UQ));
        if (!_mm256_movemask_pd(c)) continue;
        for (int k = 0; k < 8; k++)
          if (mg_better(x[j + k], ed[j + k], s[0], e[0])) { s[0] = x[j + k]; e[0] = ed[j + k]; mg_sift(s, e, K, 0); }
      }
#elif XP_NEON
      for (; j + 8 <= E; j += 8) {
        float64x2_t th = vdupq_n_f64(s[0]);
        uint64x2_t lt = vandq_u64(vandq_u64(vcltq_f64(vld1q_f64(x + j), th), vcltq_f64(vld1q_f64(x + j + 2), th)),
                                  vandq_u64(vcltq_f64(vld1q_f64(x + j + 4), th), vcltq_f64(vld1q_f64(x + j + 6), th)));
        if (vgetq_lane_u64(lt, 0) & vgetq_lane_u64(lt, 1)) continue;
        for (int k = 0; k < 8; k++)
          if (mg_better(x[j + k], ed[j + k], s[0], e[0])) { s[0] = x[j + k]; e[0] = ed[j + k]; mg_sift(s, e, K, 0); }
      }
#endif
      for (; j < E; j++)
        if (mg_better(x[j], ed[j], s[0], e[0])) { s[0] = x[j]; e[0] = ed[j]; mg_sift(s, e, K, 0); }
    }
  }
  return NULL;
}

void xp_merge_candidates(double *cs, int64_t *ce, int K, int64_t nrows, const double *sc, int64_t ld,
                         const int64_t *edges, int64_t E, int nthreads)
{
  MgJob J = { cs, ce, K, nrows, sc, ld, edges, E, 0 };
  if (K < 1 || nrows < 1) return;
  int64_t most = (nrows + 63) / 64;
  if (nthreads < 1) nthreads = 1;
  if (nthreads > most) nthreads = (int)most;
  pthread_t *th = malloc(sizeof(pthread_t) * nthreads);
  int started = 0;
  for (int k = 1; k < nthreads; k++)
    if (pthread_create(&th[k], NULL, mg_worker, &J) == 0) started = k; else break;
  mg_worker(&J);
  for (int k = 1; k <= started; k++) pthread_join(th[k], NULL);
  free(th);
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

/* Cache build, one node per call (GPLACE_CNODE=1). Placer._distal / _proximal make
 * several calls per node (xp_multi_exp, xp_evolve, xp_pairwise, xp_pull_exponent) with Python,
 * holding the GIL, in between; these do the same operations on the same operands, site by site,
 * in one call: per site, every child's evolve (evolve_site, the tensor from xp_multi_exp), the
 * products in the same order (x * y elementwise, exponents added), then the pull of
 * xp_pull_exponent. A leaf child's input is its tip vector lut[code], so its evolve depends only
 * on (rate category, code): it is computed once per pair and copied (identical inputs, identical
 * evolve_site). Children: src = distal rows (NULL for a leaf), srce = their exponents (NULL: 0),
 * codes = a leaf's tip codes (NULL otherwise). */
typedef struct { const double *src; const int64_t *srce; const uint8_t *codes; double *ten;
                 double *tab; uint8_t *have; } CnChild;

static int cn_setup(CnChild *c, int nch, const int64_t *src, const int64_t *srce, const int64_t *codes,
                    const double *bls, const double *x, const double *l, const double *xit,
                    const double *rates, int nrates, const uint8_t *occupied)
{
  for (int k = 0; k < nch; k++) {
    c[k].src = (const double *)(intptr_t)src[k]; c[k].srce = (const int64_t *)(intptr_t)srce[k];
    c[k].codes = (const uint8_t *)(intptr_t)codes[k];
    c[k].ten = (double *)malloc(sizeof(double) * (size_t)nrates * NS * NS);
    c[k].tab = NULL; c[k].have = NULL;
    if (!c[k].ten) return -1;
    xp_multi_exp(c[k].ten, x, l, xit, rates, nrates, occupied, bls[k]);
    if (c[k].codes) {
      c[k].tab = (double *)aligned_alloc(32, sizeof(double) * (size_t)nrates * 256 * NS);
      c[k].have = (uint8_t *)calloc((size_t)nrates * 256, 1);
      if (!c[k].tab || !c[k].have) return -1;
    }
  }
  return 0;
}
static void cn_free(CnChild *c, int nch)
{
  for (int k = 0; k < nch; k++) { free(c[k].ten); free(c[k].tab); free(c[k].have); }
}
/* child k's evolved distal at site s (into buf, or a pointer into its leaf table) */
static inline const double *cn_ev(CnChild *c, int s, int cs, const double *lut, double *buf)
{
  if (c->codes) {
    const size_t key = (size_t)cs * 256 + c->codes[s];
    double *t = c->tab + key * NS;
    if (!c->have[key]) { evolve_site(t, lut + (size_t)c->codes[s] * NS, c->ten + (size_t)cs * NS * NS); c->have[key] = 1; }
    return t;
  }
  evolve_site(buf, c->src + (size_t)s * NS, c->ten + (size_t)cs * NS * NS);
  return buf;
}
static inline int64_t cn_e(const CnChild *c, int s) { return c->srce ? c->srce[s] : 0; }
static inline void cn_pull(double *v, int64_t *e, int min_twoexp)   /* xp_pull_exponent, one site */
{
  double mx = v[0];
  for (int i = 1; i < NS; i++) if (v[i] > mx) mx = v[i];
  int twoexp;
  frexp(mx, &twoexp);
  if (twoexp < min_twoexp) {
    double f = ldexp(0.5, -twoexp + 1);
    for (int i = 0; i < NS; i++) v[i] *= f;
    *e += twoexp;
  }
}

/* A finished row written with non-temporal stores (AVX2; a plain copy on the other engines), so
 * the written lines are not first read into the cache: the tile walk's output rows. */
static inline void cn_store(double *dst, const double *v)
{
#if XP_AVX2
  if (!((uintptr_t)dst & 31)) {
    for (int i = 0; i < NS; i += 4) _mm256_stream_pd(dst + i, _mm256_load_pd(v + i));
    return;
  }
#endif
  memcpy(dst, v, sizeof(double) * NS);
}

/* Placer._distal: dst / de = the node's distal (pulled to -50) from its nch children */
int xp_cache_distal(double *dst, int64_t *de, int nch, const int64_t *src, const int64_t *srce,
                    const int64_t *codes, const double *bls, const double *lut, int nsites,
                    const int32_t *cat, const double *x, const double *l, const double *xit,
                    const double *rates, int nrates, const uint8_t *occupied)
{
  if (nch < 1 || nch > 64) return -1;
  CnChild c[64];
  memset(c, 0, sizeof c);
  if (cn_setup(c, nch, src, srce, codes, bls, x, l, xit, rates, nrates, occupied)) { cn_free(c, nch); return -1; }
  double bufa[NS] __attribute__((aligned(32))), bufb[NS] __attribute__((aligned(32)));
  for (int s = 0; s < nsites; s++) {
    const int cs = cat[s];
    double *d = dst + (size_t)s * NS;
    const double *a = cn_ev(&c[0], s, cs, lut, bufa);
    int64_t e = cn_e(&c[0], s);
    if (nch == 1) {
      memcpy(d, a, sizeof(double) * NS);
    } else {
      const double *b = cn_ev(&c[1], s, cs, lut, bufb);
      for (int i = 0; i < NS; i++) d[i] = a[i] * b[i];
      e = e + cn_e(&c[1], s);
      for (int k = 2; k < nch; k++) {
        const double *b2 = cn_ev(&c[k], s, cs, lut, bufb);
        for (int i = 0; i < NS; i++) d[i] = d[i] * b2[i];
        e = e + cn_e(&c[k], s);
      }
    }
    cn_pull(d, &e, -50);
    de[s] = e;
  }
  cn_free(c, nch);
  return 0;
}

/* Placer._proximal: from the node's proximal (par / pare, evolved by bl_v) and its children's
 * distals, the proximal of every child with out[k] != 0 (out = its row, oute = its exponents) */
int xp_cache_proximal(const double *par, const int64_t *pare, double bl_v, int nch, const int64_t *src,
                      const int64_t *srce, const int64_t *codes, const double *bls, const int64_t *out,
                      const int64_t *oute, const double *lut, int nsites, const int32_t *cat,
                      const double *x, const double *l, const double *xit, const double *rates,
                      int nrates, const uint8_t *occupied)
{
  if (nch < 1 || nch > 63) return -1;
  CnChild c[64];
  memset(c, 0, sizeof c);
  /* child slot nch: the node's own proximal, evolved by bl_v */
  int64_t src2[64], srce2[64], codes2[64]; double bls2[64];
  for (int k = 0; k < nch; k++) { src2[k] = src[k]; srce2[k] = srce[k]; codes2[k] = codes[k]; bls2[k] = bls[k]; }
  src2[nch] = (int64_t)(intptr_t)par; srce2[nch] = (int64_t)(intptr_t)pare; codes2[nch] = 0; bls2[nch] = bl_v;
  if (cn_setup(c, nch + 1, src2, srce2, codes2, bls2, x, l, xit, rates, nrates, occupied)) { cn_free(c, nch + 1); return -1; }
  double evp[NS] __attribute__((aligned(32))), bufs[64][NS] __attribute__((aligned(32)));
  const double *evk[64];
  for (int s = 0; s < nsites; s++) {
    const int cs = cat[s];
    const double *p = cn_ev(&c[nch], s, cs, lut, evp);
    const int64_t pe = cn_e(&c[nch], s);
    for (int k = 0; k < nch; k++) evk[k] = NULL;
    for (int t = 0; t < nch; t++) {
      if (!out[t]) continue;                        /* leaf child: rebuilt on demand by get_p */
      double *d = (double *)(intptr_t)out[t] + (size_t)s * NS;
      int64_t e;
      int first = 1;
      for (int r = 0; r < nch; r++) {
        if (r == t) continue;
        if (!evk[r]) evk[r] = cn_ev(&c[r], s, cs, lut, bufs[r]);
        if (first) { for (int i = 0; i < NS; i++) d[i] = p[i] * evk[r][i]; e = pe + cn_e(&c[r], s); first = 0; }
        else { for (int i = 0; i < NS; i++) d[i] = d[i] * evk[r][i]; e = e + cn_e(&c[r], s); }
      }
      if (first) { memcpy(d, p, sizeof(double) * NS); e = pe; }
      cn_pull(d, &e, -50);
      ((int64_t *)(intptr_t)oute[t])[s] = e;
    }
  }
  cn_free(c, nch + 1);
  return 0;
}

/* Tile build in units (GPLACE_TILEWALK=1). The cache build is bound by memory bandwidth: every
 * node's rows go to memory and its parent reads them back. Here one call builds a unit, a
 * connected group of a tile's internal nodes, in blocks of sites: within a block a node's result
 * stays in a cache-resident scratch slice that its parent (distal) or children (proximal) in the
 * unit read, and only rows entering the unit are read from memory. The transition tensors of the
 * unit's edges are made once per call and stay in cache across the blocks. Per site, the
 * operations, operands and order are those of xp_cache_distal / xp_cache_proximal, so every
 * value is identical; rows are written with non-temporal stores (cn_store).
 * Children of node j (in order) are at [off_j, off_j + nch[j]) of the child arrays:
 *   ckind 0 = a node of this unit (cref = its index, distal only), 1 = a row (cref / crefe =
 *   rows and exponents pointers), 2 = a leaf (cref = its tip codes). */
static __thread unsigned char *tw_buf = NULL;
static __thread size_t tw_cap = 0;
static void *tw_arena(size_t bytes)
{
  bytes = (bytes + 63) & ~(size_t)63;
  if (bytes > tw_cap) {
    free(tw_buf);
    tw_buf = (unsigned char *)aligned_alloc(64, bytes);
    tw_cap = tw_buf ? bytes : 0;
  }
  return tw_buf;
}

typedef struct { int kind; const double *src; const int64_t *srce; const uint8_t *codes;
                 double *ten; double *tab; uint8_t *have; } TwChild;

/* tensors (and leaf memo tables) for nc child edges, in the thread's arena after `extra` bytes */
static TwChild *tw_setup(int nc, const int32_t *ckind, const int64_t *cref, const int64_t *crefe,
                         const double *cbl, size_t extra, void **extra_out, const double *x,
                         const double *l, const double *xit, const double *rates, int nrates,
                         const uint8_t *occupied)
{
  size_t tsz = sizeof(double) * (size_t)nrates * NS * NS;
  size_t lsz = sizeof(double) * (size_t)nrates * 256 * NS, hsz = ((size_t)nrates * 256 + 63) & ~(size_t)63;
  int nleaf = 0;
  for (int k = 0; k < nc; k++) nleaf += ckind[k] == 2;
  size_t head = ((sizeof(TwChild) * (size_t)nc) + 63) & ~(size_t)63;
  extra = (extra + 63) & ~(size_t)63;
  unsigned char *b = (unsigned char *)tw_arena(head + extra + tsz * nc + (lsz + hsz) * nleaf);
  if (!b) return NULL;
  TwChild *c = (TwChild *)b;
  *extra_out = b + head;
  unsigned char *q = b + head + extra;
  for (int k = 0; k < nc; k++) {
    c[k].kind = ckind[k];
    c[k].src = ckind[k] == 1 ? (const double *)(intptr_t)cref[k] : NULL;
    c[k].srce = ckind[k] == 1 ? (const int64_t *)(intptr_t)crefe[k] : NULL;
    c[k].codes = ckind[k] == 2 ? (const uint8_t *)(intptr_t)cref[k] : NULL;
    c[k].ten = (double *)q; q += tsz;
    xp_multi_exp(c[k].ten, x, l, xit, rates, nrates, occupied, cbl[k]);
    c[k].tab = NULL; c[k].have = NULL;
  }
  for (int k = 0; k < nc; k++)
    if (ckind[k] == 2) {
      c[k].tab = (double *)q; q += lsz;
      c[k].have = (uint8_t *)q; q += hsz;
      memset(c[k].have, 0, (size_t)nrates * 256);
    }
  return c;
}
/* child k's evolved row at site s: a leaf's from its memo table (as cn_ev), a row's evolved
 * into buf from src (a local scratch slice, or the global row) */
static inline const double *tw_ev(TwChild *c, int s, int cs, const double *lut, const double *src,
                                  double *buf)
{
  if (c->codes) {
    const size_t key = (size_t)cs * 256 + c->codes[s];
    double *t = c->tab + key * NS;
    if (!c->have[key]) { evolve_site(t, lut + (size_t)c->codes[s] * NS, c->ten + (size_t)cs * NS * NS); c->have[key] = 1; }
    return t;
  }
  evolve_site(buf, src, c->ten + (size_t)cs * NS * NS);
  return buf;
}
static inline void tw_store(double *dst, const double *v) { cn_store(dst, v); }
/* cn_pull with frexp's exponent read from the bits for a positive normal maximum (the same
 * value: x = m 2^e with m in [0.5, 1)), frexp itself otherwise */
static inline void tw_pull(double *v, int64_t *e, int min_twoexp)
{
  double mx = v[0];
  for (int i = 1; i < NS; i++) if (v[i] > mx) mx = v[i];
  int twoexp;
  uint64_t bits; memcpy(&bits, &mx, 8);
  const int be = (int)((bits >> 52) & 0x7ff);
  if (be != 0 && be != 0x7ff) twoexp = be - 1022;
  else frexp(mx, &twoexp);
  if (twoexp < min_twoexp) {
    double f = ldexp(0.5, -twoexp + 1);
    for (int i = 0; i < NS; i++) v[i] *= f;
    *e += twoexp;
  }
}
static inline void tw_store_e(int64_t *dst, int64_t e)
{
#if XP_AVX2
  _mm_stream_si64((long long *)dst, (long long)e);
#else
  *dst = e;
#endif
}

/* distal rows of nu internal nodes in post-order (children before parents): dst[j] / dste[j]
 * (roots_only: only the unit's root, j = nu - 1, the one row read outside the unit) */
int xp_unit_distal(int nu, const int32_t *nch, const int32_t *ckind, const int64_t *cref,
                   const int64_t *crefe, const double *cbl, const int64_t *dst, const int64_t *dste,
                   const double *lut, int nsites, const int32_t *cat, const double *x, const double *l,
                   const double *xit, const double *rates, int nrates, const uint8_t *occupied, int block,
                   int roots_only)
{
  if (nu < 1 || block < 1) return -1;
  int nc = 0, off[nu + 1];
  for (int j = 0; j < nu; j++) {
    if (nch[j] < 1 || nch[j] > 64) return -1;
    off[j] = nc; nc += nch[j];
  }
  off[nu] = nc;
  const int B = block < nsites ? block : nsites;
  void *ex;
  TwChild *c = tw_setup(nc, ckind, cref, crefe, cbl, sizeof(double) * (size_t)nu * B * NS
                        + sizeof(int64_t) * (size_t)nu * B, &ex, x, l, xit, rates, nrates, occupied);
  if (!c) return -1;
  double *scr = (double *)ex;
  int64_t *scre = (int64_t *)(scr + (size_t)nu * B * NS);
  double bufa[NS] __attribute__((aligned(32))), bufb[NS] __attribute__((aligned(32)));
  for (int s0 = 0; s0 < nsites; s0 += B) {
    const int n = nsites - s0 < B ? nsites - s0 : B;
    for (int j = 0; j < nu; j++) {
      TwChild *cj = c + off[j];
      double *out = (double *)(intptr_t)dst[j];
      int64_t *oute = (int64_t *)(intptr_t)dste[j];
      for (int i = 0; i < n; i++) {
        const int s = s0 + i, cs = cat[s];
        double *d = scr + ((size_t)j * B + i) * NS;
#define TW_SRC(ch, k) ((ch)->kind == 0 ? scr + ((size_t)cref[off[j] + (k)] * B + i) * NS \
                       : (ch)->src ? (ch)->src + (size_t)s * NS : NULL)
#define TW_E(ch, k) ((ch)->kind == 0 ? scre[(size_t)cref[off[j] + (k)] * B + i] : (ch)->srce ? (ch)->srce[s] : 0)
        const double *a = tw_ev(&cj[0], s, cs, lut, TW_SRC(&cj[0], 0), bufa);
        int64_t e = TW_E(&cj[0], 0);
        if (nch[j] == 1) {
          memcpy(d, a, sizeof(double) * NS);
        } else {
          const double *b = tw_ev(&cj[1], s, cs, lut, TW_SRC(&cj[1], 1), bufb);
          for (int q = 0; q < NS; q++) d[q] = a[q] * b[q];
          e = e + TW_E(&cj[1], 1);
          for (int k = 2; k < nch[j]; k++) {
            const double *b2 = tw_ev(&cj[k], s, cs, lut, TW_SRC(&cj[k], k), bufb);
            for (int q = 0; q < NS; q++) d[q] = d[q] * b2[q];
            e = e + TW_E(&cj[k], k);
          }
        }
#undef TW_SRC
#undef TW_E
        tw_pull(d, &e, -50);
        scre[(size_t)j * B + i] = e;
        if (roots_only && j < nu - 1) continue;      /* read only within this unit */
        tw_store(out + (size_t)s * NS, d);
        tw_store_e(oute + s, e);
      }
    }
  }
#if XP_AVX2
  _mm_sfence();
#endif
  return 0;
}

/* proximal rows of nu internal nodes in pre-order (parents before children). Node j's own
 * proximal is the row pref[j] / prefe[j] when plocal[j] < 0 (the unit's root), otherwise the
 * scratch slice its parent wrote; it is evolved by pbl[j]. Children as above (ckind 1 or 2, the
 * child's distal); a child with cout != 0 gets its proximal there (and in its scratch slice
 * when clocal >= 0, its index in this unit). */
int xp_unit_proximal(int nu, const int32_t *plocal, const int64_t *pref, const int64_t *prefe,
                     const double *pbl, const int32_t *nch, const int32_t *ckind, const int64_t *cref,
                     const int64_t *crefe, const double *cbl, const int64_t *cout, const int64_t *coute,
                     const int32_t *clocal, const double *lut, int nsites, const int32_t *cat,
                     const double *x, const double *l, const double *xit, const double *rates,
                     int nrates, const uint8_t *occupied, int block)
{
  if (nu < 1 || block < 1) return -1;
  int nc = 0, off[nu + 1];
  for (int j = 0; j < nu; j++) {
    if (nch[j] < 1 || nch[j] > 63) return -1;
    off[j] = nc; nc += nch[j];
  }
  off[nu] = nc;
  /* the arena holds nc child tensors, then nu own-edge tensors (as rows, kind 1) */
  int32_t kk[nc + nu]; int64_t rr[nc + nu], re[nc + nu]; double bb[nc + nu];
  for (int k = 0; k < nc; k++) { kk[k] = ckind[k]; rr[k] = cref[k]; re[k] = crefe[k]; bb[k] = cbl[k]; }
  for (int j = 0; j < nu; j++) { kk[nc + j] = 1; rr[nc + j] = pref[j]; re[nc + j] = prefe[j]; bb[nc + j] = pbl[j]; }
  const int B = block < nsites ? block : nsites;
  void *ex;
  TwChild *c = tw_setup(nc + nu, kk, rr, re, bb, sizeof(double) * (size_t)nu * B * NS
                        + sizeof(int64_t) * (size_t)nu * B, &ex, x, l, xit, rates, nrates, occupied);
  if (!c) return -1;
  double *scr = (double *)ex;
  int64_t *scre = (int64_t *)(scr + (size_t)nu * B * NS);
  double evp[NS] __attribute__((aligned(32))), bufs[64][NS] __attribute__((aligned(32)));
  double tmp[NS] __attribute__((aligned(32)));
  const double *evk[64];
  for (int s0 = 0; s0 < nsites; s0 += B) {
    const int n = nsites - s0 < B ? nsites - s0 : B;
    for (int j = 0; j < nu; j++) {
      TwChild *cj = c + off[j], *own = c + nc + j;
      const int m = nch[j];
      for (int i = 0; i < n; i++) {
        const int s = s0 + i, cs = cat[s];
        const double *psrc = plocal[j] < 0 ? own->src + (size_t)s * NS : scr + ((size_t)j * B + i) * NS;
        const int64_t pe = plocal[j] < 0 ? own->srce[s] : scre[(size_t)j * B + i];
        const double *p = tw_ev(own, s, cs, lut, psrc, evp);
        for (int k = 0; k < m; k++) evk[k] = NULL;
        for (int t = 0; t < m; t++) {
          const int64_t ot = cout[off[j] + t];
          if (!ot) continue;                          /* leaf child: rebuilt on demand by get_p */
          const int lc = clocal[off[j] + t];
          double *d = lc >= 0 ? scr + ((size_t)lc * B + i) * NS : tmp;
          int64_t e;
          int first = 1;
          for (int r = 0; r < m; r++) {
            if (r == t) continue;
            if (!evk[r]) evk[r] = tw_ev(&cj[r], s, cs, lut, cj[r].src ? cj[r].src + (size_t)s * NS : NULL, bufs[r]);
            const int64_t er = cj[r].srce ? cj[r].srce[s] : 0;
            if (first) { for (int q = 0; q < NS; q++) d[q] = p[q] * evk[r][q]; e = pe + er; first = 0; }
            else { for (int q = 0; q < NS; q++) d[q] = d[q] * evk[r][q]; e = e + er; }
          }
          if (first) { memcpy(d, p, sizeof(double) * NS); e = pe; }
          tw_pull(d, &e, -50);
          if (lc >= 0) scre[(size_t)lc * B + i] = e;
          tw_store((double *)(intptr_t)ot + (size_t)s * NS, d);
          tw_store_e((int64_t *)(intptr_t)coute[off[j] + t] + s, e);
        }
      }
    }
  }
#if XP_AVX2
  _mm_sfence();
#endif
  return 0;
}

/* A tile's units in one call (GPLACE_TILEWALK=1; distal 2: only the units' root rows are
 * written, for a distal pass that keeps only boundary vectors): nthreads threads (the caller's included) run
 * the units of one pass as they become ready, distal bottom-up (a unit once every unit below it
 * is done), proximal top-down (a unit's children once it is done). Unit u's nodes are
 * [uno[u], uno[u+1]) of the node arrays and its children [uco[u], uco[u+1]) of the child arrays,
 * laid out as xp_unit_distal / xp_unit_proximal take them; uparent[u] is the unit above it
 * (-1: none), its units below are ukids[ukid_off[u] .. ukid_off[u+1]). */
typedef struct {
  int distal, nunits, block, nsites, nrates;
  const int32_t *uparent, *ukid_off, *ukids; const int64_t *uno, *uco;
  const int32_t *nch, *plocal; const int64_t *dst, *dste, *pref, *prefe; const double *pbl;
  const int32_t *ckind, *clocal; const int64_t *cref, *crefe, *cout, *coute; const double *cbl;
  const double *lut; const int32_t *cat; const double *x, *l, *xit, *rates; const uint8_t *occupied;
  pthread_mutex_t mu; pthread_cond_t cv;
  int32_t *stack; int top, remaining, err; int32_t *pending;
} TwWalk;

static int tw_unit(TwWalk *w, int u)
{
  const int64_t no = w->uno[u], co = w->uco[u];
  const int nu = (int)(w->uno[u + 1] - no);
  if (w->distal)
    return xp_unit_distal(nu, w->nch + no, w->ckind + co, w->cref + co, w->crefe + co, w->cbl + co,
                          w->dst + no, w->dste + no, w->lut, w->nsites, w->cat, w->x, w->l, w->xit,
                          w->rates, w->nrates, w->occupied, w->block, w->distal == 2);
  return xp_unit_proximal(nu, w->plocal + no, w->pref + no, w->prefe + no, w->pbl + no, w->nch + no,
                          w->ckind + co, w->cref + co, w->crefe + co, w->cbl + co, w->cout + co,
                          w->coute + co, w->clocal + co, w->lut, w->nsites, w->cat, w->x, w->l, w->xit,
                          w->rates, w->nrates, w->occupied, w->block);
}

static void *tw_worker(void *arg)
{
  TwWalk *w = (TwWalk *)arg;
  pthread_mutex_lock(&w->mu);
  for (;;) {
    while (w->top == 0 && w->remaining > 0 && !w->err) pthread_cond_wait(&w->cv, &w->mu);
    if (w->remaining == 0 || w->err) break;
    const int u = w->stack[--w->top];
    pthread_mutex_unlock(&w->mu);
    const int st = tw_unit(w, u);
    pthread_mutex_lock(&w->mu);
    if (st) { w->err = st; pthread_cond_broadcast(&w->cv); break; }
    w->remaining--;
    int pushed = 0;
    if (w->distal) {
      const int p = w->uparent[u];
      if (p >= 0 && --w->pending[p] == 0) { w->stack[w->top++] = p; pushed = 1; }
    } else {
      for (int k = w->ukid_off[u]; k < w->ukid_off[u + 1]; k++) { w->stack[w->top++] = w->ukids[k]; pushed++; }
    }
    if (w->remaining == 0 || pushed > 1) pthread_cond_broadcast(&w->cv);
    else if (pushed) pthread_cond_signal(&w->cv);
  }
  pthread_mutex_unlock(&w->mu);
  return NULL;
}

static void *tw_thread(void *arg)
{
  tw_worker(arg);
  free(tw_buf); tw_buf = NULL; tw_cap = 0;               /* this thread's arena */
  return NULL;
}

int xp_tile_walk(int distal, int nunits, const int32_t *uparent, const int32_t *ukid_off,
                 const int32_t *ukids, const int64_t *uno, const int64_t *uco, const int32_t *nch,
                 const int32_t *plocal, const int64_t *dst, const int64_t *dste, const int64_t *pref,
                 const int64_t *prefe, const double *pbl, const int32_t *ckind, const int32_t *clocal,
                 const int64_t *cref, const int64_t *crefe, const int64_t *cout, const int64_t *coute,
                 const double *cbl, const double *lut, int nsites, const int32_t *cat, const double *x,
                 const double *l, const double *xit, const double *rates, int nrates,
                 const uint8_t *occupied, int block, int nthreads)
{
  if (nunits < 1) return 0;
  TwWalk w = {distal, nunits, block, nsites, nrates, uparent, ukid_off, ukids, uno, uco, nch, plocal,
              dst, dste, pref, prefe, pbl, ckind, clocal, cref, crefe, cout, coute, cbl, lut, cat, x, l,
              xit, rates, occupied};
  w.stack = (int32_t *)malloc(sizeof(int32_t) * (size_t)nunits);
  w.pending = (int32_t *)malloc(sizeof(int32_t) * (size_t)nunits);
  if (!w.stack || !w.pending) { free(w.stack); free(w.pending); return -1; }
  w.top = 0; w.remaining = nunits; w.err = 0;
  for (int u = 0; u < nunits; u++) {
    w.pending[u] = ukid_off[u + 1] - ukid_off[u];
    if (distal ? w.pending[u] == 0 : uparent[u] < 0) w.stack[w.top++] = u;
  }
  pthread_mutex_init(&w.mu, NULL); pthread_cond_init(&w.cv, NULL);
  if (nthreads < 1) nthreads = 1;
  if (nthreads > 256) nthreads = 256;
  pthread_t th[256];
  int started = 0;
  for (int t = 1; t < nthreads; t++)
    if (pthread_create(&th[started], NULL, tw_thread, &w) == 0) started++;
  tw_worker(&w);                                         /* the caller's thread works too */
  for (int t = 0; t < started; t++) pthread_join(th[t], NULL);
  pthread_mutex_destroy(&w.mu); pthread_cond_destroy(&w.cv);
  const int err = w.err ? w.err : (w.remaining ? -2 : 0);
  free(w.stack); free(w.pending);
  return err;
}

/* The f32 prefilter kernels below are written once, on 8-lane f32 / int32 vectors (pf8f, pf8i) that
 * each SIMD engine implements: one AVX2 register, or two NEON registers (lanes 0-3, 4-7). Every
 * operation is lane-wise and IEEE-exact the same on both (the fused multiply-add included), so the
 * two engines compute the same bits. Only the 20-state engine has them: its rows are 20 floats. */
#if (XP_AVX2 || XP_NEON) && NS == 20               /* 20 states (protein); other engines use f64 */
#define XP_PF8 1
#else
#define XP_PF8 0
#endif
#if XP_AVX2
typedef __m256 pf8f;
typedef __m256i pf8i;
static inline pf8f pf8f_set1(float v) { return _mm256_set1_ps(v); }
static inline pf8f pf8f_zero(void) { return _mm256_setzero_ps(); }
static inline pf8f pf8f_load(const float *p) { return _mm256_load_ps(p); }      /* 32-byte aligned */
static inline void pf8f_store(float *p, pf8f a) { _mm256_store_ps(p, a); }
static inline void pf8f_storeu(float *p, pf8f a) { _mm256_storeu_ps(p, a); }
static inline pf8f pf8f_mul(pf8f a, pf8f b) { return _mm256_mul_ps(a, b); }
static inline pf8f pf8f_fmadd(pf8f a, pf8f b, pf8f c) { return _mm256_fmadd_ps(a, b, c); }   /* a b + c */
static inline pf8f pf8f_max(pf8f a, pf8f b) { return _mm256_max_ps(a, b); }
static inline pf8f pf8f_cmpge(pf8f a, pf8f b) { return _mm256_cmp_ps(a, b, _CMP_GE_OQ); } /* false if NaN */
static inline pf8f pf8f_cmpeq(pf8f a, pf8f b) { return _mm256_cmp_ps(a, b, _CMP_EQ_OQ); }
static inline pf8f pf8f_blend(pf8f a, pf8f b, pf8f m) { return _mm256_blendv_ps(a, b, m); }  /* m ? b : a */
static inline float pf8f_hmax(pf8f m)
{
  __m128 h = _mm_max_ps(_mm256_castps256_ps128(m), _mm256_extractf128_ps(m, 1));
  h = _mm_max_ps(h, _mm_movehl_ps(h, h)); h = _mm_max_ss(h, _mm_shuffle_ps(h, h, 1));
  return _mm_cvtss_f32(h);
}
static inline pf8i pf8i_set1(int32_t v) { return _mm256_set1_epi32(v); }
static inline pf8i pf8i_zero(void) { return _mm256_setzero_si256(); }
static inline pf8i pf8i_loadu(const int32_t *p) { return _mm256_loadu_si256((const __m256i *)p); }
static inline pf8i pf8i_add(pf8i a, pf8i b) { return _mm256_add_epi32(a, b); }
static inline pf8i pf8i_sub(pf8i a, pf8i b) { return _mm256_sub_epi32(a, b); }
static inline pf8i pf8i_and(pf8i a, pf8i b) { return _mm256_and_si256(a, b); }
static inline pf8i pf8i_or(pf8i a, pf8i b) { return _mm256_or_si256(a, b); }
static inline pf8i pf8i_srl23(pf8i a) { return _mm256_srli_epi32(a, 23); }      /* logical */
static inline pf8i pf8f_as_i(pf8f a) { return _mm256_castps_si256(a); }
static inline pf8f pf8i_as_f(pf8i a) { return _mm256_castsi256_ps(a); }
static inline void pf8_add64(int64_t *E, pf8i e) { ADD64(E, e); }   /* 8 int32 into 8 int64 (exact) */
#elif XP_NEON
typedef struct { float32x4_t lo, hi; } pf8f;
typedef struct { int32x4_t lo, hi; } pf8i;
#define PF8F2(x, y) ((pf8f){ (x), (y) })
#define PF8I2(x, y) ((pf8i){ (x), (y) })
static inline pf8f pf8f_set1(float v) { return PF8F2(vdupq_n_f32(v), vdupq_n_f32(v)); }
static inline pf8f pf8f_zero(void) { return pf8f_set1(0.0f); }
static inline pf8f pf8f_load(const float *p) { return PF8F2(vld1q_f32(p), vld1q_f32(p + 4)); }
static inline void pf8f_store(float *p, pf8f a) { vst1q_f32(p, a.lo); vst1q_f32(p + 4, a.hi); }
static inline void pf8f_storeu(float *p, pf8f a) { pf8f_store(p, a); }
static inline pf8f pf8f_mul(pf8f a, pf8f b) { return PF8F2(vmulq_f32(a.lo, b.lo), vmulq_f32(a.hi, b.hi)); }
static inline pf8f pf8f_fmadd(pf8f a, pf8f b, pf8f c)
{ return PF8F2(vfmaq_f32(c.lo, a.lo, b.lo), vfmaq_f32(c.hi, a.hi, b.hi)); }
static inline pf8f pf8f_max(pf8f a, pf8f b) { return PF8F2(vmaxq_f32(a.lo, b.lo), vmaxq_f32(a.hi, b.hi)); }
static inline pf8f pf8f_cmpge(pf8f a, pf8f b)
{ return PF8F2(vreinterpretq_f32_u32(vcgeq_f32(a.lo, b.lo)), vreinterpretq_f32_u32(vcgeq_f32(a.hi, b.hi))); }
static inline pf8f pf8f_cmpeq(pf8f a, pf8f b)
{ return PF8F2(vreinterpretq_f32_u32(vceqq_f32(a.lo, b.lo)), vreinterpretq_f32_u32(vceqq_f32(a.hi, b.hi))); }
static inline pf8f pf8f_blend(pf8f a, pf8f b, pf8f m)
{ return PF8F2(vbslq_f32(vreinterpretq_u32_f32(m.lo), b.lo, a.lo), vbslq_f32(vreinterpretq_u32_f32(m.hi), b.hi, a.hi)); }
static inline float pf8f_hmax(pf8f m) { return vmaxvq_f32(vmaxq_f32(m.lo, m.hi)); }
static inline pf8i pf8i_set1(int32_t v) { return PF8I2(vdupq_n_s32(v), vdupq_n_s32(v)); }
static inline pf8i pf8i_zero(void) { return pf8i_set1(0); }
static inline pf8i pf8i_loadu(const int32_t *p) { return PF8I2(vld1q_s32(p), vld1q_s32(p + 4)); }
static inline pf8i pf8i_add(pf8i a, pf8i b) { return PF8I2(vaddq_s32(a.lo, b.lo), vaddq_s32(a.hi, b.hi)); }
static inline pf8i pf8i_sub(pf8i a, pf8i b) { return PF8I2(vsubq_s32(a.lo, b.lo), vsubq_s32(a.hi, b.hi)); }
static inline pf8i pf8i_and(pf8i a, pf8i b) { return PF8I2(vandq_s32(a.lo, b.lo), vandq_s32(a.hi, b.hi)); }
static inline pf8i pf8i_or(pf8i a, pf8i b) { return PF8I2(vorrq_s32(a.lo, b.lo), vorrq_s32(a.hi, b.hi)); }
static inline pf8i pf8i_srl23(pf8i a)
{
  return PF8I2(vreinterpretq_s32_u32(vshrq_n_u32(vreinterpretq_u32_s32(a.lo), 23)),
               vreinterpretq_s32_u32(vshrq_n_u32(vreinterpretq_u32_s32(a.hi), 23)));
}
static inline pf8i pf8f_as_i(pf8f a) { return PF8I2(vreinterpretq_s32_f32(a.lo), vreinterpretq_s32_f32(a.hi)); }
static inline pf8f pf8i_as_f(pf8i a) { return PF8F2(vreinterpretq_f32_s32(a.lo), vreinterpretq_f32_s32(a.hi)); }
static inline void pf8_add64(int64_t *E, pf8i e) { NADD64(E, e.lo, e.hi); }
#endif
#if XP_PF8
/* The f32 kernels' query pass over one block of sites s0..s1 (table rows Tm/Te, row j = s - s0):
 * every query's 8-lane product times the table mantissas, in site order, renormalised every
 * XP_RBLK sites (exact: a power-of-two scaling; a zero product stays zero); the exponents are
 * summed exactly into Eq. Four queries at a time, their products interleaved: per query the same
 * multiplications in the same order as one at a time, so the same values. Row k of the table is
 * Tm[2k] (mantissas) next to Te[2k] (exponents): one 64-byte line per gathered row. */
#define PF_RENORM(P, E) do { \
    const pf8f z_ = pf8f_cmpeq(P, pf8f_zero()); \
    const pf8i b_ = pf8f_as_i(P); \
    E = pf8i_add(E, pf8i_sub(pf8i_srl23(b_), pf8i_set1(127))); \
    P = pf8i_as_f(pf8i_or(pf8i_and(b_, pf8i_set1(0x007FFFFF)), \
                                            pf8i_set1(0x3F800000))); \
    P = pf8f_blend(P, pf8f_zero(), z_); \
  } while (0)
static void pf32_queries(const pf8f *Tm, const pf8i *Te, int NC, int s0, int s1,
                         const uint8_t *codes, int nsites, int nq, pf8f *Pq, int64_t *Eq)
{
  int q = 0;
  for (; q + 4 <= nq; q += 4) {
    const uint8_t *c0 = codes + (size_t)q * nsites, *c1 = c0 + nsites, *c2 = c1 + nsites, *c3 = c2 + nsites;
    pf8f P0 = Pq[q], P1 = Pq[q + 1], P2 = Pq[q + 2], P3 = Pq[q + 3];
    for (int sb = s0; sb < s1; sb += XP_RBLK) {
      const int sb1 = sb + XP_RBLK < s1 ? sb + XP_RBLK : s1;
      pf8i E0 = pf8i_zero(), E1 = E0, E2 = E0, E3 = E0;
      for (int s = sb; s < sb1; s++) {
        const size_t r = (size_t)(s - s0) * NC;
        const size_t k0 = r + c0[s], k1 = r + c1[s], k2 = r + c2[s], k3 = r + c3[s];
        P0 = pf8f_mul(P0, Tm[2 * k0]); E0 = pf8i_add(E0, Te[2 * k0]);
        P1 = pf8f_mul(P1, Tm[2 * k1]); E1 = pf8i_add(E1, Te[2 * k1]);
        P2 = pf8f_mul(P2, Tm[2 * k2]); E2 = pf8i_add(E2, Te[2 * k2]);
        P3 = pf8f_mul(P3, Tm[2 * k3]); E3 = pf8i_add(E3, Te[2 * k3]);
      }
      PF_RENORM(P0, E0); PF_RENORM(P1, E1); PF_RENORM(P2, E2); PF_RENORM(P3, E3);
      pf8_add64(Eq + (size_t)q * XP_EB, E0); pf8_add64(Eq + (size_t)(q + 1) * XP_EB, E1);
      pf8_add64(Eq + (size_t)(q + 2) * XP_EB, E2); pf8_add64(Eq + (size_t)(q + 3) * XP_EB, E3);
    }
    Pq[q] = P0; Pq[q + 1] = P1; Pq[q + 2] = P2; Pq[q + 3] = P3;
  }
  for (; q < nq; q++) {
    const uint8_t *cq = codes + (size_t)q * nsites;
    pf8f P = Pq[q];
    for (int sb = s0; sb < s1; sb += XP_RBLK) {
      const int sb1 = sb + XP_RBLK < s1 ? sb + XP_RBLK : s1;
      pf8i E = pf8i_zero();
      for (int s = sb; s < sb1; s++) {
        const size_t k = (size_t)(s - s0) * NC + cq[s];
        P = pf8f_mul(P, Tm[2 * k]); E = pf8i_add(E, Te[2 * k]);
      }
      PF_RENORM(P, E);
      pf8_add64(Eq + (size_t)q * XP_EB, E);
    }
    Pq[q] = P;
  }
}
#endif

/* Prefilter in f32 product form (GPLACE_PF32=1). The prefilter only orders edges, and
 * log is monotone, so a query's score against an edge is ranked by the product of its site utils
 * as well as by their log-sum. Here every util is an f32 dot of f32-rounded inputs (all terms are
 * non-negative, so its relative error is at most ~22 * 2^-24), the product is kept as an f32
 * mantissa in [1, 2) times an integer power of two (no logs, no underflow), and the result is
 * written as ln(mantissa) + ln2 * exponent. The caller adds a per-query radius
 * n_used * 26 * 2^-24 + 1e-6 that bounds |this - the f64 score of xp_edge_scores8| and resolves
 * every edge whose interval could change the ranking with the f64 kernel. Only the codes present
 * in the batch at a site (bitmask `present`) are evaluated; sites are processed in blocks (xp_pf_bs32)
 * so the table stays in cache. A util below 2^-90 (or zero) is not bounded: the score is NaN and
 * the caller computes it in f64. Returns 1 (nothing written) if an exponent is out of range. */
#if XP_PF8
int xp_edge_scores8_f32(double *out, long ldo, long col0, int ne,
                        const double *sn8, const int64_t *se8,
                        const double *qevtab, const int32_t *cat, int nsites, int ncode, int nrates,
                        const uint8_t *codes, int nq, const uint32_t *present)
{
  if (ne < 1 || ne > XP_EB || ncode < 1 || ncode > 31) return -1;
  for (int lane = 0; lane < ne; lane++)
    for (int s = 0; s < nsites; s++) {
      const int64_t v = se8[(size_t)lane * nsites + s];
      if (v > (1 << 23) || v < -(1 << 23)) return 1;
    }
  const int NC = ncode + 1, BS = xp_pf_bs32;
  float *q32 = (float *)malloc(sizeof(float) * (size_t)nrates * ncode * NS);
  /* table row k: its 8 mantissas Tm[2k] and exponents Te[2k] side by side, one 64-byte line */
  pf8f *Tm = (pf8f *)aligned_alloc(64, sizeof(pf8f) * 2 * (size_t)BS * NC);
  pf8i *Te = (pf8i *)(Tm + 1);
  pf8f *Pq = (pf8f *)aligned_alloc(32, sizeof(pf8f) * ((size_t)nq + 1));
  int64_t *Eq = (int64_t *)calloc((size_t)nq * XP_EB + 8, sizeof(int64_t));
  if (!q32 || !Tm || !Te || !Pq || !Eq) { free(q32); free(Tm); free(Pq); free(Eq); return -1; }
  for (size_t i = 0; i < (size_t)nrates * ncode * NS; i++) q32[i] = (float)qevtab[i];
  const pf8f one = pf8f_set1(1.0f), zero = pf8f_zero(), tiny = pf8f_set1(0x1p-90f);
  const pf8i mant = pf8i_set1(0x007FFFFF), onebits = pf8i_set1(0x3F800000),
                bias = pf8i_set1(127);
  for (int q = 0; q < nq; q++) Pq[q] = one;
  float snT[NS][XP_EB] __attribute__((aligned(32)));
  for (int s0 = 0; s0 < nsites; s0 += BS) {
    const int s1 = s0 + BS < nsites ? s0 + BS : nsites;
    for (int s = s0; s < s1; s++) {
      const int j = s - s0;
      for (int lane = 0; lane < XP_EB; lane++) {
        const double *sn = sn8 + ((size_t)(lane < ne ? lane : 0) * nsites + s) * NS;
        for (int st = 0; st < NS; st++) snT[st][lane] = (float)sn[st];
      }
      int32_t sev[XP_EB];
      for (int lane = 0; lane < XP_EB; lane++) sev[lane] = lane < ne ? (int32_t)se8[(size_t)lane * nsites + s] : 0;
      const pf8i sv = pf8i_loadu(sev);
      const float *qt = q32 + (size_t)cat[s] * ncode * NS;
      uint32_t pm = present[s];
      while (pm) {
        const int c = __builtin_ctz(pm); pm &= pm - 1;
        const float *qv = qt + (size_t)c * NS;
        pf8f acc = zero;
        for (int st = 0; st < NS; st++)
          acc = pf8f_fmadd(pf8f_set1(qv[st]), pf8f_load(snT[st]), acc);
        const pf8f ok = pf8f_cmpge(acc, tiny);       /* false for NaN too */
        const pf8i b = pf8f_as_i(acc);
        const pf8i e = pf8i_sub(pf8i_srl23(b), bias);
        const pf8f m = pf8i_as_f(pf8i_or(pf8i_and(b, mant), onebits));
        Tm[2 * ((size_t)j * NC + c)] = pf8f_blend(zero, m, ok);          /* unbounded -> 0 */
        Te[2 * ((size_t)j * NC + c)] = pf8i_add(e, sv);
      }
      Tm[2 * ((size_t)j * NC + ncode)] = one; Te[2 * ((size_t)j * NC + ncode)] = pf8i_zero();
    }
    pf32_queries(Tm, Te, NC, s0, s1, codes, nsites, nq, Pq, Eq);
  }
  for (int q = 0; q < nq; q++) {
    float p8[XP_EB]; pf8f_storeu(p8, Pq[q]);
    for (int lane = 0; lane < ne; lane++)
      out[(size_t)q * ldo + col0 + lane] = p8[lane] == 0.0f ? NAN
        : log((double)p8[lane]) + M_LN2 * (double)Eq[(size_t)q * XP_EB + lane];
  }
  free(q32); free(Tm); free(Pq); free(Eq);
  return 0;
}
#else
int xp_edge_scores8_f32(double *out, long ldo, long col0, int ne, const double *sn8, const int64_t *se8,
                        const double *qevtab, const int32_t *cat, int nsites, int ncode, int nrates,
                        const uint8_t *codes, int nq, const uint32_t *present)
{ (void)out; (void)ldo; (void)col0; (void)ne; (void)sn8; (void)se8; (void)qevtab; (void)cat; (void)nsites;
  (void)ncode; (void)nrates; (void)codes; (void)nq; (void)present; return 1; }
#endif

/* f32 supernodes for the product-form prefilter (GPLACE_PF32SRC=1, with GPLACE_PF32=1). As
 * xp_edge_scores8_f32, but each lane's supernode is built here, in f32, straight from the cache:
 *   dv = evolve(d, half), pv = evolve(p, half), sn = statd * dv * pv           (as _supernodes8)
 * where d is the edge's distal vector (a leaf's tip vector) and p its proximal vector or, for a
 * leaf edge, the rebuild of get_p: p = evolve(parent proximal, bl[parent]) * prod over siblings
 * of evolve(sibling distal, bl[sibling]). Every f64 input row is scaled by a power of two to a
 * maximum in [1, 2) before rounding to f32, and every product vector is renormalised the same way,
 * so the exponents are tracked here (a power-of-two scaling is exact) instead of gPlace's pulls,
 * which only change the representation. Every operation adds or multiplies non-negative numbers,
 * so each value carries at most one relative rounding of 2^-24 per operation on its path; the
 * caller's radius counts them (pf32_rounds). Lane descriptor, 32 int64 each:
 *   [0] d  [1] d exps (NULL: 0)  [2] p mode (0 stored, 1 leaf rebuild)  [3] p  [4] p exps
 *   [5] half tensor  [6] parent tensor  [7] n siblings (<= 8)  [8..15] sibling distals
 *   [16..23] their exps (NULL: 0)  [24..31] their tensors
 * Tensors are R x NS x NS f64, transposed as from xp_multi_exp. Returns 1 if an exponent is out
 * of range (the caller then uses the f64 supernodes). */
#if XP_PF8
#define PF_W 24                                     /* a 20-float row padded to three vectors */
static inline void pf_ev(float *y, const float *x, const float *T)   /* y = x . T (T: 20 x PF_W) */
{
  pf8f a0 = pf8f_zero(), a1 = a0, a2 = a0;
  for (int j = 0; j < NS; j++) {
    const pf8f xj = pf8f_set1(x[j]);
    const float *t = T + j * PF_W;
    a0 = pf8f_fmadd(xj, pf8f_load(t), a0);
    a1 = pf8f_fmadd(xj, pf8f_load(t + 8), a1);
    a2 = pf8f_fmadd(xj, pf8f_load(t + 16), a2);
  }
  pf8f_store(y, a0); pf8f_store(y + 8, a1); pf8f_store(y + 16, a2);
}
/* scale y (20 non-negative floats, entries 20..23 zero) to a maximum in [1, 2): returns the power
 * of two taken out, or INT32_MIN if y is all zero (or not finite, or subnormal at the maximum) */
static inline int pf_norm(float *y)
{
  const pf8f a = pf8f_load(y), b = pf8f_load(y + 8), c = pf8f_load(y + 16);
  const float mx = pf8f_hmax(pf8f_max(pf8f_max(a, b), c));
  if (!(mx > 0.0f) || !(mx < INFINITY)) return INT32_MIN;
  uint32_t bits; memcpy(&bits, &mx, 4);
  const int fe = (int)(bits >> 23) & 0xFF;
  if (fe == 0) return INT32_MIN;                  /* subnormal maximum: leave it to f64 */
  const int e = fe - 127;                         /* mx = f * 2^e, f in [1, 2) */
  const uint32_t sb = (uint32_t)(127 - e) << 23;  /* 2^-e, a normal float for e in [-126, 126] */
  if (127 - e < 1 || 127 - e > 254) return INT32_MIN;
  float sc; memcpy(&sc, &sb, 4);
  const pf8f v = pf8f_set1(sc);                   /* a power of two: exact unless subnormal */
  pf8f_store(y, pf8f_mul(a, v)); pf8f_store(y + 8, pf8f_mul(b, v)); pf8f_store(y + 16, pf8f_mul(c, v));
  return e;
}
/* f64 row x (20, non-negative) to f32 with a maximum in [1, 2); entries 20..23 set to zero. The
 * maximum is exact in any order, the scaling a power of two, the conversion round to nearest. */
static inline int pf_in_scale(const double mx, double *sc)
{
  if (!(mx > 0.0) || !(mx < INFINITY)) return INT32_MIN;
  uint64_t bits; memcpy(&bits, &mx, 8);
  const int fe = (int)(bits >> 52) & 0x7FF;
  if (fe == 0) return INT32_MIN;
  const int e = fe - 1023;
  const uint64_t sb = (uint64_t)(1023 - e) << 52;
  memcpy(sc, &sb, 8);
  return e;
}
#if XP_AVX2
static inline int pf_in(float *y, const double *x)
{
  const __m256d v0 = _mm256_loadu_pd(x), v1 = _mm256_loadu_pd(x + 4), v2 = _mm256_loadu_pd(x + 8),
                v3 = _mm256_loadu_pd(x + 12), v4 = _mm256_loadu_pd(x + 16);
  __m256d m = _mm256_max_pd(_mm256_max_pd(_mm256_max_pd(v0, v1), _mm256_max_pd(v2, v3)), v4);
  __m128d h = _mm_max_pd(_mm256_castpd256_pd128(m), _mm256_extractf128_pd(m, 1));
  h = _mm_max_sd(h, _mm_unpackhi_pd(h, h));
  double sc;
  const int e = pf_in_scale(_mm_cvtsd_f64(h), &sc);
  if (e == INT32_MIN) return e;
  const __m256d v = _mm256_set1_pd(sc);
  _mm_store_ps(y, _mm256_cvtpd_ps(_mm256_mul_pd(v0, v))); _mm_store_ps(y + 4, _mm256_cvtpd_ps(_mm256_mul_pd(v1, v)));
  _mm_store_ps(y + 8, _mm256_cvtpd_ps(_mm256_mul_pd(v2, v))); _mm_store_ps(y + 12, _mm256_cvtpd_ps(_mm256_mul_pd(v3, v)));
  _mm_store_ps(y + 16, _mm256_cvtpd_ps(_mm256_mul_pd(v4, v)));
  _mm_store_ps(y + 20, _mm_setzero_ps());
  return e;
}
#else
static inline int pf_in(float *y, const double *x)
{
  float64x2_t v[10];
  for (int i = 0; i < 10; i++) v[i] = vld1q_f64(x + 2 * i);
  float64x2_t m = v[0];
  for (int i = 1; i < 10; i++) m = vmaxq_f64(m, v[i]);
  double sc;
  const int e = pf_in_scale(vmaxvq_f64(m), &sc);
  if (e == INT32_MIN) return e;
  const float64x2_t vs = vdupq_n_f64(sc);
  for (int i = 0; i < 10; i += 2)
    vst1q_f32(y + 2 * i, vcombine_f32(vcvt_f32_f64(vmulq_f64(v[i], vs)), vcvt_f32_f64(vmulq_f64(v[i + 1], vs))));
  vst1q_f32(y + 20, vdupq_n_f32(0.0f));
  return e;
}
#endif
static float *pf_tensor(const double *T, int nrates)    /* f32 copy, rows padded to PF_W */
{
  float *o = (float *)aligned_alloc(32, sizeof(float) * (size_t)nrates * NS * PF_W);
  if (!o) return NULL;
  for (int r = 0; r < nrates; r++)
    for (int j = 0; j < NS; j++) {
      for (int i = 0; i < NS; i++) o[((size_t)r * NS + j) * PF_W + i] = (float)T[((size_t)r * NS + j) * NS + i];
      for (int i = NS; i < PF_W; i++) o[((size_t)r * NS + j) * PF_W + i] = 0.0f;
    }
  return o;
}

int xp_edge_scores8_f32src(double *out, long ldo, long col0, int ne, const int64_t *desc,
                           const double *statd, const double *qevtab, const int32_t *cat, int nsites,
                           int ncode, int nrates, const uint8_t *codes, int nq, const uint32_t *present)
{
  if (ne < 1 || ne > XP_EB || ncode < 1 || ncode > 31) return -1;
  const int NC = ncode + 1, BS = xp_pf_bs32;
  const double *D[XP_EB], *P_[XP_EB], *SB[XP_EB][8];
  const int64_t *DE[XP_EB], *PE[XP_EB], *SBE[XP_EB][8];
  int pmode[XP_EB], nsib[XP_EB];
  float *TH[XP_EB], *TP[XP_EB], *TS[XP_EB][8];
  int nalloc = 0; float *allocs[XP_EB * 10];
  int fail = 0;
  for (int lane = 0; lane < ne; lane++) {
    const int64_t *dd = desc + (size_t)lane * 32;
    D[lane] = (const double *)(intptr_t)dd[0]; DE[lane] = (const int64_t *)(intptr_t)dd[1];
    pmode[lane] = (int)dd[2];
    P_[lane] = (const double *)(intptr_t)dd[3]; PE[lane] = (const int64_t *)(intptr_t)dd[4];
    nsib[lane] = (int)dd[7];
    if (nsib[lane] > 8) { fail = 1; break; }
    allocs[nalloc++] = TH[lane] = pf_tensor((const double *)(intptr_t)dd[5], nrates);
    TP[lane] = NULL;
    if (pmode[lane] == 1) {
      allocs[nalloc++] = TP[lane] = pf_tensor((const double *)(intptr_t)dd[6], nrates);
      for (int k = 0; k < nsib[lane]; k++) {
        SB[lane][k] = (const double *)(intptr_t)dd[8 + k]; SBE[lane][k] = (const int64_t *)(intptr_t)dd[16 + k];
        allocs[nalloc++] = TS[lane][k] = pf_tensor((const double *)(intptr_t)dd[24 + k], nrates);
      }
    }
  }
  float *q32 = (float *)malloc(sizeof(float) * (size_t)nrates * ncode * NS);
  float *sd32 = (float *)aligned_alloc(32, sizeof(float) * PF_W);
  /* table row k: its 8 mantissas Tm[2k] and exponents Te[2k] side by side, one 64-byte line */
  pf8f *Tm = (pf8f *)aligned_alloc(64, sizeof(pf8f) * 2 * (size_t)BS * NC);
  pf8i *Te = (pf8i *)(Tm + 1);
  pf8f *Pq = (pf8f *)aligned_alloc(32, sizeof(pf8f) * ((size_t)nq + 1));
  int64_t *Eq = (int64_t *)calloc((size_t)nq * XP_EB + 8, sizeof(int64_t));
  int ok_alloc = q32 && sd32 && Tm && Te && Pq && Eq;
  for (int i = 0; i < nalloc; i++) ok_alloc = ok_alloc && allocs[i];
  if (fail || !ok_alloc) {
    for (int i = 0; i < nalloc; i++) free(allocs[i]);
    free(q32); free(sd32); free(Tm); free(Pq); free(Eq);
    return fail ? 1 : -1;
  }
  for (size_t i = 0; i < (size_t)nrates * ncode * NS; i++) q32[i] = (float)qevtab[i];
  for (int i = 0; i < PF_W; i++) sd32[i] = i < NS ? (float)statd[i] : 0.0f;
  const pf8f one = pf8f_set1(1.0f), zero = pf8f_zero(), tiny = pf8f_set1(0x1p-90f);
  const pf8i mant = pf8i_set1(0x007FFFFF), onebits = pf8i_set1(0x3F800000),
                bias = pf8i_set1(127);
  for (int q = 0; q < nq; q++) Pq[q] = one;
  float snT[NS][XP_EB] __attribute__((aligned(32)));
  float xa[PF_W] __attribute__((aligned(32))), xb[PF_W] __attribute__((aligned(32)));
  float dv[PF_W] __attribute__((aligned(32))), pv[PF_W] __attribute__((aligned(32)));
  int rc = 0;
  for (int s0 = 0; s0 < nsites && !rc; s0 += BS) {
    const int s1 = s0 + BS < nsites ? s0 + BS : nsites;
    for (int s = s0; s < s1; s++) {
      const int j = s - s0, c_s = cat[s];
      int32_t sev[XP_EB];
      for (int lane = 0; lane < XP_EB; lane++) {
        if (lane >= ne) { for (int st = 0; st < NS; st++) snT[st][lane] = 1.0f; sev[lane] = 0; continue; }
        const size_t to = (size_t)c_s * NS * PF_W;
        long long ex = 0; int dead = 0, e;
        /* distal side */
        e = pf_in(xa, D[lane] + (size_t)s * NS); if (e == INT32_MIN) dead = 1;
        ex += e + (DE[lane] ? DE[lane][s] : 0);
        pf_ev(dv, xa, TH[lane] + to);
        /* proximal side */
        if (pmode[lane] == 0) {
          e = pf_in(xb, P_[lane] + (size_t)s * NS); if (e == INT32_MIN) dead = 1;
          ex += e + (PE[lane] ? PE[lane][s] : 0);
        } else {
          e = pf_in(xa, P_[lane] + (size_t)s * NS); if (e == INT32_MIN) dead = 1;
          ex += e + (PE[lane] ? PE[lane][s] : 0);
          pf_ev(xb, xa, TP[lane] + to);
          for (int k = 0; k < nsib[lane]; k++) {
            e = pf_in(xa, SB[lane][k] + (size_t)s * NS); if (e == INT32_MIN) dead = 1;
            ex += e + (SBE[lane][k] ? SBE[lane][k][s] : 0);
            float t[PF_W] __attribute__((aligned(32)));
            pf_ev(t, xa, TS[lane][k] + to);
            for (int i = 0; i < PF_W; i += 8) pf8f_store(xb + i, pf8f_mul(pf8f_load(xb + i), pf8f_load(t + i)));
            e = pf_norm(xb); if (e == INT32_MIN) dead = 1; ex += e;
          }
          if (nsib[lane] == 0) { e = pf_norm(xb); if (e == INT32_MIN) dead = 1; ex += e; }
        }
        pf_ev(pv, xb, TH[lane] + to);
        for (int i = 0; i < PF_W; i += 8)
          pf8f_store(xa + i, pf8f_mul(pf8f_mul(pf8f_load(sd32 + i), pf8f_load(pv + i)),
                                                pf8f_load(dv + i)));
        e = pf_norm(xa); if (e == INT32_MIN) dead = 1; ex += e;
        if (ex > (1 << 23) || ex < -(1 << 23)) { rc = 1; break; }
        for (int st = 0; st < NS; st++) snT[st][lane] = dead ? 0.0f : xa[st];
        sev[lane] = (int32_t)ex;
      }
      if (rc) break;
      const pf8i sv = pf8i_loadu(sev);
      const float *qt = q32 + (size_t)c_s * ncode * NS;
      uint32_t pm = present[s];
      while (pm) {
        const int c = __builtin_ctz(pm); pm &= pm - 1;
        const float *qv = qt + (size_t)c * NS;
        pf8f acc = zero;
        for (int st = 0; st < NS; st++)
          acc = pf8f_fmadd(pf8f_set1(qv[st]), pf8f_load(snT[st]), acc);
        const pf8f ok = pf8f_cmpge(acc, tiny);
        const pf8i b = pf8f_as_i(acc);
        const pf8i ee = pf8i_sub(pf8i_srl23(b), bias);
        const pf8f m = pf8i_as_f(pf8i_or(pf8i_and(b, mant), onebits));
        Tm[2 * ((size_t)j * NC + c)] = pf8f_blend(zero, m, ok);
        Te[2 * ((size_t)j * NC + c)] = pf8i_add(ee, sv);
      }
      Tm[2 * ((size_t)j * NC + ncode)] = one; Te[2 * ((size_t)j * NC + ncode)] = pf8i_zero();
    }
    if (rc) break;
    pf32_queries(Tm, Te, NC, s0, s1, codes, nsites, nq, Pq, Eq);
  }
  if (!rc)
    for (int q = 0; q < nq; q++) {
      float p8[XP_EB]; pf8f_storeu(p8, Pq[q]);
      for (int lane = 0; lane < ne; lane++)
        out[(size_t)q * ldo + col0 + lane] = p8[lane] == 0.0f ? NAN
          : log((double)p8[lane]) + M_LN2 * (double)Eq[(size_t)q * XP_EB + lane];
    }
  for (int i = 0; i < nalloc; i++) free(allocs[i]);
  free(q32); free(sd32); free(Tm); free(Pq); free(Eq);
  return rc;
}
#else
int xp_edge_scores8_f32src(double *out, long ldo, long col0, int ne, const int64_t *desc,
                           const double *statd, const double *qevtab, const int32_t *cat, int nsites,
                           int ncode, int nrates, const uint8_t *codes, int nq, const uint32_t *present)
{ (void)out; (void)ldo; (void)col0; (void)ne; (void)desc; (void)statd; (void)qevtab; (void)cat;
  (void)nsites; (void)ncode; (void)nrates; (void)codes; (void)nq; (void)present; return 1; }
#endif

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

/* Mixture three-taxon tree (--full_model mixtures): xp_gtt_new with category r's own eigen-system
 * (mx, ml, mxit) and frequencies (mstatd). */
TT *xp_mtt_new(int m, int R, const double *mstatd, const double *mx, const double *ml, const double *mxit,
               const double *rates, const double *query_masked)
{
  TT *t = xp_gtt_new(m, R, mstatd, mx, ml, mxit, rates, query_masked);
  if (!t) return NULL;
  t->mix = 1; t->mx = mx; t->ml = ml; t->mxit = mxit; t->mstatd = mstatd;
  return t;
}

/* ten_statd_pairwise_prod with frequencies per rate category (mixtures): mstatd[r][st] * a * b */
void xp_m_statd_pairwise(double *dst, int64_t *de, const double *mstatd, const double *a,
                         const int64_t *ae, const double *b, const int64_t *be, int nsites, int R)
{
  for (int s = 0; s < nsites; s++)
    for (int r = 0; r < R; r++) {
      const double *sd = mstatd + (size_t)r * NS;
      for (int st = 0; st < NS; st++) {
        *dst = sd[st] * (*a) * (*b);
        dst++; a++; b++;
      }
    }
  for (int s = 0; s < nsites; s++) de[s] = ae[s] + be[s];
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
  tt_tensor(t, t->tensorD, distal);
  t->prox_bl = cut_bl - distal;
  tt_tensor(t, t->tensorP, cut_bl - distal);
  tt_dist_pass(t);
  t->assert_fail = 0;
}
