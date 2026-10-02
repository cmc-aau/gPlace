//! Likelihood vectors ("glvs"): sites x NS doubles, plus one power-of-two exponent per site
//! (pplacer's Glv). CAT glvs are sites x NS; gamma (Gmix) glvs are sites x R x NS, site-major
//! (pplacer stores them rate-major; only the per-value operation order matters, and it is
//! pplacer's).

use crate::ffi::{frexp_exponent, ldexp_, log};
use crate::{simd, sl, sl_mut, NS, SIMD20};
use std::os::raw::c_int;
use std::sync::atomic::{AtomicBool, Ordering};

/// Gcat evolve_into for one site with a TRANSPOSED matrix pt[j][i] = P[i][j]:
/// y[i] = 0.0 + sum_j x[j] * P[i][j], j ascending from +0.0 (gslcblas dgemv NoTrans, alpha = 1,
/// beta = 0; 1.0 * temp == temp exactly).
#[inline(always)]
pub(crate) unsafe fn evolve_site(yv: *mut f64, xv: *const f64, pt: *const f64) {
    if SIMD20 {
        simd::evolve_site20(yv, xv, pt);
        return;
    }
    for i in 0..NS {
        let mut temp = 0.0;
        for j in 0..NS {
            temp += *xv.add(j) * *pt.add(j * NS + i);
        }
        *yv.add(i) = 0.0 + temp;
    }
}

/// Reference Gcat evolve_into with row-major matrices (self-tests).
#[no_mangle]
pub unsafe extern "C" fn xp_evolve_ref(dst: *mut f64, src: *const f64, nsites: c_int, cat: *const i32,
                                       tensor: *const f64) {
    let n = nsites as usize;
    let (dst, src, cat) = (sl_mut(dst, n * NS), sl(src, n * NS), sl(cat, n));
    for s in 0..n {
        let p = sl(tensor.add(cat[s] as usize * NS * NS), NS * NS);
        let xv = &src[s * NS..(s + 1) * NS];
        for (i, row) in p.chunks_exact(NS).enumerate() {
            let mut temp = 0.0;
            for (&x, &pij) in xv.iter().zip(row) {
                temp += x * pij;
            }
            dst[s * NS + i] = 0.0 + 1.0 * temp;
        }
    }
}

/// Gcat evolve_into for every site (transposed matrices from xp_multi_exp).
#[no_mangle]
pub unsafe extern "C" fn xp_evolve(dst: *mut f64, src: *const f64, nsites: c_int, cat: *const i32,
                                   tensor: *const f64) {
    for s in 0..nsites as usize {
        evolve_site(dst.add(s * NS), src.add(s * NS), tensor.add(*cat.add(s) as usize * NS * NS));
    }
}

/// dst = x .* y (elementwise), exponents summed.
#[no_mangle]
pub unsafe extern "C" fn xp_pairwise(dst: *mut f64, de: *mut i64, x: *const f64, xe: *const i64,
                                     y: *const f64, ye: *const i64, nsites: c_int) {
    let n = nsites as usize;
    let (dst, x, y) = (sl_mut(dst, n * NS), sl(x, n * NS), sl(y, n * NS));
    for i in 0..n * NS {
        dst[i] = x[i] * y[i];
    }
    sum_exponents(de, xe, ye, n);
}

#[inline(always)]
unsafe fn sum_exponents(de: *mut i64, xe: *const i64, ye: *const i64, n: usize) {
    let (de, xe, ye) = (sl_mut(de, n), sl(xe, n), sl(ye, n));
    for s in 0..n {
        de[s] = xe[s] + ye[s];
    }
}

/// dst = statd .* a .* b, i.e. (statd[st] * a) * b (mat_statd_pairwise_prod), exponents summed.
#[no_mangle]
pub unsafe extern "C" fn xp_statd_pairwise(dst: *mut f64, de: *mut i64, statd: *const f64, a: *const f64,
                                           ae: *const i64, b: *const f64, be: *const i64, nsites: c_int) {
    let n = nsites as usize;
    let (dst, statd, a, b) = (sl_mut(dst, n * NS), sl(statd, NS), sl(a, n * NS), sl(b, n * NS));
    for ((d_site, a_site), b_site) in dst.chunks_exact_mut(NS).zip(a.chunks_exact(NS)).zip(b.chunks_exact(NS)) {
        for (((d, &w), &x), &y) in d_site.iter_mut().zip(statd).zip(a_site).zip(b_site) {
            *d = w * x * y;
        }
    }
    sum_exponents(de, ae, be, n);
}

/// Glv.perhaps_pull_exponent: per site, if the frexp exponent of the largest entry is below
/// min_twoexp, scale the site by 2^-twoexp (ldexp(0.5, -twoexp + 1)) and add twoexp to its exponent.
#[no_mangle]
pub unsafe extern "C" fn xp_pull_exponent(a: *mut f64, e: *mut i64, nsites: c_int, min_twoexp: c_int) {
    let n = nsites as usize;
    let (a, e) = (sl_mut(a, n * NS), sl_mut(e, n));
    for s in 0..n {
        let v = &mut a[s * NS..(s + 1) * NS];
        let mut mx = v[0];
        for &x in &v[1..] {
            if x > mx {
                mx = x;
            }
        }
        let twoexp = frexp_exponent(mx);
        if twoexp < min_twoexp {
            let f = ldexp_(0.5, -twoexp + 1);
            for x in v.iter_mut() {
                *x *= f;
            }
            e[s] += twoexp as i64;
        }
    }
}

/// Sum of the exponents as a double, in site order.
#[no_mangle]
pub unsafe extern "C" fn xp_total_twoexp(e: *const i64, nsites: c_int) -> f64 {
    let mut tot = 0.;
    for &x in sl(e, nsites as usize) {
        tot += x as f64;
    }
    tot
}

/// mat_log_like3: sum over sites of log(sum_st ((statd * x) * y) * z).
#[no_mangle]
pub unsafe extern "C" fn xp_log_like3(statd: *const f64, x: *const f64, y: *const f64, z: *const f64,
                                      nsites: c_int) -> f64 {
    let n = nsites as usize;
    let (statd, x, y, z) = (sl(statd, NS), sl(x, n * NS), sl(y, n * NS), sl(z, n * NS));
    let mut ll_tot = 0.;
    for ((x_site, y_site), z_site) in x.chunks_exact(NS).zip(y.chunks_exact(NS)).zip(z.chunks_exact(NS)) {
        let mut util = 0.;
        for (((&w, &xv), &yv), &zv) in statd.iter().zip(x_site).zip(y_site).zip(z_site) {
            util += w * xv * yv * zv;
        }
        ll_tot += log(util);
    }
    ll_tot
}

// --------------------------------------------------------------- gamma mixture (Gmix) ---

/// ten_pairwise_prod + iba1_pairwise_sum.
#[no_mangle]
pub unsafe extern "C" fn xp_g_pairwise(dst: *mut f64, de: *mut i64, x: *const f64, xe: *const i64,
                                       y: *const f64, ye: *const i64, nsites: c_int, rr: c_int) {
    let n = nsites as usize;
    let m = n * rr as usize * NS;
    let (dst, x, y) = (sl_mut(dst, m), sl(x, m), sl(y, m));
    for i in 0..m {
        dst[i] = x[i] * y[i];
    }
    sum_exponents(de, xe, ye, n);
}

/// ten_statd_pairwise_prod: (statd[st] * a) * b for every site and rate.
#[no_mangle]
pub unsafe extern "C" fn xp_g_statd_pairwise(dst: *mut f64, de: *mut i64, statd: *const f64, a: *const f64,
                                             ae: *const i64, b: *const f64, be: *const i64, nsites: c_int,
                                             rr: c_int) {
    let n = nsites as usize;
    let k = n * rr as usize;
    let (dst, statd, a, b) = (sl_mut(dst, k * NS), sl(statd, NS), sl(a, k * NS), sl(b, k * NS));
    for ((d_site, a_site), b_site) in dst.chunks_exact_mut(NS).zip(a.chunks_exact(NS)).zip(b.chunks_exact(NS)) {
        for (((d, &w), &x), &y) in d_site.iter_mut().zip(statd).zip(a_site).zip(b_site) {
            *d = w * x * y;
        }
    }
    sum_exponents(de, ae, be, n);
}

/// Optional fix (NOT pplacer's behaviour; gplace --fix_gamma_underflow): frexp(0.0) gives exponent
/// 0, so in pplacer a rate whose values have all underflowed to 0.0 blocks the rescaling of the
/// whole site; on large trees the other rates then underflow too (site likelihood 0 -> -inf). With
/// the fix, all-zero rates are skipped when the exponent is chosen. Per process (forked workers
/// inherit it).
static GAMMA_FIX: AtomicBool = AtomicBool::new(false);

#[no_mangle]
pub extern "C" fn xp_set_gamma_fix(on: c_int) {
    GAMMA_FIX.store(on != 0, Ordering::Relaxed);
}

/// Gmix Glv.perhaps_pull_exponent: per site, the largest frexp exponent of the rates' maxima; if
/// below min_twoexp, every rate is scaled by ldexp(0.5, -max + 1).
#[no_mangle]
pub unsafe extern "C" fn xp_g_pull_exponent(a: *mut f64, e: *mut i64, nsites: c_int, rr: c_int, min_twoexp: c_int) {
    let (n, r) = (nsites as usize, rr as usize);
    let (a, e) = (sl_mut(a, n * r * NS), sl_mut(e, n));
    let fix = GAMMA_FIX.load(Ordering::Relaxed);
    for s in 0..n {
        let site = &mut a[s * r * NS..(s + 1) * r * NS];
        let mut max_twoexp = -i32::MAX;
        for v in site.chunks_exact(NS) {
            let mut mx = v[0]; // gsl_vector_max
            for &x in &v[1..] {
                if x > mx {
                    mx = x;
                }
            }
            if fix && mx == 0.0 {
                continue;
            }
            let twoexp = frexp_exponent(mx);
            if twoexp > max_twoexp {
                max_twoexp = twoexp;
            }
        }
        if max_twoexp == -i32::MAX {
            continue; // only with the fix: the whole site is zero
        }
        if max_twoexp < min_twoexp {
            let f = ldexp_(0.5, -max_twoexp + 1);
            for x in site.iter_mut() {
                *x *= f;
            }
            e[s] += max_twoexp as i64;
        }
    }
}
