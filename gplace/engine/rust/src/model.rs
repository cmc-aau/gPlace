//! The substitution model: eigendecomposition (GSL 1.16), transition matrices
//! P(t) = X exp(rate * t * Lambda) X^-1 for every rate category (pplacer's Diagd.multi_exp,
//! linear_c.c dediagonalize), and the discrete gamma rates (Gamma.discrete_gamma).

use crate::ffi::{self, exp};
use crate::{simd, sl, sl_mut, NS, SIMD20};
use std::ffi::CStr;
use std::os::raw::{c_char, c_int};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Mutex;

/// Gsl_eigen.symmv: unsorted eigendecomposition of a symmetric n x n matrix (row-major).
#[no_mangle]
pub unsafe extern "C" fn xp_symmv(m_in: *const f64, eval: *mut f64, evec: *mut f64, n: c_int) -> c_int {
    let n = n as usize;
    let m = ffi::gsl_matrix_alloc(n, n);
    let ev = ffi::gsl_vector_alloc(n);
    let vec = ffi::gsl_matrix_alloc(n, n);
    let w = ffi::gsl_eigen_symmv_alloc(n);
    std::ptr::copy_nonoverlapping(m_in, (*m).data, n * n);
    let st = ffi::gsl_eigen_symmv(m, ev, vec, w);
    std::ptr::copy_nonoverlapping((*ev).data, eval, n);
    std::ptr::copy_nonoverlapping((*vec).data, evec, n * n);
    ffi::gsl_eigen_symmv_free(w);
    ffi::gsl_matrix_free(vec);
    ffi::gsl_vector_free(ev);
    ffi::gsl_matrix_free(m);
    st
}

/// Reference dediagonalize (linear_c.c), row-major: dst[i][j] = sum_k lambda[k] * u[i][k] * uit[j][k].
/// Kept for the self-tests.
#[no_mangle]
pub unsafe extern "C" fn xp_dediag_ref(dst: *mut f64, u: *const f64, lambda: *const f64, uit: *const f64,
                                       n: c_int) {
    let n = n as usize;
    let (dst, u, lambda, uit) = (sl_mut(dst, n * n), sl(u, n * n), sl(lambda, n), sl(uit, n * n));
    for i in 0..n {
        for j in 0..n {
            let mut d = 0.0;
            for k in 0..n {
                d += lambda[k] * u[i * n + k] * uit[j * n + k];
            }
            dst[i * n + j] = d;
        }
    }
}

/// Reference Diagd.multi_exp: tensor[r] = X exp(rates[r] * bl * Lambda) X^-1 (row-major) for the
/// rates with mask[r]. Kept for the self-tests.
#[no_mangle]
pub unsafe extern "C" fn xp_multi_exp_ref(tensor: *mut f64, x: *const f64, l: *const f64, xit: *const f64,
                                          rates: *const f64, nrates: c_int, mask: *const u8, bl: f64) {
    let nrates = nrates as usize;
    let (rates, mask, l) = (sl(rates, nrates), sl(mask, nrates), sl(l, NS));
    sl_mut(tensor, nrates * NS * NS).fill(0.0);
    let mut util = [0.0; NS];
    for r in 0..nrates {
        if mask[r] == 0 {
            continue;
        }
        for i in 0..NS {
            util[i] = exp(rates[r] * bl * l[i]);
        }
        xp_dediag_ref(tensor.add(r * NS * NS), x, util.as_ptr(), xit, NS as c_int);
    }
}

/// Dediagonalize with a TRANSPOSED result: dst_t[j][i] = P[i][j]. Each value is the scalar
/// expression of xp_dediag_ref; lambda[k] * u[i][k] does not depend on j, so the SIMD kernel
/// computes it once (the same value).
#[no_mangle]
pub unsafe extern "C" fn xp_dediag(dst_t: *mut f64, u: *const f64, lambda: *const f64, uit: *const f64,
                                   n: c_int) {
    let n = n as usize;
    if n != NS || !SIMD20 {
        let (dst_t, u, lambda, uit) = (sl_mut(dst_t, n * n), sl(u, n * n), sl(lambda, n), sl(uit, n * n));
        for i in 0..n {
            for j in 0..n {
                let mut d = 0.0;
                for k in 0..n {
                    d += lambda[k] * u[i * n + k] * uit[j * n + k];
                }
                dst_t[j * n + i] = d;
            }
        }
        return;
    }
    let mut a = [0.0; NS * NS]; // a[k][i] = lambda[k] * u[i][k]
    let (u, lambda) = (sl(u, NS * NS), sl(lambda, NS));
    for k in 0..NS {
        for i in 0..NS {
            a[k * NS + i] = lambda[k] * u[i * NS + k];
        }
    }
    simd::dediag20(dst_t, a.as_ptr(), uit);
}

/// Diagd.multi_exp with transposed matrices (the layout the evolve kernels read).
#[no_mangle]
pub unsafe extern "C" fn xp_multi_exp(tensor: *mut f64, x: *const f64, l: *const f64, xit: *const f64,
                                      rates: *const f64, nrates: c_int, mask: *const u8, bl: f64) {
    let nrates = nrates as usize;
    let (rates, mask, l) = (sl(rates, nrates), sl(mask, nrates), sl(l, NS));
    sl_mut(tensor, nrates * NS * NS).fill(0.0);
    let mut util = [0.0; NS];
    for r in 0..nrates {
        if mask[r] == 0 {
            continue;
        }
        for i in 0..NS {
            util[i] = exp(rates[r] * bl * l[i]);
        }
        xp_dediag(tensor.add(r * NS * NS), x, util.as_ptr(), xit, NS as c_int);
    }
}

/// pplacer turns every GSL error into an OCaml exception (`Gsl_error.init ()` in ppatteries.ml), and
/// its discrete gamma catches none: the first GSL error raised while the rates are computed ends the
/// run (e.g. `Gsl_exn(6, "inverse failed to converge")` for alpha 1e6). The engine records that first
/// error so that gPlace stops with the same message; outside the rates GSL errors are ignored (the
/// optimiser uses GSL's return codes, as pplacer's handlers of Gsl_exn do).
static RECORDING: AtomicBool = AtomicBool::new(false);
static FIRST_ERROR: Mutex<Option<(c_int, String)>> = Mutex::new(None);

/// The GSL error handler installed by xp_init.
pub(crate) unsafe extern "C" fn gsl_error_handler(reason: *const c_char, _file: *const c_char, _line: c_int,
                                                  gsl_errno: c_int) {
    if RECORDING.load(Ordering::Relaxed) {
        let mut first = FIRST_ERROR.lock().unwrap_or_else(|p| p.into_inner());
        if first.is_none() {
            let r = if reason.is_null() { String::new() } else { CStr::from_ptr(reason).to_string_lossy().into_owned() };
            *first = Some((gsl_errno, r));
        }
    }
}

/// The first GSL error of the last xp_discrete_gamma call: its gsl_errno (0 if none), with its message
/// copied into buf (NUL-terminated, truncated to len bytes).
#[no_mangle]
pub unsafe extern "C" fn xp_gamma_gsl_error(buf: *mut c_char, len: usize) -> c_int {
    let first = FIRST_ERROR.lock().unwrap_or_else(|p| p.into_inner());
    let Some((errno, reason)) = first.as_ref() else { return 0 };
    if len > 0 {
        let n = reason.len().min(len - 1);
        std::ptr::copy_nonoverlapping(reason.as_ptr() as *const c_char, buf, n);
        *buf.add(n) = 0;
    }
    *errno
}

/// Gamma.discrete_gamma: inverse of the gamma CDF, below a tiny epsilon -> 0.
fn inverse_gamma_cdf(alpha: f64, beta: f64, y: f64) -> f64 {
    let epsilon = 1e-7;
    let beta = 1. / beta;
    unsafe {
        if y < ffi::gsl_cdf_gamma_P(epsilon, alpha, beta) {
            return 0.;
        }
        ffi::gsl_cdf_gamma_Pinv(y, alpha, beta)
    }
}

/// Gamma.discrete_gamma n_cat alpha (Yang 1994, mean rates), then make_mean_one. Returns 0, -1 for
/// n_cat <= 0, or 2 if GSL raised an error (where pplacer stops; xp_gamma_gsl_error tells which).
#[no_mangle]
pub unsafe extern "C" fn xp_discrete_gamma(n_cat: c_int, alpha: f64, out: *mut f64) -> c_int {
    if n_cat <= 0 {
        return -1;
    }
    let n = n_cat as usize;
    let out = sl_mut(out, n);
    let mut freq = vec![0.0; n + 1];
    *FIRST_ERROR.lock().unwrap_or_else(|p| p.into_inner()) = None;
    RECORDING.store(true, Ordering::Relaxed);
    for (i, f) in freq.iter_mut().enumerate() {
        *f = if i == 0 {
            0.
        } else if i == n {
            1.
        } else {
            let cutpoint = inverse_gamma_cdf(alpha, 1., i as f64 / n as f64);
            ffi::gsl_sf_gamma_inc_P(alpha + 1., cutpoint)
        };
    }
    RECORDING.store(false, Ordering::Relaxed);
    if FIRST_ERROR.lock().unwrap_or_else(|p| p.into_inner()).is_some() {
        return 2;
    }
    let mut tot = 0.;
    for i in 0..n {
        out[i] = freq[i + 1] - freq[i];
        tot += out[i];
    }
    let f_n = n as f64;
    for o in out.iter_mut() {
        *o = f_n * *o / tot;
    }
    0
}
