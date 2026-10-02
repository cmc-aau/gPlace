//! The C functions the engine calls: GSL 1.16, CORE-MATH (correctly rounded exp/log/pow/sin, as
//! glibc 2.19's are) and the few libm functions that are exact everywhere (frexp, ldexp). All of
//! them are linked statically into the engine (build.rs); none comes from the system libm except
//! frexp/ldexp (and sqrt, which Rust compiles to an instruction).

use std::os::raw::{c_int, c_void};

/// `gsl_function`: the objective of the Brent minimiser.
#[repr(C)]
pub struct GslFunction {
    pub function: extern "C" fn(f64, *mut c_void) -> f64,
    pub params: *mut c_void,
}

/// `gsl_matrix` (only `data` is used: matrices are allocated n x n, so tda == n).
#[repr(C)]
pub struct GslMatrix {
    pub size1: usize,
    pub size2: usize,
    pub tda: usize,
    pub data: *mut f64,
    pub block: *mut c_void,
    pub owner: c_int,
}

/// `gsl_vector`
#[repr(C)]
pub struct GslVector {
    pub size: usize,
    pub stride: usize,
    pub data: *mut f64,
    pub block: *mut c_void,
    pub owner: c_int,
}

/// Opaque GSL types.
#[repr(C)]
pub struct GslMinFminimizer {
    _private: [u8; 0],
}
#[repr(C)]
pub struct GslMinFminimizerType {
    _private: [u8; 0],
}
#[repr(C)]
pub struct GslEigenSymmvWorkspace {
    _private: [u8; 0],
}

extern "C" {
    // CORE-MATH (gplace/engine/csrc/core-math)
    pub fn cr_exp(x: f64) -> f64;
    pub fn cr_log(x: f64) -> f64;
    pub fn cr_pow(x: f64, y: f64) -> f64;
    pub fn cr_sin(x: f64) -> f64;
    // glibc 2.19's e_hypot.c (gplace/engine/csrc/glibc), compiled as `hypot`
    pub fn hypot(x: f64, y: f64) -> f64;
    // libm, exact on every platform
    pub fn frexp(x: f64, exp: *mut c_int) -> f64;
    pub fn ldexp(x: f64, exp: c_int) -> f64;

    // GSL 1.16
    pub fn gsl_set_error_handler_off() -> *mut c_void;
    pub fn gsl_matrix_alloc(n1: usize, n2: usize) -> *mut GslMatrix;
    pub fn gsl_matrix_free(m: *mut GslMatrix);
    pub fn gsl_vector_alloc(n: usize) -> *mut GslVector;
    pub fn gsl_vector_free(v: *mut GslVector);
    pub fn gsl_eigen_symmv_alloc(n: usize) -> *mut GslEigenSymmvWorkspace;
    pub fn gsl_eigen_symmv_free(w: *mut GslEigenSymmvWorkspace);
    pub fn gsl_eigen_symmv(a: *mut GslMatrix, eval: *mut GslVector, evec: *mut GslMatrix,
                           w: *mut GslEigenSymmvWorkspace) -> c_int;
    pub static gsl_min_fminimizer_brent: *const GslMinFminimizerType;
    pub fn gsl_min_fminimizer_alloc(t: *const GslMinFminimizerType) -> *mut GslMinFminimizer;
    pub fn gsl_min_fminimizer_free(s: *mut GslMinFminimizer);
    pub fn gsl_min_fminimizer_set(s: *mut GslMinFminimizer, f: *mut GslFunction, x_minimum: f64,
                                  x_lower: f64, x_upper: f64) -> c_int;
    pub fn gsl_min_fminimizer_iterate(s: *mut GslMinFminimizer) -> c_int;
    pub fn gsl_min_fminimizer_x_minimum(s: *const GslMinFminimizer) -> f64;
    pub fn gsl_min_fminimizer_x_lower(s: *const GslMinFminimizer) -> f64;
    pub fn gsl_min_fminimizer_x_upper(s: *const GslMinFminimizer) -> f64;
    pub fn gsl_cdf_gamma_P(x: f64, a: f64, b: f64) -> f64;
    pub fn gsl_cdf_gamma_Pinv(p: f64, a: f64, b: f64) -> f64;
    pub fn gsl_sf_gamma_inc_P(a: f64, x: f64) -> f64;
}

/// Correctly rounded natural exponential (as glibc 2.19's, which pplacer's binary contains).
#[inline(always)]
pub fn exp(x: f64) -> f64 {
    unsafe { cr_exp(x) }
}

/// Correctly rounded natural logarithm.
#[inline(always)]
pub fn log(x: f64) -> f64 {
    unsafe { cr_log(x) }
}

/// C's frexp: the binary exponent e with x = m * 2^e, 0.5 <= |m| < 1 (0 for x = 0).
#[inline(always)]
pub fn frexp_exponent(x: f64) -> i32 {
    let mut e: c_int = 0;
    unsafe { frexp(x, &mut e) };
    e
}

/// C's ldexp(x, e) = x * 2^e.
#[inline(always)]
pub fn ldexp_(x: f64, e: i32) -> f64 {
    unsafe { ldexp(x, e) }
}
