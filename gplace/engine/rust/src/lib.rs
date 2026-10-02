//! gPlace placement engine (libgplace_core*.so).
//!
//! Exact re-implementation of the numerical core of pplacer v1.1.alpha19 (commit 807f6f3, the
//! build GTDB-Tk uses) for its Gcat (FastTree CAT) and Gmix (discrete gamma) models, protein and
//! DNA. Every floating-point result keeps pplacer's operands, operation order and rounding, so
//! placements are bit-identical:
//!
//! - dediagonalize / multi_exp                (pplacer_src/linear_c.c, diagd.ml)    -> `model`
//! - evolve_into = per-site gsl dgemv         (gcat_model.ml, gslcblas loop order)  -> `glv`
//! - pairwise / statd-pairwise products, perhaps_pull_exponent, mat_log_like3        -> `glv`
//! - the prefilter (h_ranking scores)                                                -> `prefilter`
//! - Minimization.brent + bisection start, Three_tax.optimize                        -> `three_tax`
//! - the gamma mixture model (gmix_model.ml, gamma.ml)                               -> all of them
//!
//! Rust never fuses a multiply and an add (no FMA contraction), so every `a * b + c` below is the
//! two roundings pplacer's C does. The SIMD kernels (`simd`) process 4 (AVX2) or 2 (NEON)
//! *different* results side by side, each with the same operations in the same order as the
//! scalar loop; no sum is reordered.
//!
//! The library is built per state count and kernel set (gplace/engine/build.py):
//! `--cfg xp_dna` gives 4 states (DNA), otherwise 20 (protein); `-C target-feature=+avx2,+fma`
//! gives the AVX2 kernels on x86-64; on ARM64 NEON is always available, and `--cfg
//! xp_force_scalar` builds the scalar ("portable") engine.
//!
//! This program is free software: you can redistribute it and/or modify it under the terms of the
//! GNU General Public License as published by the Free Software Foundation, either version 3 of
//! the License, or (at your option) any later version.

#![allow(clippy::missing_safety_doc, clippy::too_many_arguments)]

use std::os::raw::{c_char, c_int};

pub mod ffi;
pub mod glv;
pub mod model;
pub mod prefilter;
mod simd;
pub mod three_tax;

/// Number of character states: 20 amino acids, or 4 nucleotides with `--cfg xp_dna`.
#[cfg(xp_dna)]
pub const NS: usize = 4;
#[cfg(not(xp_dna))]
pub const NS: usize = 20;

/// Which kernels this build uses: 0 scalar, 1 AVX2, 2 NEON.
pub const SIMD_KIND: c_int = if simd::AVX2 { 1 } else if simd::NEON { 2 } else { 0 };

/// SIMD kernels written for 20 states (other state counts use the scalar forms).
pub const SIMD20: bool = NS == 20 && (simd::AVX2 || simd::NEON);

const fn kernels() -> &'static str {
    if simd::AVX2 { "AVX2" } else if simd::NEON { "NEON" } else { "portable, scalar" }
}

const fn arch() -> &'static str {
    if cfg!(target_arch = "x86_64") { "x86-64" } else if cfg!(target_arch = "aarch64") { "arm64" } else { "other" }
}

/// NUL-terminated version string: the C side of ctypes reads it as `char *`.
static VERSION: std::sync::OnceLock<String> = std::sync::OnceLock::new();

#[no_mangle]
pub extern "C" fn xp_engine_version() -> *const c_char {
    VERSION
        .get_or_init(|| {
            format!("gplace-core 3 ({NS} states): pplacer v1.1.alpha19 exact (GSL 1.16, glibc-2.19 hypot, \
                     CORE-MATH exp/log, {} {}; Rust)\0", arch(), kernels())
        })
        .as_ptr() as *const c_char
}

#[no_mangle]
pub extern "C" fn xp_simd() -> c_int {
    SIMD_KIND
}

#[no_mangle]
pub extern "C" fn xp_n_states() -> c_int {
    NS as c_int
}

/// GSL must report errors as return codes (the optimiser relies on them), not abort.
#[no_mangle]
pub extern "C" fn xp_init() -> c_int {
    unsafe { ffi::gsl_set_error_handler_off() };
    0
}

// ------------------------------------------------------------------------------- libm ---

#[no_mangle]
pub extern "C" fn xp_exp(v: f64) -> f64 {
    ffi::exp(v)
}

#[no_mangle]
pub extern "C" fn xp_log(v: f64) -> f64 {
    ffi::log(v)
}

/// pow and sin are only used by GSL (its gamma functions); exported for the self-tests.
#[no_mangle]
pub extern "C" fn xp_pow(a: f64, b: f64) -> f64 {
    unsafe { ffi::cr_pow(a, b) }
}

#[no_mangle]
pub extern "C" fn xp_sin(v: f64) -> f64 {
    unsafe { ffi::cr_sin(v) }
}

#[no_mangle]
pub extern "C" fn xp_hypot(a: f64, b: f64) -> f64 {
    unsafe { ffi::hypot(a, b) }
}

/// Targets of the libm functions GSL references but gPlace never calls (build.py LIBM_TRAPS):
/// abort instead of silently using an inexact system implementation.
macro_rules! libm_traps {
    ($($name:ident => $f:literal),* $(,)?) => {$(
        #[no_mangle]
        pub extern "C" fn $name() {
            eprintln!(concat!("gPlace engine: unexpected call to libm ", $f,
                              "; results would not be bit-exact. Please report this."));
            std::process::abort();
        }
    )*};
}

libm_traps!(
    xp_libm_trap_cos => "cos", xp_libm_trap_cosh => "cosh", xp_libm_trap_sinh => "sinh",
    xp_libm_trap_tan => "tan", xp_libm_trap_tanh => "tanh", xp_libm_trap_log1p => "log1p",
    xp_libm_trap_sincos => "sincos", xp_libm_trap_expm1 => "expm1", xp_libm_trap_atan => "atan",
    xp_libm_trap_atan2 => "atan2", xp_libm_trap_asin => "asin", xp_libm_trap_acos => "acos",
    xp_libm_trap_log10 => "log10", xp_libm_trap_log2 => "log2", xp_libm_trap_exp2 => "exp2",
    xp_libm_trap_cbrt => "cbrt", xp_libm_trap_erf => "erf", xp_libm_trap_erfc => "erfc",
    xp_libm_trap_lgamma => "lgamma", xp_libm_trap_tgamma => "tgamma", xp_libm_trap_fma => "fma",
    xp_libm_trap_acosh => "acosh", xp_libm_trap_atanh => "atanh", xp_libm_trap_asinh => "asinh",
);

/// Slice views of the arrays Python passes (numpy buffers; the caller guarantees the lengths).
#[inline(always)]
pub(crate) unsafe fn sl<'a, T>(p: *const T, n: usize) -> &'a [T] {
    if n == 0 { &[] } else { std::slice::from_raw_parts(p, n) }
}

#[inline(always)]
pub(crate) unsafe fn sl_mut<'a, T>(p: *mut T, n: usize) -> &'a mut [T] {
    if n == 0 { &mut [] } else { std::slice::from_raw_parts_mut(p, n) }
}
