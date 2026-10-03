//! Link the third-party C the engine needs for bit-exactness (see gplace/engine/build.py, which
//! compiles it and passes the paths through environment variables):
//!   XP_GSL_PREFIX  static, -fPIC GSL 1.16 (lib/libgsl.a, lib/libgslcblas.a)
//!   XP_DEPS_OBJS   CORE-MATH exp/log/pow/sin and glibc 2.19 hypot/log1p object files (space-separated)
//!   XP_LINK_ARGS   linker arguments binding GSL's libm calls (exp, log, pow, sin) to CORE-MATH and
//!                  the other libm functions GSL references to aborting traps (space-separated)
use std::env;

fn main() {
    for var in ["XP_GSL_PREFIX", "XP_DEPS_OBJS", "XP_LINK_ARGS"] {
        println!("cargo:rerun-if-env-changed={var}");
    }
    let gsl = env::var("XP_GSL_PREFIX").expect("XP_GSL_PREFIX (GSL 1.16 prefix) is not set");
    for obj in env::var("XP_DEPS_OBJS").unwrap_or_default().split_whitespace() {
        println!("cargo:rustc-link-arg={obj}");
    }
    println!("cargo:rustc-link-search=native={gsl}/lib");
    println!("cargo:rustc-link-lib=static=gsl");
    println!("cargo:rustc-link-lib=static=gslcblas");
    // bind the library's own references (GSL -> hypot, cr_*) inside the library
    println!("cargo:rustc-link-arg=-Wl,-Bsymbolic");
    for arg in env::var("XP_LINK_ARGS").unwrap_or_default().split_whitespace() {
        println!("cargo:rustc-link-arg={arg}");
    }
    println!("cargo:rustc-link-lib=m");
}
