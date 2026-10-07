/*
 * macOS only (gplace/engine/build.py): binds GSL's libm calls inside the engine.
 *
 * GSL's own calls to exp, log, pow and sin must give glibc 2.19's results, which are correctly
 * rounded (as CORE-MATH's are), and the other inexact libm functions GSL references must never be
 * taken from the system. On Linux the linker does this with --defsym (build.py LIBM_MAP); Apple's
 * linker has no --defsym, so on macOS the engine defines those functions itself. A definition in
 * the link wins over the system library for every reference in the link, so GSL's calls bind
 * here; the symbols are hidden, so nothing outside the engine sees them (other libraries keep the
 * system's functions). hypot and log1p are glibc 2.19's own code (csrc/glibc), linked the same way.
 * The traps abort: gPlace never calls them, and if an input ever made GSL call one, a silently
 * inexact result would break bit-exactness. Compiled with -fno-builtin.
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 */
#include <stdio.h>
#include <stdlib.h>

#define XP_HIDDEN __attribute__((visibility("hidden")))

double cr_exp(double); double cr_log(double); double cr_pow(double, double); double cr_sin(double);

XP_HIDDEN double exp(double x) { return cr_exp(x); }
XP_HIDDEN double log(double x) { return cr_log(x); }
XP_HIDDEN double pow(double x, double y) { return cr_pow(x, y); }
XP_HIDDEN double sin(double x) { return cr_sin(x); }

static void xp_trap(const char *f)
{
  fprintf(stderr, "gPlace engine: unexpected call to libm %s; results would not be bit-exact. "
          "Please report this.\n", f);
  abort();
}

/* build.py LIBM_TRAPS, and the functions clang may call instead of them on macOS (sin and cos of
 * the same argument -> __sincos_stret, pow(10, x) -> __exp10) */
#define XP_TRAP1(f) XP_HIDDEN double f(double x) { (void)x; xp_trap(#f); return 0.; }
#define XP_TRAP2(f) XP_HIDDEN double f(double x, double y) { (void)x; (void)y; xp_trap(#f); return 0.; }
XP_TRAP1(cos) XP_TRAP1(cosh) XP_TRAP1(sinh) XP_TRAP1(tan) XP_TRAP1(tanh) XP_TRAP1(expm1)
XP_TRAP1(atan) XP_TRAP1(asin) XP_TRAP1(acos) XP_TRAP1(acosh) XP_TRAP1(atanh) XP_TRAP1(asinh)
XP_TRAP1(log10) XP_TRAP1(log2) XP_TRAP1(exp2) XP_TRAP1(cbrt) XP_TRAP1(erf) XP_TRAP1(erfc)
XP_TRAP1(lgamma) XP_TRAP1(tgamma) XP_TRAP1(__exp10)
XP_TRAP2(atan2)
XP_HIDDEN void sincos(double x, double *s, double *c) { (void)x; (void)s; (void)c; xp_trap("sincos"); }
XP_HIDDEN void __sincos_stret(void) { xp_trap("__sincos_stret"); }
XP_HIDDEN void __sincospi_stret(void) { xp_trap("__sincospi_stret"); }
/* The avx2 engine computes fma() with the FMA instruction; any call would come from elsewhere.
 * The x86-64 portable engine calls the system's fma, which is correctly rounded by definition. */
#ifdef XP_TRAP_FMA
XP_HIDDEN double fma(double x, double y, double z) { (void)x; (void)y; (void)z; xp_trap("fma"); return 0.; }
#endif
