/* minimal shim for compiling glibc 2.19 dbl-64 e_hypot.c and s_log1p.c standalone */
#include <stdint.h>
typedef uint32_t u_int32_t;
typedef union { double value; struct { uint32_t lsw, msw; } p; } ieee_double_shape_type;
#define GET_HIGH_WORD(i,d) do { ieee_double_shape_type u_; u_.value = (d); (i) = u_.p.msw; } while (0)
#define GET_LOW_WORD(i,d)  do { ieee_double_shape_type u_; u_.value = (d); (i) = u_.p.lsw; } while (0)
#define SET_HIGH_WORD(d,v) do { ieee_double_shape_type u_; u_.value = (d); u_.p.msw = (v); (d) = u_.value; } while (0)
#define SET_LOW_WORD(d,v)  do { ieee_double_shape_type u_; u_.value = (d); u_.p.lsw = (v); (d) = u_.value; } while (0)
#define INSERT_WORDS(d,ix0,ix1) do { ieee_double_shape_type u_; u_.p.msw = (ix0); u_.p.lsw = (ix1); (d) = u_.value; } while (0)
#define EXTRACT_WORDS(ix0,ix1,d) do { ieee_double_shape_type u_; u_.value = (d); (ix0) = u_.p.msw; (ix1) = u_.p.lsw; } while (0)
#include <math.h>
#define __ieee754_sqrt sqrt
#define strong_alias(a,b)
#define weak_alias(a,b)
/* sysdeps/generic/math_private.h: evaluate x (for its floating-point flags); no effect on results */
#define math_force_eval(x) ({ __typeof (x) __x = (x); __asm __volatile__ ("" : : "m" (__x)); })
