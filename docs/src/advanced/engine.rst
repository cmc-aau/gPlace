.. _advanced/engine:

Engine internals
================

How gPlace's engine runs on different CPUs, the rules that keep every result identical to
pplacer's, and what makes it fast. :ref:`how_it_works` describes what the engine computes.

Platforms
---------

The engine is built in variants from the same C source: ``avx2`` (x86-64; the kernels process 4
*different* results per instruction), ``neon`` (ARM64; 2 per instruction) and ``portable``
(scalar loops performing the same operations in the same order; any 64-bit CPU). gPlace loads
``avx2`` on x86-64 CPUs that have AVX2 and FMA, ``neon`` on ARM64 and ``portable`` elsewhere. The
SIMD kernels never fuse a multiply and an add: each lane does exactly what the scalar loop does. Neither depends on the CPU's floating-point details beyond IEEE
754 double precision: no fused multiply-add is contracted implicitly, and the only libm
functions taken from the system are exact by definition (``sqrt``, ``frexp``, ``ldexp``, and
``fma`` for the portable engine on x86-64). On macOS
the engine is built the same way with Apple's clang; Apple's linker has no ``--defsym``, so GSL's
libm calls are bound to CORE-MATH and the traps by definitions inside the engine
(``csrc/libm_bind.c``) instead.

.. figure:: /_static/figures/engines.svg
   :width: 100%
   :alt: Three rows of 16 boxes: avx2 fills them four at a time, neon two at a time and portable one at a time; all three end with identical bits.

   Each AVX2 instruction computes 4 different results and each NEON instruction 2, every one with the
   operations of the portable loop in the same order (no fused multiply-add, no reordered sum): only the
   speed differs.

Exactness rules
---------------

Every floating-point result keeps pplacer's operands, operation order and rounding:

* the engine (C, ``gplace/engine/csrc/xp.c``) mirrors pplacer's ``linear_c.c`` loops and GSL's
  ``dgemv`` loop order; SIMD lanes only process *different* results, no reduction is reordered,
  and no multiply-add is fused (``-ffp-contract=off``);
* the numerical libraries match pplacer's statically linked binary: GSL 1.16 (GSL 2.x changes the
  eigendecomposition in its last bits), glibc 2.19's ``hypot`` (glibc >= 2.35 changed it) and
  ``log1p`` (GSL's inverse gamma CDF calls it for 21 or more gamma categories), and
  correctly-rounded ``exp``/``log``/``pow``/``sin`` (CORE-MATH gives the same results as glibc
  2.19's; GSL's own calls in the gamma functions are bound to them as well, and no inexact system
  libm function is imported);
* GSL errors behave as in pplacer, which turns every GSL error into an exception: the optimiser's
  are handled as pplacer handles them, and a GSL error while the discrete gamma rates are computed
  (e.g. alpha 10^6: "inverse failed to converge") stops gPlace with pplacer's message;
* pplacer's control-flow quirks are reproduced: the optimiser keeps the *last evaluated* branch
  lengths; the bisection start finder evaluates f(right) before f(left); an assertion failing
  inside a GSL callback returns NaN (ocamlgsl's ``caml_callback_exn``), so GSL reports an error
  and pplacer retries from the default branch lengths; ties in the ranking go to the lower edge
  number;
* the jplace file prints numbers with ``%.12g`` and writes the tree exactly as
  ``Newick_bark.to_newick_string`` does.

Speed
-----

The cache is built by a thread pool that starts a node as soon as its children (distal) or its
parent (proximal) are done. The prefilter and the placements run in forked worker processes; the
large arrays live in shared memory, so they are neither copied nor re-mapped per process.
Placement uses a memory-lean three-taxon kernel: the query's evolved vectors come from a small
table of its distinct characters, and ``(pi x proximal) x distal`` is formed once per branch-length
change, with pplacer's association.
