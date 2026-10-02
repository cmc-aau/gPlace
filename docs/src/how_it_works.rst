.. _how_it_works:

How gPlace works
================

What pplacer computes
---------------------

For each query, pplacer v1.1.alpha19 (``pplacer -m WAG -c <refpkg>``, as run by GTDB-Tk):

1. **Pre-masks** the alignment: keeps columns informative in at least one reference and one query.
2. **Builds the model**: WAG exchangeabilities and frequencies, FastTree's 20 CAT rates and
   per-site categories from the refpkg; one eigendecomposition gives P(t) for every rate.
3. **Caches likelihood vectors**: for every node, a *distal* vector (the tree below it,
   post-order) and a *proximal* vector (the rest of the tree, pre-order), each *sites x 20*
   doubles with per-site power-of-two exponents, plus a *supernode* vector for the prefilter.
4. **Prefilters**: scores the query against every edge (query evolved by 0.1, dotted with the
   supernode) and ranks the edges.
5. **Baseball**: fully optimises edges in ranked order (pendant and distal branch lengths,
   alternating GSL Brent, tolerance 0.01) until 6 edges are more than 3 log-units worse than the
   best so far, or 40 edges have been tried.
6. **Keeps and refines**: up to 7 edges within log(0.01) of the best are re-optimised at 1e-5;
   the like-weight ratios are normalised likelihoods.

What gPlace changes
-------------------

Only what is stored and when it is computed:

* distal/proximal vectors are stored for **internal nodes only**, before pplacer's final
  "pull exponents (-10)" pass, which is re-applied to a copy when a vector is used. Rescaling by
  powers of two is exact, so every value is bit-identical;
* **leaf vectors are rebuilt on demand**: a leaf's distal vector is its sequence, and its
  proximal vector is its parent's proximal vector times its siblings' distal vectors, in
  pplacer's operation order;
* **supernodes are computed during the prefilter**, eight edges at a time, and discarded;
* all queries of a reference package share one cache; a query's result does not depend on the
  other queries.

The memory needed drops from *3 x nodes x sites x 20 x 8 bytes* to about
*2 x internal nodes x sites x 20 x 8 bytes* plus the prefilter scores.

Tiled placement
---------------

For trees whose vectors still do not fit (``--tile_edges``, ``--max_memory``), the tree is cut
into connected **tiles** of at most N edges. A tile touches the rest of the tree only through
the edge above its root and the roots of the tiles below it, and everything it needs from the
outside passes through those points:

* the *distal* vector of a tile root depends only on the tree below it: a bottom-up pass builds
  the tiles one at a time and stores the distal vector of each tile root;
* the *proximal* vector of a tile root depends only on the tree outside it: a top-down pass
  stores the proximal vector of each tile root.

With those boundary vectors (two per tile), a tile's vectors are rebuilt with the same
operations, operands and order as the untiled cache, so every value is bit-identical, and
every edge's likelihood is a full-tree likelihood, comparable across tiles. Only one tile is in
memory at a time. The cut: walking up the tree, the part of a node's subtree not yet in a tile
stays whole while it has at most N nodes; beyond that, its internal children with the largest
parts become tile roots (leaves never root a tile).

pplacer's algorithm is unchanged; only the order of the work is:

1. **Prefilter**, tile by tile and in batches of 20,000 queries: each query keeps its 160 best
   edges in pplacer's ranking order (score descending, ties by ascending edge number), which is
   pplacer's global ranking without a *queries x edges* score table.
2. **Baseball in rounds**: each first-stage optimisation starts from pplacer's default branch
   lengths, so its result does not depend on the edges tried before. The next 8 ranked edges of
   every query still searching are optimised, grouped by tile, and pplacer's stopping rule is
   replayed on the stored results; a few extra edges are evaluated and discarded. pplacer's
   search needs at most 41 successful edges; a query with many failed optimisations that runs
   out of candidates is re-ranked with more.
3. **Refinement** of the kept edges at 1e-5, tile by tile, then the like-weight ratios.

Each pass rebuilds the tiles, so tiled placement is 1.7-2.5x slower than untiled placement on a
tree that fits in memory; ``--max_memory`` therefore places untiled whenever it can.

Platforms
---------

The engine is built in variants from the same C source: ``avx2`` (x86-64; the kernels process 4
*different* results per instruction), ``neon`` (ARM64; 2 per instruction) and ``portable``
(scalar loops performing the same operations in the same order; any 64-bit CPU). gPlace loads
``avx2`` on x86-64 CPUs that have AVX2 and FMA, ``neon`` on ARM64 and ``portable`` elsewhere. The
SIMD kernels never fuse a multiply and an add: each lane does exactly what the scalar loop does. Neither depends on the CPU's floating-point details beyond IEEE
754 double precision: no fused multiply-add is contracted implicitly, and the only libm
functions taken from the system are exact by definition (``sqrt``, ``frexp``, ``ldexp``, and
``fma`` for the portable engine on x86-64).

Exactness rules
---------------

Every floating-point result keeps pplacer's operands, operation order and rounding:

* the engine (Rust, ``gplace/engine/rust``) mirrors pplacer's ``linear_c.c`` loops and GSL's
  ``dgemv`` loop order; Rust never fuses a multiply and an add, and SIMD lanes
  only process *different* results, no reduction is reordered, and no multiply-add is fused
  (``-ffp-contract=off``);
* the numerical libraries match pplacer's statically linked binary: GSL 1.16 (GSL 2.x changes the
  eigendecomposition in its last bits), glibc 2.19's ``hypot`` (glibc >= 2.35 changed it) and
  correctly-rounded ``exp``/``log``/``pow``/``sin`` (CORE-MATH gives the same results as glibc
  2.19's; GSL's own calls in the gamma functions are bound to them as well, and no inexact system
  libm function is imported);
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
