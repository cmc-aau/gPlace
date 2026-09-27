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

Exactness rules
---------------

Every floating-point result keeps pplacer's operands, operation order and rounding:

* the C kernels mirror pplacer's ``linear_c.c`` loops and GSL's ``dgemv`` loop order; SIMD lanes
  only process *different* results, no reduction is reordered, and no multiply-add is fused
  (``-ffp-contract=off``);
* the numerical libraries match pplacer's statically linked binary: GSL 1.16 (GSL 2.x changes the
  eigendecomposition in its last bits), glibc 2.19's ``hypot`` (glibc >= 2.35 changed it) and
  correctly-rounded ``exp``/``log`` (CORE-MATH gives the same results as glibc 2.19's);
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
