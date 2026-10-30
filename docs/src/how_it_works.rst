.. _how_it_works:

How gPlace works
================

gPlace places each query, an aligned genome, on a fixed reference tree and writes the result as
a jplace file, as pplacer does. It runs pplacer's algorithm unchanged and changes only what is kept
in memory and when it is computed.

.. figure:: /_static/figures/overview.svg
   :width: 100%
   :alt: Three aligned query genomes, a reference tree with each query attached to one edge by a dashed orange branch, and the jplace file listing genome A on edge 3, genome B on edge 9 and genome C on edge 12, identical to pplacer.

   What placement produces. Each query is attached to the edge where it fits best (the orange dots);
   the jplace file lists, for every query, the edges kept with their likelihoods and branch lengths.

What pplacer computes
---------------------

.. figure:: /_static/figures/pipeline.svg
   :width: 100%
   :alt: Six steps in two rows: pre-mask, model and cache once per run; prefilter, baseball and keep and refine for each query. gPlace changes the cache (internal nodes only) and the prefilter (supernodes on the fly).

   pplacer's six steps. Steps 1 to 3 run once per reference package and steps 4 to 6 for each
   query. gPlace changes only the cache (step 3) and the prefilter (step 4), shown in green.

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

.. figure:: /_static/figures/baseball.svg
   :width: 100%
   :alt: Log-likelihoods of edges in prefilter order with a strike box 3 units below the best; six edges below it are strikes and the search stops at rank 11; seven edges are kept and their like-weight ratios are shown as bars.

   Steps 5 and 6 for one query (illustrative values). A strike is an edge more than 3 log-likelihood
   units below the best so far (the strike box follows the best); the search stops at the 6th strike, or
   after 40 edges, and the edges after it are never optimised. At most 7 edges within a factor 0.01 of
   the best likelihood are kept, refined and normalised into like-weight ratios.

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

.. figure:: /_static/figures/storage.svg
   :width: 100%
   :alt: The same six-leaf tree twice: pplacer keeps distal, proximal and supernode vectors for all 11 nodes (33 vectors); gPlace keeps distal and proximal vectors for the 5 internal nodes (10 vectors) and rebuilds leaf vectors when needed.

   What is kept for a six-leaf tree. pplacer keeps three vectors for every node; gPlace keeps two for
   each internal node, rebuilds a leaf's vectors when an edge next to it is used, and computes
   supernodes during the prefilter.

The memory needed drops from *3 x nodes x sites x 20 x 8 bytes* to about
*2 x internal nodes x sites x 20 x 8 bytes* plus the prefilter scores (every query against every
edge).

With ``--ranked_edges N`` there is no queries x edges score table: the prefilter scores the queries
in batches of 20,000 and each keeps only its first N ranked edges (164 covers pplacer's default
search, which uses at most 41 successful ones), as tiled placement does. A query whose search goes
past them, after many failed optimisations, is scored against every edge again and its ranking
continues in the same order, so results are identical; memory then barely grows with the number of
queries.

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

.. figure:: /_static/figures/tiling.svg
   :width: 100%
   :alt: A 24-leaf tree cut into five coloured tiles of at most 12 edges. Each tile root except the tree root keeps a distal vector D from the bottom-up pass and a proximal vector P from the top-down pass; placement then holds one tile at a time.

   gPlace's cut (``tile_layout``) of a 24-leaf tree with ``--tile_edges 12``. With its two boundary
   vectors (D and P), each tile is rebuilt on its own, exactly as in the untiled cache.

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

How the engine runs on different CPUs, the rules that keep it exact, and what makes it fast are
described in :ref:`advanced/engine`.

.. _how_it_works/comparison:

Comparison with other tools
---------------------------

Other placement tools make large trees affordable by changing *what* is computed: they place each
genome on a part of the tree (BSCAMPP, Splicer), use distances instead of likelihoods (APPLES-2),
or use another model (EPA-ng). Their placements therefore differ from pplacer's for some genomes.
gPlace changes only what is stored, so it gives pplacer's placements.

In short:

* Only gPlace gave pplacer's placement for every genome, and it was also among the fastest: 5
  minutes against pplacer's 6 hours on tree 2.
* BSCAMPP with pplacer came closest to pplacer (99.2-99.9% of the genomes on the same edge), but
  took 1 to 8 hours and, with 32 threads, about as much memory as pplacer.
* APPLES-2 was as fast as gPlace and used the least memory, but put over a quarter of the
  genomes on another edge, a few hundred of them in another genus.
* Splicer and BSCAMPP with EPA-ng put about 15% of the genomes on another edge; Splicer put
  152-221 in another phylum.

Every tool ran on the same server (64 vCPUs, 245 GB of memory), one at a time, with 32 threads and
at most 200 GiB of memory. The genomes are GlobDB r226 genomes aligned by GTDB-Tk:

* **Tree 2**: GTDB-Tk's R226 class-level reference tree 2 (21,566 reference genomes) and the 22,061
  genomes GTDB-Tk placed on it, compared with GTDB-Tk's own pplacer output;
* **Full tree**: the full GTDB R226 bacterial tree (136,646 reference genomes) and 10,000 random
  genomes, compared with pplacer's output from a 1 TB server (pplacer: 17 h 13 min and 654.9 GiB;
  gPlace on that server: 13 min 47 s).

gPlace ran with ``--max_memory 64G``: untiled on tree 2, in 5 tiles on the full tree.

.. list-table:: Tree 2: 22,061 genomes
   :header-rows: 1
   :widths: 34 14 12 16 12 12

   * - Tool
     - Time
     - Peak memory
     - Same edge as pplacer
     - Different genus
     - Different phylum
   * - pplacer v1.1.alpha19
     - 6 h 2 min
     - 131.2 GiB
     - 22,061 (100%)
     - 0
     - 0
   * - **gPlace**
     - **5 min 2 s**
     - **42.8 GiB**
     - **22,061 (100%)**
     - **0**
     - **0**
   * - APPLES-2
     - 6 min 46 s
     - 2.6 GiB
     - 16,195 (73.4%)
     - 385
     - 5
   * - BSCAMPP + pplacer (b = 2000)
     - 1 h 29 min
     - 157.4 GiB
     - 22,009 (99.8%)
     - 10
     - 1
   * - BSCAMPP + pplacer (b = 12000, V = 15)
     - 4 h 55 min
     - 127.0 GiB
     - 22,041 (99.9%)
     - 10
     - 0
   * - BSCAMPP + EPA-ng (b = 2000) [#cpj]_
     - 14 min 16 s
     - 176.9 GiB
     - 18,951 (85.9%)
     - 236
     - 9
   * - Splicer + pplacer
     - 1 h 32 min
     - 17.8 GiB
     - 18,673 (84.6%)
     - 1,436
     - 221

.. list-table:: Full tree: 10,000 genomes
   :header-rows: 1
   :widths: 34 14 12 16 12 12

   * - Tool
     - Time
     - Peak memory
     - Same edge as pplacer
     - Different genus
     - Different phylum
   * - **gPlace**
     - **24 min 15 s**
     - **56.6 GiB**
     - **10,000 (100%)**
     - **0**
     - **0**
   * - APPLES-2
     - 23 min 8 s
     - 16.3 GiB
     - 7,080 (70.8%)
     - 216
     - 6
   * - BSCAMPP + pplacer (b = 2000)
     - 1 h 12 min
     - 156.7 GiB
     - 9,918 (99.2%)
     - 11
     - 5
   * - BSCAMPP + pplacer (b = 12000, V = 15)
     - 8 h 14 min
     - 126.9 GiB
     - 9,945 (99.5%)
     - 13
     - 4
   * - BSCAMPP + EPA-ng (b = 2000) [#cpj]_
     - 50 min 37 s
     - 173.2 GiB
     - 8,472 (84.7%)
     - 115
     - 6
   * - Splicer + pplacer
     - 1 h 13 min
     - 28.9 GiB
     - 8,353 (83.5%)
     - 905
     - 152

.. [#cpj] With 8 EPA-ng jobs at a time (``--cpus-per-job 4``); with its default (16) it ran out of
   memory.

* *Same edge as pplacer*: the genome's best placement is on the same edge of the reference tree
  as pplacer's (tools that re-root the tree or re-estimate its branch lengths are compared on the
  same edges, matched by their leaves).
* *Different genus*, *Different phylum*: the GTDB genus or phylum that contains the placement is
  not the one that contains pplacer's. This is not a full GTDB-Tk classification, which also uses
  RED and ANI, but it counts the genomes placed in another named group.
* The tools: pplacer v1.1.alpha19 as run by GTDB-Tk (``-m WAG``); BSCAMPP 1.0.8, which places
  each genome on a subtree of *b* reference genomes around its closest ones (*b* = 2000 is its
  default; *b* = 12000 with *V* = 15 votes are GTDB-Tk's settings); Splicer (git 224b604), which
  splits the tree into pieces and places each genome in one of them; APPLES-2 2.0.11, which places
  by distances. With pplacer, BSCAMPP and Splicer used pplacer v1.1.alpha19 and GTDB-Tk's model
  (WAG + CAT); BSCAMPP with EPA-ng used EPA-ng's model (below).
* EPA-ng 0.3.8 is not in the tables. It has no CAT model, so it used WAG with 4 gamma rates, and
  with its default settings it ran out of memory while loading tree 2. Without pre-masking it
  placed 1,000 random tree-2 genomes in 14 min 50 s with 192.5 GiB: 892 (89.2%) on pplacer's edge,
  9 in another genus, none in another phylum (gPlace: 1,000 on pplacer's edge).
* UShER was not tested: it places nucleotide differences from one reference genome, not protein
  alignments.
