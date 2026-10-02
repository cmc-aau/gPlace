FAQ
===

**Are the results really identical to pplacer's?**
  Yes, for the WAG + CAT reference packages used by GTDB-Tk. On 127,343 genomes (8 GTDB-Tk R226
  class-level trees) every placement matched GTDB-Tk's jplace files: edges, likelihoods,
  like-weight ratios and branch lengths, with every printed number textually identical.
  ``gplace test`` repeats this check on bundled data.

**How was identity established?**
  By comparing against GTDB-Tk's own jplace files and against fresh runs of GTDB-Tk's pplacer
  binary, and, where results differed, by tracing that binary under gdb (it is statically linked
  and ships debug symbols) and comparing every likelihood evaluation bit for bit. The scripts are
  in ``scripts/validation``.

**Why does gPlace build its own GSL?**
  pplacer's binary contains GSL 1.16. GSL 2.x computes the WAG eigenvectors differently in the
  last bits, which is enough to change some secondary placements.

**Which CPUs are supported?**
  64-bit Linux on x86-64 or ARM64, the platforms pplacer runs on (macOS is not supported yet).
  x86-64 CPUs with AVX2 and FMA use AVX2 kernels, ARM64 CPUs NEON kernels, and other x86-64 CPUs
  a portable engine with scalar kernels. All perform the same operations in the same order, so
  the results are identical. ``gplace check_install`` verifies the engine on your machine, and
  ``GPLACE_ENGINE_VARIANT=portable`` forces the portable engine (e.g. to check it).

**How many CPUs should I use?**
  Run ``gplace benchmark`` once per machine type, or use ``--cpus auto``. Placement is limited
  by memory bandwidth and shared cache: on a 64-vCPU AMD EPYC server 32-40 worker processes were
  fastest and more made every evaluation slower. Memory barely grows with ``--cpus``.

**How much memory is needed?**
  ``gplace benchmark --refpkg <refpkg> --queries <N>`` prints the estimate (about
  *2 x internal nodes x sites x 168 bytes + edges x queries x 8 bytes*, plus a little per
  worker): 43.8 GiB peak for GTDB-Tk R226 class-level tree 2 with 22,061 queries. To stay within
  a limit, pass ``--max_memory`` (e.g. ``64G``): gPlace places untiled if the estimate fits, and
  tile by tile otherwise, with identical results.

**Can gPlace place on the full GTDB bacterial tree?**
  Yes, with ``--max_memory`` or ``--tile_edges``. Untiled, the R226 bacterial tree (273,290 nodes)
  would need about 490 GiB for 127,343 genomes (pplacer: about 615 GiB). With ``--max_memory 64G``
  gPlace placed all 127,343 GlobDB genomes in 2 h 19 min (6 tiles, 59.2 GiB peak) and matched
  GTDB-Tk's full-tree classification tree for every genome it contains (edge, pendant and distal
  lengths). A lower limit gives smaller tiles: tiles of at most 14,000 edges need about 20 GiB.
  See :ref:`validation` and :ref:`how_it_works`.

**Can I place other kinds of data?**
  No. gPlace reproduces the WAG + FastTree CAT protein refpkgs used by GTDB-Tk. Other models
  (e.g. GTR, empirical frequencies, gamma rates) are rejected.

**Does the output file look exactly like pplacer's?**
  Every value and the tree string are identical. JSON whitespace, key order and the order of the
  query records can differ, which JSON readers (including GTDB-Tk) ignore.
