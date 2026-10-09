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
  64-bit Linux on x86-64 or ARM64 and macOS on Apple Silicon or Intel, the platforms pplacer runs on.
  x86-64 CPUs with AVX2 and FMA use AVX2 kernels, ARM64 CPUs (including Apple Silicon) NEON kernels, and other x86-64 CPUs
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
  worker): 43.8 GiB peak for GTDB-Tk R226 class-level tree 2 with 22,061 queries. With many queries
  in one run, ``--ranked_edges 164`` keeps only each query's first ranked edges instead of its score
  for every edge (identical results): all 127,343 GTDB-Tk genomes on class-level tree 8 then peak
  at 17.5 GiB instead of 29.5 GiB. To stay within
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
  Protein and DNA alignments with every model pplacer v1.1.alpha19 has (WAG, LG, JTT, GTR; CAT or
  gamma rates; model or empirical frequencies), with pplacer's results; see :ref:`models`. With
  ``--full_model``, also IQ-TREE's and RAxML-NG's other models (next question). Other data types (codons,
  binary, morphology) are not supported.

**Which IQ-TREE and RAxML-NG models can gPlace use?**
  By default gPlace places as pplacer does: a reference package built with ``taxit create`` from an
  IQ-TREE or RAxML-NG fit, or the ``.iqtree`` file or RAxML-NG log given with ``--stats``, gives the
  substitution matrix (WAG, LG, JTT or GTR), the gamma model and GTR's rates; the rest of the fit is not
  used, and gPlace says so (see :ref:`commands/place`). ``--full_model`` places with the fit's complete
  model instead (validated against the fitting program's likelihoods; results then differ from pplacer's):

  * IQ-TREE: every reversible DNA model and protein matrix of IQ-TREE 3 (including the Q matrices), the
    fitted frequencies, ``+G`` (mean or median), ``+I``, ``+I+G``, ``+R``, ``+I+R``, and the protein mixtures
    C10 to C60, EX2,
    EX3, EHO, UL2, UL3, LG4M, LG4X and ``MIX{...}``;
  * RAxML-NG: every DNA model and protein matrix of RAxML-NG 2 (including the Q matrices and ``PROTGTR``),
    ``+FC``, ``+FE``, ``+FU`` (and ``+FO`` fits), ``+G`` with mean or median rates, ``+I``, ``+IC``, ``+R``, and
    the mixtures LG4M and LG4X.

  Not supported yet, with ``--full_model``:

  * invariable sites inside a mixture (e.g. ``LG+C20+I+G4``, ``LG4M+I``);
  * partitioned models (IQ-TREE ``-p``, ``-q``, ``-Q``; RAxML-NG partition files): linked or proportional
    branch lengths are feasible (each site would carry every partition's rate categories, so cost grows with
    their total); unlinked branch lengths have no single placement position;
  * IQ-TREE's fitted protein exchangeabilities (GTR20), user-defined matrices (``-mdef``, model files) and
    DNA mixtures; RAxML-NG's custom protein symmetries (``PROT...``) and frequencies of 1e-6 or less;
  * models outside reversible DNA and protein substitution: codon, binary, multi-state, genotype and
    morphology models, non-reversible models (Lie Markov, UNREST, NONREV, NQ matrices), heterotachy (``+H``,
    GHOST), PoMo, ascertainment bias (``+ASC``) and site frequency models (PMSF, ``--site-freq``).

**Does the output file look exactly like pplacer's?**
  Every value and the tree string are identical. JSON whitespace, key order and the order of the
  query records can differ, which JSON readers (including GTDB-Tk) ignore.
