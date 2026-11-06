FAQ
===

**Is gPlace ready to replace pplacer in my work?**
  Not without checking. gPlace is still being tested: it reproduces pplacer's placements on the
  data sets in :ref:`validation`, but it has not been used widely beyond them. Compare its output
  with pplacer's on some of your own genomes (``gplace compare``) before relying on it, and
  report problems on the GitHub issue tracker.

**Are the results really identical to pplacer's?**
  Yes. On 127,343 genomes placed on the 8 GTDB R226 class-level trees, and 10,000 genomes on the
  full GTDB R226 bacterial tree, every placement matched pplacer v1.1.alpha19's output: the same
  edges, likelihoods, weights and branch lengths, down to the last printed digit. The same holds
  for every other model and option pplacer supports (see :ref:`validation`). ``gplace test``
  repeats this check on bundled data.

**Can I use my pplacer command lines?**
  Yes: replace ``pplacer`` with ``gplace place``. gPlace accepts pplacer's options (``-c``, ``-o``,
  ``-j``, ``--out-dir``, ``-t``/``-r``/``-s``, ``-p``, ``--keep-at-most``, ...) and query files, and
  names its jplace files as pplacer does. A few rarely used pplacer options are not implemented and
  stop the run with a message (see :ref:`advanced/pplacer_options`).

**Does the output file look exactly like pplacer's?**
  Every value and the tree are identical. The spacing and the order of the genomes in the file
  can differ, which programs that read jplace files ignore.

**How much memory do I need?**
  It depends on the tree: at most 44 GiB for the GTDB class-level trees (up to 21,566 reference
  genomes). For larger trees, such as the full GTDB bacterial tree, tell gPlace how much it
  may use, e.g. ``--max_memory 64G``, and it stays within it. To get an
  estimate for your data before running, use
  ``gplace benchmark --refpkg <refpkg> --queries <number of genomes>``. See :ref:`advanced/memory`.

**Can gPlace place on the full GTDB bacterial tree?**
  Yes, with ``--max_memory``. pplacer would need about 615 GiB for the R226 bacterial tree; with
  ``--max_memory 64G`` gPlace placed all 127,343 GlobDB genomes in 2 h 19 min, and 10,000 genomes in
  14 minutes on the same server where pplacer took 17 hours.

**How many CPUs should I use?**
  Use ``--cpus auto``, or run ``gplace benchmark`` once on each type of machine. More is not
  always faster: on a 64-core server, 32-40 workers were fastest. Memory barely grows with
  ``--cpus``.

**Which computers does gPlace run on?**
  64-bit Linux (x86-64 or ARM64) and macOS (Apple Silicon or Intel), the platforms pplacer runs
  on. The results are the same on all of them; ``gplace check_install`` verifies this on yours.

**Can I place on my own tree, or other kinds of data?**
  Yes: protein and DNA alignments, with every model pplacer has (WAG, LG, JTT, GTR; see
  :ref:`models`), on trees from FastTree, RAxML, PhyML, IQ-TREE or RAxML-NG
  (see :ref:`advanced/own_trees`). Codon, binary and morphology data are not supported.


Technical questions
-------------------

**How was identity established?**
  By comparing against jplace files written by pplacer v1.1.alpha19 and against fresh runs of its
  release binary, and, where results differed, by tracing that binary under gdb (it is statically linked
  and ships debug symbols) and comparing every likelihood evaluation bit for bit. The scripts are
  in ``scripts/validation``.

**Why does gPlace build its own GSL?**
  pplacer's binary contains GSL 1.16. GSL 2.x computes the WAG eigenvectors differently in the
  last bits, which is enough to change some secondary placements.

**Which CPU features are used?**
  x86-64 CPUs with AVX2 and FMA use AVX2 kernels, ARM64 CPUs (including Apple Silicon) NEON
  kernels, and other x86-64 CPUs a portable engine with scalar kernels. All perform the same
  operations in the same order, so the results are identical. ``GPLACE_ENGINE_VARIANT=portable``
  forces the portable engine (e.g. to check it).

**Which IQ-TREE and RAxML-NG models can gPlace use?**
  By default, the part of the fit pplacer reads (substitution matrix, gamma model, GTR's rates).
  With ``--full_model``, the fit's complete model, including IQ-TREE's and RAxML-NG's other
  matrices, ``+I``, ``+R`` and protein mixtures; the placements then differ from pplacer's. The
  full lists, and what is not supported yet, are in :ref:`advanced/own_trees`.
