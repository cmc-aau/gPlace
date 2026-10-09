Change log
==========

0.1.0
-----

Major Changes:

* Initial release: ``place``, ``place_batch``, ``compare``, ``test`` and ``check_install``
  commands; the engine is compiled during ``pip install`` (``build_engine`` rebuilds it).
* Placement engine reproducing pplacer v1.1.alpha19 (GTDB-Tk's build) bit for bit for WAG +
  FastTree CAT reference packages, with internal-node-only likelihood storage, a SIMD prefilter
  and placement in forked worker processes.
* ``benchmark`` command (fastest ``--cpus``, memory and time estimates, SLURM request) and
  ``--cpus auto``.
* Tiled placement for trees too large for memory, such as the full GTDB bacterial tree, with
  results identical to untiled placement: ``--tile_edges N`` (tiles of at most N edges) or
  ``--max_memory SIZE`` (untiled if it fits, otherwise the largest tiles that fit, chosen per
  reference package from its tree and number of queries), for ``place``, ``place_batch``,
  ``test`` and ``benchmark``. Without these options gPlace warns when the estimated peak memory
  exceeds the memory available.
* ``benchmark`` times tiled jobs tile by tile, within the job's own memory: fastest ``--cpus``
  from one rebuilt tile, and a time estimate per phase (within about 10% of measured runs on the
  full GTDB bacterial tree).
* Prefilter kernel: exponent totals in 64-bit integers, so very large trees stay on the SIMD path
  (7-37x faster prefilter on the full GTDB bacterial tree; identical scores).
* Runs where pplacer runs: Linux on x86-64 (AVX2 kernels with AVX2 + FMA, otherwise a portable
  engine with identical results) and ARM64 (NEON kernels), and macOS on Apple Silicon (NEON) and
  Intel (AVX2), compiled from source with Apple's clang, with the same results.
* pplacer's taxonomic ``classification`` field for reference packages with a taxonomy (edge
  painting, or ``--mrca_class``); identical query sequences placed once and written as one record;
  queries that do not overlap the reference alignment stop the run, or are dropped with
  ``--discard_nonoverlapped``.
* Prebuilt ``manylinux2014`` wheels for x86-64 and ARM64 (one wheel per architecture for every
  Linux distribution and Python >= 3.12).
* jplace output with pplacer's number formatting (``%.12g``) and a byte-identical tree string.
* The engine is written in C (``gplace/engine/csrc/xp.c``), like the libraries it links (GSL 1.16,
  CORE-MATH, glibc 2.19's ``hypot`` and ``log1p``): installing from source needs only a C compiler and
  ``make``.
* ``tog`` (the same output as ``guppy tog``, byte for byte, Newick or phyloXML) and ``taxit create``,
  ``taxit check``, ``taxit info`` (as taxtastic's ``taxit``): a placement workflow needs neither
  pplacer nor taxtastic.
* pplacer's settings as options of ``place`` and ``place_batch``, under gPlace's and pplacer's
  names, with pplacer's results: ``--ml_tolerance``, ``--start_pend``, ``--max_pend``,
  ``--max_strikes``, ``--strike_box``, ``--max_pitches``, ``--keep_at_most``, ``--keep_factor``,
  ``--no_pre_mask`` and ``--always_refine``.
* Model refinement as in pplacer: when the query alignment holds the reference sequences, or with
  ``--always_refine``, CAT models get new site categories on the reference tree.
* ``place`` without a reference package: ``--tree``, ``--ref_aln``, ``--stats``, ``--ref_dir``
  (pplacer's ``-t``, ``-r``, ``-s``, ``-d``), with RAxML 7 info, PhyML 3.0 stats and JSON model files
  read as pplacer reads them, and pplacer's model settings ``--model_freqs``, ``--gamma_cats``,
  ``--gamma_alpha`` (``--model`` is accepted).
* Stockholm query files (``.sto``, ``.sth``) and reference packages with ``aln_sto``; gzipped
  alignments (``.gz``).
* Source installs with older C compilers (gcc 4.8 to 9, e.g. CentOS 7): the vendored CORE-MATH
  ``sin`` and ``pow`` no longer need ``__builtin_roundeven`` (gcc 10) or ``__builtin_mul_overflow``
  (gcc 5), and gcc < 5 compiles them as C11. Before, gcc < 10 built an engine that failed to load
  (``undefined symbol: __builtin_roundeven``); the build now stops on any unknown builtin and checks
  that every engine loads. Rebuilds relink the engine whenever the compiled C changes (before, a link
  cached in the build directory could be reused). A failed engine build now ends ``pip install`` with one
  clear message (editable installs used to hide it behind a setuptools ``KeyError``).
* ``--ranked_edges N`` (``place``, ``place_batch``): untiled placement without a queries x edges score
  table; queries are scored in batches of 20,000 and keep their first N ranked edges, as in tiled
  placement, so memory barely grows with the number of queries (all 127,343 GTDB-Tk genomes on
  class-level tree 8: 17.5 GiB instead of 29.5 GiB with ``--ranked_edges 164``). Results are
  identical; the default is unchanged.
* ``gplace test`` places a test set tile by tile when its untiled placement would not fit in the
  memory available (identical results); ``gplace test --models`` and the unit tests now run on
  8 GB machines such as GitHub's standard runners.
* Gamma models with 21 or more categories (GSL calls ``log1p``: glibc 2.19's, as in pplacer) are
  supported; a GSL error while the gamma rates are computed (e.g. alpha 10^6) stops gPlace with
  pplacer's message, as pplacer stops.
* Trees fitted with IQ-TREE or RAxML-NG: packages built from them with ``taxit create`` place as in
  pplacer (identical results), and the ``.iqtree`` file or RAxML-NG log can be given directly with
  ``--stats``. gPlace warns about the parts of the fitted model pplacer does not use (``+I``, ``+R``,
  ML-estimated frequencies, ...) and explains the fits pplacer cannot use.
* ``place --full_model``: places with an IQ-TREE fit's complete model at its full precision (every
  reversible DNA model and protein matrix of IQ-TREE 3, fitted frequencies, ``+G``, ``+I``, ``+R``, protein
  mixtures C10-C60, EX2/EX3/EHO, UL2/UL3, LG4M, LG4X, ``MIX{...}``). Not pplacer's model, so results differ
  from pplacer's; placement likelihoods equal IQ-TREE's own (83 fits). The default is unchanged.
* ``place --full_model`` with RAxML-NG fits: reads ``<prefix>.raxml.bestModel`` as RAxML-NG reads it (every
  DNA model and protein matrix of RAxML-NG 2, ``PROTGTR``, ``+FC``/``+FE``/``+FU``, ``+G`` with mean or median
  rates computed with RAxML-NG's gamma algorithm, ``+I``/``+IC``, ``+R``, LG4M, LG4X). Placement likelihoods
  equal RAxML-NG's own to the jplace file's precision (89 fits). ``--full_model`` also takes the
  ``.raxml.bestModel`` file, or finds it next to a ``--stats`` RAxML-NG log.
* ``place --full_model`` with IQ-TREE fits: gamma rates computed with IQ-TREE's own algorithm (Yang's, as in
  IQ-TREE's ``RateGamma``; bit-identical) instead of pplacer's, so placement likelihoods equal IQ-TREE's to
  its printed precision; median gamma rates (``--gamma-median``) are now supported. Reference packages and
  pplacer's models keep pplacer's gamma rates.
* ``place --full_model``: an invariable-sites category (``+I``) now evolves with the exact identity matrix
  (it was the eigen-decomposition's, exact only to rounding, which at variable sites added noise to tiny
  likelihoods): IQ-TREE LG+I fits now agree with IQ-TREE to 1e-7 instead of 3.4e-5. pplacer models are
  unchanged.
* Every file gPlace opens is closed explicitly (``gplace.fileio`` helpers or ``with`` blocks): no
  ``ResourceWarning: unclosed file`` in the unit-test log.
