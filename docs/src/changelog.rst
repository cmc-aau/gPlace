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
* Runs where pplacer runs on Linux: x86-64 (AVX2 kernels with AVX2 + FMA, otherwise a portable
  engine with identical results) and ARM64 (NEON kernels).
* pplacer's taxonomic ``classification`` field for reference packages with a taxonomy (edge
  painting, or ``--mrca_class``); identical query sequences placed once and written as one record;
  queries that do not overlap the reference alignment stop the run, or are dropped with
  ``--discard_nonoverlapped``.
* Prebuilt ``manylinux2014`` wheels for x86-64 and ARM64 (one wheel per architecture for every
  Linux distribution and Python >= 3.12).
* jplace output with pplacer's number formatting (``%.12g``) and a byte-identical tree string.
* The engine is written in Rust (``gplace/engine/rust``, no dependencies; installing from source
  needs Rust 1.75 or newer); GSL 1.16, CORE-MATH and glibc 2.19's ``hypot`` and ``log1p`` stay C.
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
* Gamma models with 21 or more categories (GSL calls ``log1p``: glibc 2.19's, as in pplacer) are
  supported; a GSL error while the gamma rates are computed (e.g. alpha 10^6) stops gPlace with
  pplacer's message, as pplacer stops.
