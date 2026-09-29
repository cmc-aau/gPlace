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
* Prebuilt ``manylinux2014`` wheel (one wheel for every Linux distribution and Python >= 3.12).
* jplace output with pplacer's number formatting (``%.12g``) and a byte-identical tree string.
