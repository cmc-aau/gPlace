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
* Prebuilt ``manylinux2014`` wheel (one wheel for every Linux distribution and Python >= 3.12).
* jplace output with pplacer's number formatting (``%.12g``) and a byte-identical tree string.
