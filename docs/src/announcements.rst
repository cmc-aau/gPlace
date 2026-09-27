Announcements
=============

gPlace 0.1.0
------------

*September, 2026*

* First release. gPlace reproduces GTDB-Tk's pplacer placements (v1.1.alpha19) bit for bit on
  the WAG + CAT reference packages, with under 44 GiB of memory for the largest GTDB-Tk R226
  class-level tree.
* The C engine is compiled automatically by ``pip install``; verify it with
  ``gplace check_install`` and ``gplace test``.
