Announcements
=============

gPlace 0.1.0
------------

*September, 2026*

* First release. gPlace reproduces GTDB-Tk's pplacer placements (v1.1.alpha19) bit for bit on
  the WAG + CAT reference packages, with under 44 GiB of memory for the largest GTDB-Tk R226
  class-level tree.
* The engine (Rust) is compiled automatically by ``pip install``; verify it with
  ``gplace check_install`` and ``gplace test``.
* Trees too large for memory, such as the full GTDB bacterial tree, can be placed tile by tile
  with identical results: give ``--max_memory`` (e.g. ``64G``) and gPlace picks the tile size,
  or set it with ``--tile_edges``.
