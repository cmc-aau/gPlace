.. _advanced/memory:

Large trees and memory limits
=============================

In short: add ``--max_memory`` (e.g. ``64G``) and gPlace stays within it, with the same result.
This page explains how, and the other memory options.

By default gPlace holds the likelihood vectors of the whole reference tree in memory, which is
the fastest way to place. For trees whose vectors do not fit, such as the full GTDB bacterial
tree (about 490 GiB untiled for 127,343 genomes), gPlace can place **tile by tile**: the tree is
cut into connected tiles and only one tile is held in memory at a time. The results are
identical to untiled placement (see :ref:`how_it_works`); tiled runs take about 1.7-2.5x longer.

Give at most one of:

* ``--max_memory 64G``: the memory gPlace may use (``K``, ``M``, ``G``, ``T`` are binary units,
  as in SLURM; a plain number is GiB). gPlace estimates the peak memory from the reference tree
  and the number of queries, keeping 5% of the limit as headroom. If untiled placement fits,
  it places without tiles; otherwise it uses the largest tiles that fit, because larger tiles
  are faster. If even tiles of 1,000 edges do not fit, it stops and prints the limit needed.
* ``--tile_edges 20000``: tiles of at most this many edges. Memory then depends on the tile
  size rather than the tree size: about 8 GiB of likelihood vectors per 10,000 edges (5,036
  sites), plus the alignments and the queries' working data.

With :ref:`commands/place_batch`, untiled or tiled placement, and the tile size, are chosen
separately for each reference package, from its tree and the number of queries pooled on it.


Estimating memory before a run
------------------------------

``gplace benchmark --refpkg <refpkg> --queries <N> --max_memory 64G`` prints the tile size and
peak memory gPlace would use, without placing anything (see :ref:`commands/benchmark`). Without
``--max_memory`` or ``--tile_edges``, gPlace warns when the estimated peak is more than the
memory available to it (the machine's memory, or a SLURM or cgroup limit).

The estimate is about *2 x internal nodes x sites x 168 bytes + edges x queries x 8 bytes*, plus
a little per worker: 43.8 GiB for GTDB-Tk R226 class-level tree 2 with 22,061 queries.


Many queries in one run
-----------------------

Untiled placement keeps every query's prefilter score for every edge (8 bytes per edge and query),
which matters with very many queries in one run. ``--ranked_edges 164`` keeps only each query's first
164 ranked edges instead (enough for pplacer's default search), with identical results: all 127,343
GTDB-Tk genomes on class-level tree 8 peak at 17.5 GiB instead of 29.5 GiB, about 7% slower.
