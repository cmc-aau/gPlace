.. _commands/place:

place
=====

The ``place`` command places the aligned queries of one FASTA file on one pplacer reference package and writes a jplace file. It computes exactly what GTDB-Tk's ``pplacer -m WAG -j <cpus> -c <refpkg> -o <out.json> <msa>`` computes: same placements, same printed numbers, same tree string. The model comes from the reference package (WAG, LG, JTT or GTR; CAT or gamma; see :ref:`models`); ``--fix_gamma_underflow`` avoids pplacer's gamma underflow on large trees (:ref:`models/gamma_underflow`).

Arguments
---------

.. argparse::
   :module: gplace.cli
   :func: get_main_parser
   :prog: gplace
   :path: place
   :nodefaultconst:

Memory limits and large trees
-----------------------------

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

``gplace benchmark --refpkg <refpkg> --queries <N> --max_memory 64G`` prints the tile size and
peak memory gPlace would use, without placing anything. Without either option, gPlace warns
when the estimated peak is more than the memory available to it (the machine's memory, or a
SLURM or cgroup limit).

Files output
------------

* :ref:`jplace file <files/jplace>` (``--out_file``)
* :ref:`gplace.log <files/gplace.log>` (next to the jplace file)

Example
-------

Input
^^^^^

.. code-block:: bash

    gplace place --msa_file user_msa_file.fasta --refpkg mini.refpkg --out_file gplace.jplace --cpus 4

For a GTDB-Tk class-level tree, pass e.g. ``--refpkg release226/split/class_level/pplacer/gtdbtk.package.2.refpkg``.
On the full GTDB bacterial tree within a 64 GB job:

.. code-block:: bash

    gplace place --msa_file gtdbtk.bac120.user_msa.fasta --refpkg release226/pplacer/gtdb_r226_bac120.refpkg \
        --out_file gplace.jplace --cpus 40 --max_memory 64G

Output
^^^^^^

.. code-block:: text

    [2026-10-01 19:23:40] INFO: gPlace v0.1.0
    [2026-10-01 19:23:40] INFO: gplace place --msa_file user_msa_file.fasta --refpkg mini.refpkg --out_file gplace.jplace --cpus 4
    [2026-10-01 19:23:42] INFO: Read 20 sequences from user_msa_file.fasta
    [2026-10-01 19:23:44] INFO: Model: WAG + CAT (20 rates), model frequencies (protein)
    [2026-10-01 19:23:44] INFO: Placing 20 queries on mini.refpkg (493 nodes, 4,891 sites, 4 CPUs; estimated peak memory 0.9 GiB).
    [2026-10-01 19:23:47] INFO: Cached likelihood vectors for 246 internal nodes (0.4 GB).
    [2026-10-01 19:23:48] INFO: Scored 20 queries against 492 edges.
    [2026-10-01 19:23:51] INFO: Wrote 20 placements to gplace.jplace (8.0 s).
