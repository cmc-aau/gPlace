.. _commands/place:

place
=====

The ``place`` command places the aligned queries of one FASTA file on one pplacer reference package and writes a jplace file. It computes exactly what GTDB-Tk's ``pplacer -m WAG -j <cpus> -c <refpkg> -o <out.json> <msa>`` computes: same placements, same printed numbers, same tree string.

Arguments
---------

.. argparse::
   :module: gplace.cli
   :func: get_main_parser
   :prog: gplace
   :path: place
   :nodefaultconst:

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

Output
^^^^^^

.. code-block:: text

    [2026-09-26 21:05:55] INFO: gPlace v0.1.0
    [2026-09-26 21:05:55] INFO: gplace place --msa_file user_msa_file.fasta --refpkg mini.refpkg --out_file gplace.jplace --cpus 4
    [2026-09-26 21:05:55] INFO: Read 20 sequences from user_msa_file.fasta
    [2026-09-26 21:05:55] INFO: Placing 20 queries on mini.refpkg (493 nodes, 4,891 sites, 4 CPUs).
    [2026-09-26 21:05:56] INFO: Cached likelihood vectors for 246 internal nodes (0.4 GB).
    [2026-09-26 21:05:56] INFO: Scored 20 queries against 492 edges.
    [2026-09-26 21:05:58] INFO: Wrote 20 placements to gplace.jplace (3.2 s).
