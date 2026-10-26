.. _commands/place:

place
=====

The ``place`` command places the genomes of one alignment file on one reference package and
writes a :ref:`jplace file <files/jplace>`. It computes exactly what GTDB-Tk's
``pplacer -m WAG -j <cpus> -c <refpkg> -o <out.json> <msa>`` computes: the same placements and the
same numbers.

.. code-block:: bash

    gplace place --msa_file user_msa_file.fasta --refpkg gtdbtk.package.2.refpkg \
                 --out_file gplace.jplace --cpus 40


The options you will usually need
---------------------------------

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Option
     - What it does
   * - ``--msa_file``
     - The alignment of your genomes (FASTA), e.g. GTDB-Tk's ``user_msa_file.fasta``.
   * - ``--refpkg``
     - The reference package (a ``.refpkg`` folder), e.g. ``gtdbtk.package.2.refpkg``.
   * - ``--out_file``
     - The jplace file to write.
   * - ``--cpus``
     - Number of worker processes, or ``auto`` to pick a good number for this computer.
   * - ``--max_memory``
     - Memory gPlace may use, e.g. ``64G``. Needed for large trees such as the full GTDB
       bacterial tree; the result is the same. See :ref:`advanced/memory`.


More options
------------

Everything else is for special cases, and each option keeps pplacer's default:

* :ref:`advanced/memory`: tile size, memory estimates, many queries in one run.
* :ref:`advanced/pplacer_options`: pplacer's search settings, pre-masking, Stockholm query files,
  taxonomic classification, and reference files instead of a reference package (pplacer's
  ``-t``/``-r``/``-s``).
* :ref:`advanced/own_trees`: trees fitted with IQ-TREE or RAxML-NG, including ``--full_model``.
* :ref:`models`: the evolutionary models, and ``--fix_gamma_underflow``.

A pplacer command line can be reused as it is: gPlace also accepts pplacer's spelling of every
option (``-c`` for ``--refpkg``, ``-o`` for ``--out_file``, ``--keep-at-most`` for
``--keep_at_most``, ...).


Files output
------------

* :ref:`jplace file <files/jplace>` (``--out_file``)
* :ref:`gplace.log <files/gplace.log>` (next to the jplace file)


Example
-------

.. code-block:: bash

    gplace place --msa_file user_msa_file.fasta --refpkg mini.refpkg --out_file gplace.jplace --cpus 4

.. code-block:: text

    [2026-10-01 19:23:40] INFO: gPlace v0.1.0
    [2026-10-01 19:23:40] INFO: gplace place --msa_file user_msa_file.fasta --refpkg mini.refpkg --out_file gplace.jplace --cpus 4
    [2026-10-01 19:23:42] INFO: Read 20 sequences from user_msa_file.fasta
    [2026-10-01 19:23:44] INFO: Model: WAG + CAT (20 rates), model frequencies (protein)
    [2026-10-01 19:23:44] INFO: Placing 20 queries on mini.refpkg (493 nodes, 4,891 sites, 4 CPUs; estimated peak memory 0.9 GiB).
    [2026-10-01 19:23:47] INFO: Cached likelihood vectors for 246 internal nodes (0.4 GB).
    [2026-10-01 19:23:48] INFO: Scored 20 queries against 492 edges.
    [2026-10-01 19:23:51] INFO: Wrote 20 placements to gplace.jplace (8.0 s).

On the full GTDB bacterial tree within a 64 GB job:

.. code-block:: bash

    gplace place --msa_file gtdbtk.bac120.user_msa.fasta --refpkg release226/pplacer/gtdb_r226_bac120.refpkg \
        --out_file gplace.jplace --cpus 40 --max_memory 64G


All options
-----------

.. argparse::
   :module: gplace.cli
   :func: get_main_parser
   :prog: gplace
   :path: place
   :nodefaultconst:
