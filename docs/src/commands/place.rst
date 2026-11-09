.. _commands/place:

place
=====

The ``place`` command places aligned sequences on a reference package and writes a
:ref:`jplace file <files/jplace>`. It takes pplacer's command line and computes exactly what pplacer
computes: the same placements and the same numbers. Replace ``pplacer`` with ``gplace place``:

.. code-block:: bash

    pplacer      -j 40 -c my.refpkg -o queries.jplace queries.fasta
    gplace place -j 40 -c my.refpkg -o queries.jplace queries.fasta    # the same jplace file

or use gPlace's option names:

.. code-block:: bash

    gplace place --msa_file queries.fasta --refpkg my.refpkg \
                 --out_file queries.jplace --cpus 40

The options you will usually need
---------------------------------

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Option
     - What it does
   * - query file(s), or ``--msa_file``
     - The alignment of your sequences (FASTA or Stockholm, gzipped or not), aligned to the
       reference alignment, at the end of the command line as in pplacer.
   * - ``-c``, ``--refpkg``
     - The reference package (a ``.refpkg`` folder).
   * - ``-o``, ``--out_file``
     - The jplace file to write. Without it, as in pplacer, ``queries.fasta`` is written to
       ``queries.jplace`` (the file name without ``.gz`` and its last extension) in ``--out-dir``.
   * - ``--out-dir``, ``--out_dir``
     - The directory for those jplace files (default: the current directory).
   * - ``-j``, ``--cpus``
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

pplacer's command line works as it is, with pplacer's results:

* every option is accepted under pplacer's name and under gPlace's (``--keep-at-most`` or
  ``--keep_at_most``, ``-p``, ``-t``/``-r``/``-s``/``-d``, ``-m``, ...);
* several query files are placed one after the other, each into its own jplace file (``-o`` then
  stops the run, as in pplacer), and query files may come before or after the options;
* without a query file gPlace warns and stops with exit status 0, as pplacer does;
* ``--mmap-file``, ``--verbosity`` and ``--timing`` are accepted and change nothing (gPlace needs no
  scratch file and logs its timings in ``gplace.log``);
* the pplacer options gPlace does not implement stop the run with a message: see
  :ref:`advanced/pplacer_options`.

Two differences: gPlace creates a missing ``--out-dir`` (pplacer stops), and it reads a query file
without a ``.fasta``/``.fa``/``.sto``/``.sth`` suffix as FASTA (pplacer stops with "unfamiliar
suffix").


Files output
------------

* :ref:`jplace file <files/jplace>` (``--out_file``)
* :ref:`gplace.log <files/gplace.log>` (next to the jplace file)


Example
-------

.. code-block:: bash

    gplace place --msa_file queries.fasta --refpkg mini.refpkg --out_file gplace.jplace --cpus 4

.. code-block:: text

    [2026-10-01 19:23:40] INFO: gPlace v0.1.0
    [2026-10-01 19:23:40] INFO: gplace place --msa_file queries.fasta --refpkg mini.refpkg --out_file gplace.jplace --cpus 4
    [2026-10-01 19:23:42] INFO: Read 20 sequences from queries.fasta
    [2026-10-01 19:23:44] INFO: Model: WAG + CAT (20 rates), model frequencies (protein)
    [2026-10-01 19:23:44] INFO: Placing 20 queries on mini.refpkg (493 nodes, 4,891 sites, 4 CPUs; estimated peak memory 0.9 GiB).
    [2026-10-01 19:23:47] INFO: Cached likelihood vectors for 246 internal nodes (0.4 GB).
    [2026-10-01 19:23:48] INFO: Scored 20 queries against 492 edges.
    [2026-10-01 19:23:51] INFO: Wrote 20 placements to gplace.jplace (8.0 s).

On the full GTDB bacterial tree within a 64 GB job:

.. code-block:: bash

    gplace place --msa_file queries.fasta --refpkg gtdb_r226_bac120.refpkg \
        --out_file queries.jplace --cpus 40 --max_memory 64G


All options
-----------

.. argparse::
   :module: gplace.cli
   :func: get_main_parser
   :prog: gplace
   :path: place
   :nodefaultconst:
