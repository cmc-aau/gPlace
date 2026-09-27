.. _commands/test:

test
====

The ``test`` command places 20 GlobDB genomes on a small reference package bundled with gPlace (one 247-genome clade of the GTDB-Tk R226 class-level reference package 8) and compares the result with the output of GTDB-Tk's pplacer v1.1.alpha19 on the same input. It passes only if every printed value and the tree string are identical.

Arguments
---------

.. argparse::
   :module: gplace.cli
   :func: get_main_parser
   :prog: gplace
   :path: test
   :nodefaultconst:

Files output
------------

* ``test/gplace_test.jplace``
* ``test/gplace_test_report.tsv``

Example
-------

.. code-block:: bash

    gplace test --out_dir /tmp/gplace_test --cpus 8

.. code-block:: text

    [2026-09-26 20:20:13] INFO: Placing 20 queries on gplace/tests/data/mini.refpkg (493 nodes, 4,891 sites, 8 CPUs).
    [2026-09-26 20:20:16] INFO: 20/20 genomes: every printed value identical to GTDB-Tk's pplacer v1.1.alpha19; reference tree string identical: True.
    [2026-09-26 20:20:16] INFO: Test passed (2.2 s). Output: /tmp/gplace_test
