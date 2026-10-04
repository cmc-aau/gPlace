.. _commands/test:

test
====

The ``test`` command places 20 GlobDB genomes on a small reference package bundled with gPlace (one 247-genome clade of the GTDB-Tk R226 class-level reference package 8) and compares the result with the output of GTDB-Tk's pplacer v1.1.alpha19 on the same input. It passes only if every printed value and the tree string are identical.

With ``--models`` it also places the bundled test sets of pplacer's other models (see :ref:`models`)
and compares them with the output of the stock pplacer v1.1.alpha19: the same reference package
with LG, JTT, empirical frequencies and gamma (4 and 8 rates, alpha 0.05 to 1.3), pplacer's own
PhyML test data (HIV DNA with GTR + gamma, actin with WAG + gamma) and a FastTree GTR + CAT DNA
package. One more set checks ``--fix_gamma_underflow`` against gPlace's validated output.

A test set whose untiled placement would need more memory than is available (the gamma sets with
21 to 32 rate categories need 8 to 12 GiB) is placed tile by tile, with identical results, so
``gplace test --models`` also runs on an 8 GB machine.

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
