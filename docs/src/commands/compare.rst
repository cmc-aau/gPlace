.. _commands/compare:

compare
=======

The ``compare`` command compares two jplace files genome by genome. It reports how many genomes have the same best edge, the same list of placement edges, and identical values in every field, plus the largest difference per field on the best placement.

Arguments
---------

.. argparse::
   :module: gplace.cli
   :func: get_main_parser
   :prog: gplace
   :path: compare
   :nodefaultconst:

Example
-------

.. code-block:: bash

    gplace compare --jplace_a gplace.jplace --jplace_b pplacer.jplace

.. code-block:: text

    [2026-09-26 21:05:59] INFO: Genomes: 20 (only in A: 0, only in B: 0)
    [2026-09-26 21:05:59] INFO: Same reference tree string: True
    [2026-09-26 21:05:59] INFO: Best edge identical:   20 (100.0000%)
    [2026-09-26 21:05:59] INFO: Edge list identical:   20 (100.0000%)
    [2026-09-26 21:05:59] INFO: All values identical:  20 (100.0000%)
    [2026-09-26 21:05:59] INFO: Largest |difference| on the best placement: likelihood 0, like_weight_ratio 0, distal_length 0, pendant_length 0
