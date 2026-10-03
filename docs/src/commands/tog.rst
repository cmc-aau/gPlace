.. _commands/tog:

tog
===

The ``tog`` command makes a tree with every placed query added as a pendant edge at its best placement, as ``guppy tog`` (pplacer v1.1.alpha19) does, so that no pplacer installation is needed. The output is the same file, byte for byte: Newick (``.tog.tre``), or phyloXML with ``--xml`` (``.tog.xml``), with guppy's options and output names (``-o``, ``--out-dir``, ``--prefix``, one tree per placefile, ``placefile:groups.csv`` to split one jplace into one tree per group).

What guppy does, and so what ``tog`` does:

* A query becomes a leaf named after it (branch length = its pendant length), attached at its distal length on the edge of its best placement; a query written as one record with several names becomes a node with one zero-length leaf per name.
* The best placement is the one with the highest like-weight ratio (the first in pplacer's order), also with ``--pp``: in pplacer's code the comparison never looks past the first placement, so ``--pp`` only requires posterior probabilities to be present.
* ``--node-numbers`` writes the reference tree's node numbers where bootstraps go (``[n]``) and prefixes its labels with ``@``.

Arguments
---------

.. argparse::
   :module: gplace.cli
   :func: get_main_parser
   :prog: gplace
   :path: tog
   :nodefaultconst:

Example
-------

.. code-block:: bash

    gplace tog -o gtdbtk.bac120.classify.tree pplacer.bac120.json      # as GTDB-Tk runs guppy tog
    gplace tog --xml --out-dir trees/ placements_*.jplace

Validation
----------

Compared with the stock guppy on 607 cases, every output file identical: GTDB-Tk's 112 class-level
jplaces and gPlace's test jplaces with each option set, a 4,040-query jplace on the full GTDB
bacterial tree, every combination of ``-o`` / ``--out-dir`` / ``--prefix``, several placefiles,
gzipped input, ``placefile:groups.csv`` splitting, and synthetic placements (labels that need quoting
or XML escaping, ties, placements at the ends of an edge, multi-name records, float-like labels,
posterior probabilities).
