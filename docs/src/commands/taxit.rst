.. _commands/taxit:

taxit
=====

``taxit create``, ``taxit check`` and ``taxit info`` from taxtastic (v1.0.1), the tool that builds pplacer's reference packages, so that no taxtastic installation is needed. A refpkg comes out as taxit writes it: the same files, the same ``CONTENTS.json`` (contents, order, log, md5 sums) and the same ``phylo_model.json``, parsed from the tree program's log or stats file (RAxML, RAxML-NG, IQ-TREE, FastTree, PhyML). As with taxit, the creation date and the file name of ``phylo_model`` (a temporary file name) differ from run to run.

* ``create``: with a taxonomy, seq_info and tree, taxit reroots the tree with pplacer's ``rppr reroot`` when ``rppr`` is installed and otherwise warns and skips it; ``gplace taxit create`` does the same (``--no-reroot`` to skip it).
* ``check``: whether the refpkg is a valid pplacer input (its files, md5 sums, and that the alignment, Stockholm alignment, seq_info and tree name the same sequences), with taxit's messages.
* ``info``: number of sequences and components, ``-n`` sequence names, ``-t`` tally per taxon, ``-l`` sequence lengths (taxit 1.0.1 crashes on ``-l``; gPlace prints the table).

Arguments
---------

.. argparse::
   :module: gplace.cli
   :func: get_main_parser
   :prog: gplace
   :path: taxit
   :nodefaultconst:

Example
-------

.. code-block:: bash

    gplace taxit create -P my.refpkg -l bac120 -f msa.faa -t tree.nwk -s fasttree.log
    gplace taxit check my.refpkg
    gplace taxit info my.refpkg

Validation
----------

Compared with taxit on 54 cases: refpkgs created from real logs of every tree program (19 files:
FastTree JTT/WAG/LG/GTR, RAxML GTR/WAG/LG+F, RAxML-NG, IQ-TREE, PhyML DNA and protein with model or
empirical frequencies), GTDB-Tk's 8 class-level refpkgs rebuilt from their own files (their
``phylo_model`` identical to the one GTDB-Tk ships), every option (rerooting with and without
``rppr``, ``--no-reroot``, ``--clobber``, ``--stats-type``, errors), ``check`` on valid and 11 kinds of
broken refpkgs, and ``info``.
