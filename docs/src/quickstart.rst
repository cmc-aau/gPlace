.. _quickstart:

Quick start
===========

gPlace takes the same two inputs as pplacer does inside GTDB-Tk:

* **an alignment of your genomes** (FASTA), as written by GTDB-Tk, for example
  ``user_msa_file.fasta`` for a class-level tree or ``gtdbtk.bac120.user_msa.fasta``;
* **a reference package** (a ``.refpkg`` folder) from the GTDB-Tk reference data, for example
  ``release226/split/class_level/pplacer/gtdbtk.package.2.refpkg`` or
  ``release226/pplacer/gtdb_r226_bac120.refpkg``.

It writes a :ref:`jplace file <files/jplace>`, the same file pplacer writes.


1. Check the installation
-------------------------

.. code-block:: bash

    gplace check_install     # gPlace computes exactly as pplacer on this computer
    gplace test              # places 20 bundled genomes and compares them with pplacer's output

Both should end with a success message. See :ref:`installing` if gPlace is not installed yet.


2. Place one alignment
----------------------

.. code-block:: bash

    gplace place --msa_file user_msa_file.fasta --refpkg gtdbtk.package.2.refpkg \
                 --out_file gplace.jplace --cpus 40

This gives the same result as GTDB-Tk's
``pplacer -m WAG -j 40 -c gtdbtk.package.2.refpkg -o pplacer.json user_msa_file.fasta``.
Not sure how many CPUs to use? ``--cpus auto`` picks a good number for your computer.


3. Stay within a memory limit
-----------------------------

Large trees, such as the full GTDB bacterial tree, need more memory than most computers have.
Give gPlace a limit and it splits the tree into pieces when it has to; the result is the same.

.. code-block:: bash

    gplace place --msa_file gtdbtk.bac120.user_msa.fasta --refpkg gtdb_r226_bac120.refpkg \
                 --out_file gplace.jplace --cpus 40 --max_memory 64G

To see how much memory a job needs before running it:

.. code-block:: bash

    gplace benchmark --refpkg gtdb_r226_bac120.refpkg --queries 10000 --max_memory 64G

More in :ref:`advanced/memory`.


4. Place many alignments at once
--------------------------------

List the alignments in a tab-separated :ref:`batchfile <files/batchfile>`, one line per
alignment: the alignment, its reference package and the output file. Alignments that share a
reference package are placed together, which is faster than placing them one by one.

.. code-block:: bash

    gplace place_batch --batchfile batch.tsv --out_dir gplace_out --cpus 40


5. Compare with pplacer
-----------------------

.. code-block:: bash

    gplace compare --jplace_a gplace.jplace --jplace_b pplacer.jplace

.. code-block:: text

    INFO: Genomes: 20 (only in A: 0, only in B: 0)
    INFO: Same reference tree string: True
    INFO: Best edge identical:   20 (100.0000%)
    INFO: Edge list identical:   20 (100.0000%)
    INFO: All values identical:  20 (100.0000%)


6. See the placements on a tree
-------------------------------

.. code-block:: bash

    gplace tog gplace.jplace     # writes gplace.tog.tre: the reference tree with your genomes added

This is the same tree ``guppy tog`` makes, and can be opened in any tree viewer.


Using your own tree
-------------------

gPlace also places on trees you built yourself, with FastTree, IQ-TREE, RAxML or PhyML: see
:ref:`advanced/own_trees`.
