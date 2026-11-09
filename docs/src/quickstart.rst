.. _quickstart:

Quick start
===========

gPlace takes the same two inputs as pplacer:

* **a reference package** (a ``.refpkg`` folder), as pplacer reads with ``-c``: the reference
  tree, its alignment and its model, for example one of the GTDB reference packages or one built
  with :ref:`gplace taxit create <commands/taxit>`;
* **an alignment of your sequences** (FASTA or Stockholm), aligned to the reference alignment, for
  example ``queries.fasta``.

It writes a :ref:`jplace file <files/jplace>`, the same file pplacer writes.


1. Check the installation
-------------------------

.. code-block:: bash

    gplace check_install     # gPlace computes exactly as pplacer on this computer
    gplace test              # places 20 bundled sequences and compares them with pplacer's output

Both should end with a success message. See :ref:`installing` if gPlace is not installed yet.


2. Place one alignment
----------------------

Take a pplacer command and replace ``pplacer`` with ``gplace place``:

.. code-block:: bash

    gplace place -j 40 -c my.refpkg -o queries.jplace queries.fasta

This gives the same result as ``pplacer -j 40 -c my.refpkg -o queries.jplace queries.fasta``.
gPlace's own option names work too:

.. code-block:: bash

    gplace place --refpkg my.refpkg --msa_file queries.fasta --out_file queries.jplace --cpus 40

Not sure how many CPUs to use? ``--cpus auto`` (or ``-j auto``) picks a good number for your computer.


3. Stay within a memory limit
-----------------------------

Large trees, such as the full GTDB bacterial tree, need more memory than most computers have.
Give gPlace a limit and it splits the tree into pieces when it has to; the result is the same.

.. code-block:: bash

    gplace place --msa_file queries.fasta --refpkg gtdb_r226_bac120.refpkg \
                 --out_file queries.jplace --cpus 40 --max_memory 64G

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

    gplace compare --jplace_a queries.jplace --jplace_b pplacer.jplace

.. code-block:: text

    INFO: Genomes: 20 (only in A: 0, only in B: 0)
    INFO: Same reference tree string: True
    INFO: Best edge identical:   20 (100.0000%)
    INFO: Edge list identical:   20 (100.0000%)
    INFO: All values identical:  20 (100.0000%)


6. See the placements on a tree
-------------------------------

.. code-block:: bash

    gplace tog queries.jplace    # writes queries.tog.tre: the reference tree with your sequences added

This is the same tree ``guppy tog`` makes, and can be opened in any tree viewer.


Using your own tree
-------------------

gPlace also places on trees you built yourself, with FastTree, IQ-TREE, RAxML or PhyML: see
:ref:`advanced/own_trees`.
