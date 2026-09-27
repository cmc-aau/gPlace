.. _files/batchfile:

batchfile
=========

Input of :ref:`commands/place_batch`: a tab-separated file with three columns and no header.

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Column
     - Content
   * - 1
     - Query alignment (FASTA, aligned to the reference package's alignment).
   * - 2
     - Reference package directory (contains ``CONTENTS.json``).
   * - 3
     - Output jplace file.

Rows sharing a reference package are placed together. Query names must be unique within each
reference package. Empty lines and lines starting with ``#`` are ignored.

Example
-------

.. code-block:: text

    runs/part00/tree_2/user_msa_file.fasta	pplacer/gtdbtk.package.2.refpkg	runs/part00/tree_2/pplacer.class_level.bac120.json
    runs/part01/tree_2/user_msa_file.fasta	pplacer/gtdbtk.package.2.refpkg	runs/part01/tree_2/pplacer.class_level.bac120.json
    runs/part00/tree_8/user_msa_file.fasta	pplacer/gtdbtk.package.8.refpkg	runs/part00/tree_8/pplacer.class_level.bac120.json
