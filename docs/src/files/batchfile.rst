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

    runs/sample1.fasta	refpkgs/tree_a.refpkg	out/sample1.tree_a.jplace
    runs/sample2.fasta	refpkgs/tree_a.refpkg	out/sample2.tree_a.jplace
    runs/sample3.fasta	refpkgs/tree_b.refpkg	out/sample3.tree_b.jplace
