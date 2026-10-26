.. _commands:

####################
  Commands
####################

Every command prints its options with ``-h``, e.g. ``gplace place -h``.

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Command
     - What it is for
   * - **Placing genomes**
     -
   * - :ref:`place <commands/place>`
     - Place the genomes of one alignment on one reference package.
   * - :ref:`place_batch <commands/place_batch>`
     - Place many alignments at once, each reference package's alignments together.
   * - :ref:`compare <commands/compare>`
     - Compare two jplace files, e.g. gPlace's with pplacer's.
   * - **pplacer tools, without installing pplacer**
     -
   * - :ref:`tog <commands/tog>`
     - Make a tree with your genomes added at their placements (``guppy tog``).
   * - :ref:`taxit <commands/taxit>`
     - Build or check a reference package (taxtastic's ``taxit create``, ``check``, ``info``).
   * - **Setting up and checking**
     -
   * - :ref:`check_install <commands/check_install>`
     - Check that gPlace computes exactly as pplacer on this computer.
   * - :ref:`test <commands/test>`
     - Place bundled genomes and compare them with pplacer's output.
   * - :ref:`benchmark <commands/benchmark>`
     - Find the fastest ``--cpus`` and estimate a job's memory and time.
   * - :ref:`build_engine <commands/build_engine>`
     - Rebuild the compiled engine without reinstalling (advanced).

.. toctree::
   :hidden:

   place
   place_batch
   compare
   tog
   taxit
   check_install
   test
   benchmark
   build_engine (advanced) <build_engine>
