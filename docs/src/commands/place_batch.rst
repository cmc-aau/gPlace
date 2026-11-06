.. _commands/place_batch:

place_batch
===========

The ``place_batch`` command places many alignments in one run, for example the alignments of many
samples on a few reference trees. List them in a :ref:`batchfile <files/batchfile>`,
one line per alignment: the alignment, its reference package and the output jplace file.

Alignments that share a reference package are placed together, so the work on each reference
tree is done once. The results are exactly the same as placing each file separately with
:ref:`place <commands/place>`.

The options of :ref:`place <commands/place>` work here too. ``--max_memory`` is applied to each
reference package separately (see :ref:`advanced/memory`). Reference files without a reference
package, and query files holding the reference sequences, need ``place``
(see :ref:`advanced/pplacer_options`).

Arguments
---------

.. argparse::
   :module: gplace.cli
   :func: get_main_parser
   :prog: gplace
   :path: place_batch
   :nodefaultconst:

Files output
------------

* one :ref:`jplace file <files/jplace>` per batchfile row
* :ref:`gplace.batch_summary.tsv <files/gplace.batch_summary.tsv>`
* :ref:`gplace.log <files/gplace.log>`

Example
-------

.. code-block:: bash

    gplace place_batch --batchfile batch.tsv --out_dir gplace_out --cpus 40

    # within 32 GB: tiles only for the reference packages that need them
    gplace place_batch --batchfile batch.tsv --out_dir gplace_out --cpus 40 --max_memory 32G
