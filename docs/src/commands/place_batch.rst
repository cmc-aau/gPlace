.. _commands/place_batch:

place_batch
===========

The ``place_batch`` command places the alignments listed in a :ref:`batchfile <files/batchfile>`. All alignments that share a reference package are pooled and placed in one pass, so the reference likelihood vectors are computed once per reference package. A query's placement does not depend on which other queries are placed with it, so pooling gives exactly the same results as placing each file separately.

``--max_memory`` and ``--tile_edges`` work as for :ref:`place <commands/place>`. With
``--max_memory``, untiled or tiled placement, and the tile size, are chosen separately for each
reference package, from its tree and the number of queries pooled on it.

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

    python scripts/make_gtdbtk_batchfile.py --msa_root gtdbtk_runs \
        --refpkg_root release226/split/class_level/pplacer > batch.tsv
    gplace place_batch --batchfile batch.tsv --out_dir gplace_out --cpus 40

    # within 32 GB: tiles only for the reference packages that need them
    gplace place_batch --batchfile batch.tsv --out_dir gplace_out --cpus 40 --max_memory 32G
