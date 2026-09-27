.. _commands/place_batch:

place_batch
===========

The ``place_batch`` command places the alignments listed in a :ref:`batchfile <files/batchfile>`. All alignments that share a reference package are pooled and placed in one pass, so the reference likelihood vectors are computed once per reference package. A query's placement does not depend on which other queries are placed with it, so pooling gives exactly the same results as placing each file separately.

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
