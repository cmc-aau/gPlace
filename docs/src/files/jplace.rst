.. _files/jplace:

jplace file
===========

Placements in the jplace format (version 3), as written by pplacer, e.g. GTDB-Tk's
``pplacer.class_level.bac120.json``.

* ``tree``: the reference tree with ``{edge}`` numbers, byte-identical to pplacer's (labels
  quoted only when needed, branch lengths ``%g``).
* ``placements``: one record per query, ``{"p": [[...], ...], "nm": [[name, 1]]}``, placements
  sorted by like-weight ratio.
* ``fields``: ``distal_length, edge_num, like_weight_ratio, likelihood, pendant_length``.

Numbers are printed as pplacer prints them (``%.12g``), so with gPlace's bit-identical results
every number is textually identical to pplacer's. Whitespace, key order and the order of the
query records are JSON-insignificant and may differ. As in pplacer, a non-finite value would be
written as ``nan``/``inf``.

Produced by
-----------

* :ref:`commands/place`
* :ref:`commands/place_batch`
* :ref:`commands/test`

Example
-------

.. code-block:: text

    {"tree": "(GB_GCA_016183025.1:0.31573{0},GB_GCA_036381125.1:0.34958{1},...)d__Bacteria{13713};", "placements":
      [
      {"p": [[0.0758007809161, 11051, 1, -13036489.6202, 0.0629997563924]], "nm": [["SPIREOTU_01109497", 1]]},
      ...
      ], "metadata": {"invocation": "gplace place ..."}, "version": 3, "fields":
      ["distal_length", "edge_num", "like_weight_ratio", "likelihood", "pendant_length"]
    }
