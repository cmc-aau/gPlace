Announcements
=============

gPlace 0.1.2
------------

*November, 2026*

* Peak memory is the same from run to run. Before, identical runs could peak at different levels
  (for example 13.6 or 14.2 GiB on CheckM's reference tree with 64 CPUs).
* ``place --full_model`` gives the same gamma rates, and so the same placements, on every machine.
* Placements with pplacer's models are unchanged.

gPlace 0.1.1
------------

*November, 2026*

* gPlace is on Bioconda and PyPI as ``gtdb-gplace``:
  ``conda create -n gplace -c conda-forge -c bioconda gtdb-gplace`` or
  ``python -m pip install gtdb-gplace``. Bioconda's packages are prebuilt for Linux and macOS;
  pip installs a prebuilt package on Linux and compiles the engine on macOS. Run
  ``gplace check_install`` afterwards (see :ref:`installing`).
* Tiled placement (``--max_memory``, ``--tile_edges``) is much faster on trees cut into many tiles,
  with identical results.
* Large reference packages load with much less memory (1.3 GiB instead of 4.8 GiB for the GTDB R232
  bacterial tree), and a ``--max_memory`` that is too small stops with a message before loading
  instead of being killed by a hard memory limit.

Work in progress
----------------

*October, 2026*

* gPlace is under active development and its results are not yet fully tested beyond the data
  sets described in :ref:`validation`. Check it on your own data before relying on it for
  published results (for example with :ref:`commands/compare`), and report problems on the
  `issue tracker <https://github.com/cmc-aau/gPlace/issues>`_.
