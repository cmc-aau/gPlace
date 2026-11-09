******
gPlace
******

**gPlace** is a replacement for **pplacer**: it places sequences on a fixed reference tree and
gives **exactly the same result as pplacer**, while needing far less memory and time.

.. note::

   **Work in progress.** gPlace is still being tested. Before you use its results in published
   work, compare them with pplacer's on some of your own genomes (:ref:`commands/compare`), and
   please report problems on the `issue tracker <https://github.com/cmc-aau/gPlace/issues>`_.


Why gPlace?
===========

* **Same results as pplacer.** Every placement, likelihood and branch length is identical to
  pplacer v1.1.alpha19, for every model pplacer supports.
* **pplacer's command line.** Replace ``pplacer`` with ``gplace place``: the same reference
  packages, query alignments and options, and the same jplace file.
* **Much less memory.** The full GTDB R226 bacterial tree fits in 64 GB (pplacer: about
  655 GiB), and the largest GTDB class-level tree needs under 44 GiB (pplacer: about 104 GB).
* **Faster.** On the same 64-core server, gPlace placed 10,000 genomes on the full GTDB R226
  bacterial tree in 14 minutes; pplacer took 17 hours.
* **Checks itself.** ``gplace check_install`` confirms that gPlace computes exactly as pplacer
  does on your computer.

Quick start
===========

.. code-block:: bash

    # install with conda (see Installing gPlace for pip)
    conda create -n gplace -c conda-forge -c bioconda gtdb-gplace
    conda activate gplace
    gplace check_install

    # place aligned sequences on a pplacer reference package: pplacer's command line
    gplace place -j 40 -c my.refpkg -o queries.jplace queries.fasta

The :ref:`quickstart` page shows the other common tasks: very large trees, many alignments at
once, and comparing with pplacer.


Citing
======

gPlace reproduces pplacer's placements. If you use gPlace, please cite **pplacer**, and the
source of your reference tree, for example **GTDB** (see :ref:`references`).

gPlace is open source, released under the
`GNU General Public License (Version 3) <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.
Questions and problems: the `issue tracker <https://github.com/cmc-aau/gPlace/issues>`_.


.. toctree::
   :caption: Getting started
   :maxdepth: 1

   installing/index
   quickstart
   faq
   announcements


.. toctree::
   :caption: Using gPlace
   :maxdepth: 1

   commands/index
   files/index


.. toctree::
   :caption: Advanced
   :maxdepth: 1

   advanced/memory
   advanced/pplacer_options
   advanced/own_trees
   models
   advanced/engine


.. toctree::
   :caption: About
   :maxdepth: 1

   how_it_works
   validation
   changelog
   references
