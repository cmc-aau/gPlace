******
gPlace
******

**gPlace** places genomes on the GTDB reference trees and gives **exactly the same result as
pplacer**, the placement program used by GTDB-Tk, while needing far less memory.

.. note::

   **Work in progress.** gPlace is still being tested. Before you use its results in published
   work, compare them with pplacer's on some of your own genomes (:ref:`commands/compare`), and
   please report problems on the `issue tracker <https://github.com/cmc-aau/gPlace/issues>`_.


Why gPlace?
===========

* **Same results as pplacer.** Every placement, likelihood and branch length is identical to
  pplacer v1.1.alpha19, the version GTDB-Tk uses.
* **Much less memory.** The largest GTDB-Tk class-level tree needs under 44 GiB (pplacer: about
  104 GB), and the full GTDB bacterial tree fits in 64 GB (pplacer: about 615 GiB).
* **Faster.** On the same 64-core server, gPlace placed 10,000 genomes on the full GTDB R226
  bacterial tree in 14 minutes; pplacer took 17 hours.
* **Checks itself.** ``gplace check_install`` confirms that gPlace computes exactly as pplacer
  does on your computer.


Quick start
===========

.. code-block:: bash

    # install (see Installing gPlace)
    git clone https://github.com/cmc-aau/gPlace.git
    cd gPlace
    pip install .
    gplace check_install

    # place the genomes of a GTDB-Tk alignment on a GTDB-Tk reference package
    gplace place --msa_file user_msa_file.fasta --refpkg gtdbtk.package.2.refpkg \
                 --out_file gplace.jplace --cpus 40

The :ref:`quickstart` page shows the other common tasks: the full GTDB tree, many alignments at
once, and comparing with pplacer.


Citing
======

gPlace reproduces pplacer's placements on the GTDB reference trees. If you use gPlace, please
cite **pplacer** and **GTDB** (see :ref:`references`).

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


.. toctree::
   :caption: About
   :maxdepth: 1

   how_it_works
   validation
   changelog
   references
